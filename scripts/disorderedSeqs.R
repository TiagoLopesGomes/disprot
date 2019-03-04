library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(ggplot2)
library(Rtsne)
source(file.path(projectDir, 'scripts', 'featurize.R'))

# Variables
# - disorderPred: data.frame
#     Columns: predictor and/or predictor_id, start, end, seqID
# - disorderDB: GRanges
#     Assumed to have metadata columns analogous to disorderPred
# - proteome: data.frame
#     Columns: seq, d2p2_id
# - proteinDB: AAStringSet
#     

aaOrder = c('P', 'G', 'A', 'E', 'D', 'K', 'R', 'H', 'Q', 'N', 'S', 'T', 'V', 'L', 'I', 'M', 'C', 'Y', 'F', 'W')
aaOther = c("U", "O", "B", "J", "Z", "X", "*", "-", "+", ".", "other")

setops_multi <- function(rangeList, fun = GenomicRanges::intersect) {
  # Apply set operations on multiple XRanges (e.g., IRanges or GRanges) objects
  #
  # Args
  # - rangeList: GRangesList or IRangesList
  # - fun: function. default = GenomicRanges::intersect
  # 
  # Returns: XRanges object
  
  range = rangeList[[1]]
  if (length(rangeList) > 1) {
    for (i in 2:length(rangeList)) {
      range = fun(range, rangeList[[i]])
    }
  }
  invisible(range)
}

disorderRanges.prop <- function(disorderDB, predictors = NULL, thresh = 0.75) {
  # Returns GRanges of disordered regions predicted by proportion of predictors
  # greater than some threshold.
  # 
  # Args
  # - disorderDB: GRanges
  #     Must contain "predictor" metadata column
  # - predictors: vector, character. default = NULL
  #     Predictors to consider. If NULL, considers all predictors
  # - thresh: numeric. default = 0.75
  #     Threshold proportion of predictors that indicate a letter (e.g., base or amino acid) to be disordered
  # 
  # Returns: GRanges

  if (!is.null(predictors)) {
    disorderDB = disorderDB[disorderDB$predictor %in% predictors]
  }
  
  numPredictors = length(unique(disorderDB$predictor))
  cov = coverage(disorderDB)
  ir = unlist(IRanges::slice(cov, lower = numPredictors * thresh, rangesOnly = TRUE))
  gr = GRanges(seqnames = names(ir), ranges = ir)
  names(gr) <- NULL
  invisible(gr)
}

disorderRanges.intersect <- function(disorderDB, predictors = c("VLXT", "VSL2b"), fun = GenomicRanges::intersect) {
  # Returns GRanges where for each sequence, the set operation defined by the fun argument
  # is applied to subsets of the GRanges defined by predictors.
  # 
  # Args
  # - disorderDB: GRanges
  #     Must contain "predictor" metadata column
  # - predictors: vector, character. default = c("VLXT", "VSL2b")
  #     Predictors by which to separate and then merge (via set operations) GRanges for a given sequence
  # - fun: function. default = GenomicRanges::intersect
  #     Set operation (e.g., intersect, union, setdiff) to perform on subsets (defined by predictors)
  #     of disorderDB for a given sequence
  # 
  # Returns: GRanges
  
  disorderDB = disorderDB[disorderDB$predictor %in% predictors]
  grList = split(disorderDB, seqnames(disorderDB))
  
  disorderRanges = NULL
  for (i in 1:length(grList)) {
    gr = grList[[i]]
    grList_predList = split(gr, gr$predictor)
    gr = setops_multi(grList_predList, fun)
    if (is.null(disorderRanges)) {
      disorderRanges = gr
    } else {
      disorderRanges = c(disorderRanges, gr)
    }
  }
  
  invisible(disorderRanges)
}

keepValidDisorder <- function(disorderPred, proteinDB, disorderIdCol = "d2p2_id") {
  # Remove invalid disorder entries.
  # 
  # Args
  # - disorderPred: data.frame
  #     Columns: start, end, `d2p2_id`
  # - proteinDB: XStringSet
  #     Each protein should be identified by its names attribute
  # - disorderIdCol: character. default="d2p2_id"
  #     Column in `disorderPred` with disorder protein ID
  #
  # Returns: GRanges
  
  num_unfiltered = nrow(disorderPred)
  
  # remove any predictions where the start or end position extends beyond the known length of the protein
  disorderPred = disorderPred %>%
    dplyr::filter(start > 0,
                  end > 0,
                  start <= width(proteinDB)[match(disorderPred[[disorderIdCol]], names(proteinDB))],
                  end <= width(proteinDB)[match(disorderPred[[disorderIdCol]], names(proteinDB))])
  
  print(paste("Removed", num_unfiltered - nrow(disorderPred), "invalid predictions"))
  
  idx = which(disorderPred$start - disorderPred$end > 0)
  
  if (length(idx) > 0) {
    warning("Following indices have flipped start/end")
    print(idx)
  }
  
  return(disorderPred)
}

extractDisorder <- function(disorderPred, proteome,
                            disorder_idCol = "d2p2_id", disorder_predCol = "predictor",
                            proteome_seqCol = "seq", proteome_idCol = "d2p2_id",
                            uniqueNames = TRUE) {
  # Convert data frames into Bioconductor objects
  # - disordered regions -> GRanges
  # - protein sequences -> AAStringSet
  # 
  # Args
  # - disorderPred: data.frame
  #     Disordered ranges, e.g., as provided by D2P2.
  #     Columns: start, end, `disorder_idCol`, `disorder_predCol`
  # - proteome: data.frame
  #     Proteome sequences from which to extract disordered regions.
  #     Columns: `proteome_seqCol`, `proteome_idCol`
  #       The `proteome_idCol` column must consist of unique entries.
  # - disorder_idCol, disorder_predCol, proteome_seqCol, proteome_idCol: character
  #     Relevant columns from data frames
  # - uniqueNames: logical. default = FALSE
  #     Give unique names "[D2P2 ID]_[start position in protein sequence]" to subsequences
  #     in disorderSeqs. If FALSE, each disordered region extracted from the same original
  #     proteome sequence will have the same "[D2P2 ID]" name.
  # 
  # Returns: list(GRanges, AAStringSet)
  
  proteinDB = Biostrings::AAStringSet(proteome[[proteome_seqCol]])
  names(proteinDB) = proteome[[proteome_idCol]]
  
  disorderDB = GRanges(seqnames = disorderPred[[disorder_idCol]],
                       ranges = IRanges(start = disorderPred$start,
                                        end = disorderPred$end))
  disorderDB$predictor = disorderPred[[disorder_predCol]]
  seqlengths(disorderDB) = width(proteinDB)[runValue(match(seqnames(disorderDB), names(proteinDB)))]
  
  disorderRanges = disorderRanges.prop(disorderDB)
  seqlengths(disorderRanges) = width(proteinDB)[runValue(match(seqnames(disorderRanges), names(proteinDB)))]
  # disorderRanges = disorderRanges.intersect(disorderDB)
  
  # Subset AAStringSet by GRanges - requires that names of proteinDB are unique
  disorderSeqs = proteinDB[disorderRanges]
  
  if (uniqueNames) {
    names(disorderSeqs) = paste(seqnames(disorderRanges), start(disorderRanges), sep = '_')
    if (sum(duplicated(names(disorderSeqs))) > 0) {
      warning('Duplicated sequence names given even with argument `uniqueNames = TRUE`')
    }
  }
  
  return(invisible(list(disorderRanges=disorderRanges, disorderSeqs=disorderSeqs)))
}

featurizeGL_precomputed <- function(gl, featuresDB, features = NULL, exact_match = FALSE, verbose = TRUE) {
  # Return feature matrix for input sequence identifiers by looking up precomputed features
  # from a database.
  # 
  # Args
  # - gl: vector, character
  #     sequence identifiers
  # - featuresDB: matrix
  #     Features database: rows = sequence identifiers, columns = features
  # - features: vector, character. default = NULL
  #     Features (e.g., AAC, PAAC, Geary, etc.) to extract by matching up to 
  #     match up to first underscore '_' in colnames of `featuresDB`.
  #     If NULL, all features are extracted.
  # - exact_match: logical. default = FALSE
  #     Match sequence identifiers in `gl` to rownames of `featuresDB` exactly.
  #     If FALSE, match up to first underscore '_' in rownames of `featuresDB`.
  #     
  #     Example 1: `gl` consists of D2P2 IDs
  #       Set to FALSE to extract features of all disordered regions of each D2P2 ID sequence
  #       --> nrow of returned matrix may be greater than number of entries in `gl`
  #     
  #     Example 2: `gl` consists of [D2P2 ID]_[start position]
  #       Set to TRUE to extract features of exact disordered regions as identified in `gl`
  # - verbose: logical. default = TRUE
  # 
  # Returns: matrix
  #   Feature matrix: rows = sequence identifiers, columns = features
  
  gl = unique(gl)
  
  # determine features (columns) to subset
  if (!is.null(features)) {
    featureClass = do.call(rbind, strsplit(x=colnames(featuresDB), split='_'))[,1]
    featuresToMatch = featureClass %in% features
  } else {
    featuresToMatch = 1:ncol(featuresDB)
  }
  
  # determine sequence identifiers (rows) to subset
  if (exact_match) {
    seqsToMatch = gl
  } else {
    # truncate rownames of featuresDB to everything before '_'
    rownames(featuresDB) = do.call(rbind, strsplit(x=rownames(featuresDB), split='_'))[,1]
    seqsToMatch = rownames(featuresDB) %in% gl
  }
  
  if (verbose) {
    nSeq = length(gl)
    seqsInFeaturesDB = gl %in% rownames(featuresDB)
    nSeqInFeaturesDB = sum(seqsInFeaturesDB)
    print(paste('Features available for', nSeqInFeaturesDB, 'out of', nSeq,
                'input sequence identifiers.'))
    print('Features not found for the following sequence identifiers:')
    print(gl[!seqsInFeaturesDB])
  }
  
  invisible(featuresDB[seqsToMatch, featuresToMatch])
}

featurizeGL_disorder <- function(gl, disorderPred, proteome,
                                 minLength = 0, lambda = minLength, nlag = minLength,
                                 ..., verbose = TRUE) {
  # Compute feature matrix for disordered regions of input sequence identifiers
  # 
  # Args: see featurizeGL()
  # - disorderPred: data.frame
  #     Disordered ranges, e.g., as provided by D2P2.
  #     Columns: start, end, d2p2_id
  # 
  # Returns: matrix
  #   Feature matrix: rows = sequence identifiers, columns = features

  gl = unique(gl)
  disorderPred_gl = dplyr::filter(disorderPred, d2p2_id %in% gl)
  
  if (verbose) {
    nSeq = length(gl)
    nSeqInDisorderPred = sum(gl %in% disorderPred_gl[['d2p2_id']])
    print(paste('Disorder predictions available for', nSeqInDisorderPred, 'out of', nSeq,
                'input sequence identifiers.'))
  }
  
  # make sure disorderPred_gl has 'predictor' column to satisfy extractDisorder()
  # - if a 'predictor_id' column is provided, rename that to 'predictor'
  # - if no predictor/predictor_id column is provided, create a column with values 1
  if (!'predictor' %in% names(disorderPred_gl)) {
    if ('predictor_id' %in% names(disorderPred_gl)) {
      disorderPred_gl = dplyr::rename(disorderPred, predictor = predictor_id)
    } else {
      disorderPred_gl['predictor'] = 1
    }
  }
  
  tmp = extractDisorder(disorderPred = disorderPred_gl, proteome = proteome, uniqueNames = FALSE)
  disorderRanges = tmp$disorderRanges
  disorderSeqs = tmp$disorderSeqs
  filter = width(disorderSeqs) > minLength
  names(disorderSeqs) = paste(seqnames(disorderRanges), start(disorderRanges), sep = '_')
  seqs = as.vector(disorderSeqs[filter])
  
  if (verbose) {
    stopifnot(identical(length(unique(names(tmp$disorderSeqs))), nSeqInDisorderPred))
    nLongSeqs = length(unique(names(tmp$disorderSeqs)[filter]))
    nDroppedProteins = nLongSeqs - nLongSeqs
    print(paste("Excluding", nDroppedProteins, "proteins without any disordered regions longer than", minLength, "aa."))
  }
  
  mat = featurizeSeqs(seqs, lambda=lambda, nlag=nlag, ...)
  invisible(mat)
}

featurizeGL <- function(gl, proteome,
                        minLength = 0, lambda = minLength, nlag = minLength,
                        ..., verbose = TRUE) {
  # Compute feature matrix for input sequence identifiers
  # 
  # Args
  # - gl: vector, character
  #     sequence identifiers
  # - proteome: data.frame
  #     Proteome sequences from which to extract disordered regions.
  #     Columns: seq, d2p2_id
  #       The d2p2_id column must consist of unique entries.
  # - lambda, nlag, features, format, save, cores: see featurizeSeqs()
  #     lambda, nlag default to `minLength`
  # - verbose: logical. default = TRUE
  # 
  # Returns: matrix
  #   Feature matrix: rows = sequence identifiers, columns = features
  
  gl = unique(gl)
  proteome = dplyr::filter(proteome, d2p2_id %in% gl)
  
  if (verbose) {
    nSeq = length(gl)
    nSeqInProteome = sum(gl %in% proteome[['d2p2_id']])
    print(paste('Sequences available for', nSeqInProteome, 'out of', nSeq,
                'input sequence identifiers.'))
  }
  
  filter = nchar(proteome[['seq']]) > minLength
  seqs = proteome[['seq']][filter]
  names(seqs) = proteome[['d2p2_id']][filter]
  
  if (verbose) {
    print(paste("Excluding", sum(!filter), "proteins of length <=", minLength, "aa."))
  }
  
  mat = featurizeSeqs(seqs, lambda=lambda, nlag=nlag, ...)
  invisible(mat)
}

compareGeneListsAllFeatures <- function(gls, featuresDB, seqtypes, lambdas, transformations, plotsDir=NULL, verbose=TRUE) {
  # Compare 2 gene lists using all available features databases.
  #
  # Args
  # - gls: list of vector, character
  #     Named list of vectors of sequence identifiers (e.g., d2p2_id).
  #     The names of this list will be used as classes in plotting.
  # - featuresDB: list of matrix
  #     Mapping from features database (character: [seqtype][lambda]) to features matrix (rows = proteins, columns = features)
  # - seqtypes: vector, character
  #     'full' (entire sequence); 'disorder' (only disordered regions)
  # - lambdas: vector, numeric
  #     sequence length thresholds
  # - transformations: vector, character
  #     'PCA'; 't-SNE'
  # - plotsDir: character. default=NULL
  #     Directory in which to store plots
  # 
  # Returns: NULL
  
  for (seqtype in seqtypes) {
    for (lambda in lambdas) {
      db_name = paste0(seqtype, lambda)
      
      if (verbose) {
        print(paste('Visualizing features from', db_name))
      }
      
      if (lambda == 0) {
        # only use length-independent features
        db = featuresDB[[db_name]]
      } else {
        # merge length-independent features with length-dependent features
        db_rownames = rownames(featuresDB[[db_name]])
        db0_rownames = rownames(featuresDB[[paste0(seqtype, 0)]])
        common_rownames = intersect(db_rownames, db0_rownames)
        db = rbind(featuresDB[[paste0(seqtype, 0)]][common_rownames,], featuresDB[[db_name]][common_rownames,])
      }
      
      for (transformation in transformations) {
        # get plot without showing it
        result = compareGeneLists(
          gls = gls,
          featurizeFun = featurizeGL_precomputed,
          featurizeArgs = list(featuresDB = db),
          plotFun = transformAndPlot,
          plotArgs = list(transformation = transformation, show = FALSE, alpha = 0.5))
        
        # add title and subtitle to plot
        result$plot = results$plot +
          labs(title = paste(gl_names_short, collapse=' v. '),
               subtitle = paste(paste('seqtype:', seqtype),
                                paste('minLength:', lambda),
                                paste('transformation:', transformation),
                                sep = ', '))
        
        if (!is.null(plotsDir)) {
          ggsave(filename = file.path(paste0(paste(transformation, db_name, 'human', sep = '_'), '.png')),
                 path = plotsDir,
                 plot = results$plot,
                 width = 7, height = 7, units = "in")
        }
      }
    }
  }
}

compareGeneLists <- function(gls, featurizeFun, featurizeArgs = NULL, plotFun = NULL, plotArgs = NULL) {
  # Compare featurization of gene lists.
  # 
  # Args
  # - gls: list of vector, character
  #     Named list of vectors of sequence identifiers (e.g., d2p2_id).
  #     The names of this list will be used as classes in plotting.
  # - featurizeFun: function
  #     Featurization function for each of the "gene lists" supplied in `gls`.
  #     Must take a "gene list" as its first argument.
  #     Examples: featurizeGL, featurizeGL_disorder, featurizeGL_precomputed
  # - featurizeArgs: list. default = NULL
  #     Arguments to pass to featurizeFun.
  # - plotFun: function. default = NULL
  #     Visualization function.
  #     Must take a feature matrix as its first argument.
  #     Examples: transformAndPlot
  # - plotArgs: list. default = NULL
  #     Arguments to pass to plotFun
  # 
  # Returns: list
  # - plot: ggplot or plotly
  # - mat_list: list
  #     named list of feature matrices corresponding to each "gene list"
  
  # list of feature matrices
  mat_list = sapply(gls, function(gl) do.call(featurizeFun, c(list(gl), featurizeArgs)))
  
  plot = NULL
  
  if (!is.null(plotFun)) {
    # row-bind feature matrices together, with rownames indicating the gene list from which
    # each row came from
    mat = do.call(rbind, mat_list)
    rownames(mat) = unlist(sapply(1:length(gls), function(idx) rep(names(gls)[idx], nrow(mat_list[[idx]]))))
  
    plot = do.call(plotFun, c(list(mat), plotArgs))
  }
  
  invisible(list(plot = plot, mat_list = mat_list))
}

transformAndPlot <- function(mat, transformation = 'PCA', normalize = TRUE,
                             dims = if (identical(transformation, 'PCA')) c('PC1', 'PC2') else c('V1', 'V2'),
                             show = TRUE, ...) {
  # Create PCA plot of matrix.
  # 
  # Args
  # - mat: matrix
  #     Data matrix: rows = samples, columns = features.
  #     Rownames of matrix are the sample labels.
  # - transformation: character. default = 'PCA'
  #     'PCA' or 't-SNE'
  # - normalize: logical. default = TRUE
  #     Normalize (0-mean, 1-variance) features before PCA transformation.
  # - dims: vector, character. default = c('PC1', 'PC2')
  #     Principal components ('PC1', 'PC2', ...) or t-SNE output dimensions ('V1', 'V2', ...)
  #     to plot. If `length(dims) == 3`, makes a 3D plot (requires plotly library).
  # - show: logical. default = TRUE
  #     Show plot, as opposed to only returning the plot object.
  # - ...
  #     Additional arguments to pass to <geom_function>() or plot_ly().
  #     Example: alpha = 0.5
  # 
  # Returns: ggplot or plotly
  
  stopifnot(length(dims) %in% 1:3)
  stopifnot(transformation %in% c('PCA', 't-SNE'))
  
  if (transformation == 'PCA') {
    p = prcomp(mat, center = normalize, scale. = normalize, tol = ifelse(normalize, 0, NULL))
    data = as.data.frame(p$x)
  } else {
    p = Rtsne(mat, dims = length(dims), normalize = normalize, check_duplicates = FALSE)
    data = as.data.frame(p$Y)
  }
  data[['label']] = rownames(mat)
  
  if (length(dims) == 1) {
    plot = ggplot(data = data) + 
      geom_dotplot(aes_string(x = dims[1], color = 'label'), ...)
  } else if (length(dims) == 2) {
    plot = ggplot(data = data) +
      geom_point(aes_string(x = dims[1], y = dims[2], color = 'label'), ...)
  } else {
    library(plotly)
    plot = plotly::plot_ly(x = data[[dims[1]]], y = data[[dims[2]]], z = data[[dims[3]]],
                           type = "scatter3d", mode = "markers", color = data[['label']], ...)
  }
  
  if (show) {
    print(plot)
  }
  invisible(plot)
}

plotCharFreq <- function(disorderSeqsList, thresh=0, scale=TRUE, plot=TRUE, aaClassMap=NULL, verbose=TRUE) {
  disorderSeqsList = sapply(disorderSeqsList, function(x) x[width(x) > thresh])
  labels = unlist(sapply(names(disorderSeqsList),
                         function(x) rep(x, length(disorderSeqsList[[x]]))))
  
  freqs_all = NULL
  for (i in 1:length(disorderSeqsList)) {
    freq_mat = alphabetFrequency(disorderSeqsList[[i]])
    freq_tb = as_tibble(data.frame(AA = colnames(freq_mat), counts = colSums(freq_mat)))
    freq_tb[["AA"]] = factor(freq_tb[["AA"]], levels = c(aaOrder))
    freq_tb[["label"]] = names(disorderSeqsList)[i]
    if (scale) {
      freq_tb[["counts"]] = freq_tb[["counts"]] / sum(freq_tb[["counts"]])
    }
    freqs_all = rbind(freqs_all, freq_tb)
  }
  
  if (!is.null(aaClassMap)) {
    freqs_all = dplyr::mutate(freqs_all, class = aaClassMap[as.character(AA)])
  } else {
    freqs_all = dplyr::mutate(freqs_all, class = '')
  }
  
  g = ggplot(data = freqs_all) +
    geom_bar(aes(x = AA, y = counts, fill = class), stat = "identity") +
    xlab("Amino Acid") + 
    labs(title="Distribution of amino acid frequencies",
         subtitle=paste('minimum disordered region length:', thresh)) +
    facet_grid(. ~ label)
  
  if (plot) {
    print(g)
  }
  invisible(g)
}

plotHist <- function(disorderRangesList, fun, varname, thresh=0, plot=TRUE) {
  disorderRangesList = sapply(disorderRangesList, function(x) x[width(x) > thresh])

  data = NULL
  for (i in 1:length(disorderRangesList)) {
    values = fun(disorderRangesList[[i]])
    df = as.data.frame(values)
    names(df) = varname
    df[["label"]] = names(disorderRangesList)[i]
    data = rbind(data, df)
  }
  
  g = ggplot(data = data) +
    geom_histogram(aes_string(x = varname, y = "..density..")) +
    facet_grid(. ~ label) +
    labs(subtitle=paste('minimum disordered region length:', thresh))
  
  if (plot) {
    print(g)
  }
  invisible(g)
}

calcProps <- function(disorderRanges, prop=NULL) {
  df = as.data.frame(disorderRanges) %>%
    group_by(seqnames) %>%
    summarise(sum_dr_length = sum(width), num_dr = length(width)) %>%
    mutate(seqlength = seqlengths(disorderRanges)[match(seqnames, names(seqlengths(disorderRanges)))]) %>%
    mutate(num_dr_norm = num_dr / seqlength) %>%
    mutate(prop_dr = sum_dr_length / seqlength)
  if (is.null(prop)) {
    return(invisible(df))
  } else {
    return(invisible(df[[prop]]))
  }
}