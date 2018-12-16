library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(ggplot2)

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
                            proteome_seqCol = "seq", proteome_idCol = "d2p2_id") {
  # Convert data frames into Bioconductor objects
  # - disordered regions -> GRanges
  # - protein sequences -> AAStringSet
  # 
  # Args
  # - disorderPred: data.frame
  # - proteome: data.frame
  # - disorder_idCol, disorder_predCol, proteome_seqCol, proteome_idCol: character
  #     Relevant columns from data frames
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
  
  # Subset AAStringSet by GRanges
  disorderSeqs = proteinDB[disorderRanges]
  
  return(invisible(list(disorderRanges=disorderRanges, disorderSeqs=disorderSeqs)))
}

plotPCA <- function(disorderSeqsList, thresh=0, pcs=c('PC1', 'PC2'), plot=TRUE, verbose=TRUE) {
  if (verbose) {
    print('Number of unique proteins (unfiltered)')
    print(sapply(disorderSeqsList, function(x) length(unique(names(x)))))
    print('Number of disordered regions (unfiltered)')
    print(sapply(disorderSeqsList, length))
  }
  disorderSeqsList = sapply(disorderSeqsList, function(x) x[width(x) > thresh])
  labels = unlist(sapply(names(disorderSeqsList),
                         function(x) rep(x, length(disorderSeqsList[[x]]))))
  if (verbose) {
    print('Number of disordered regions (filtered)')
    print(sapply(disorderSeqsList, length))
  }
  freq_mat = do.call(rbind, sapply(disorderSeqsList, alphabetFrequency))
  p = prcomp(freq_mat)
  data = as.data.frame(p$x)
  data$label = labels
  pca = ggplot(data = data) +
    geom_point(aes_string(x = pcs[1], y = pcs[2], color='label')) + 
    labs(title="PCA of amino acid frequencies",
         subtitle=paste('minimum disordered region length:', thresh))
  if (plot) {
    print(pca)
  }
  invisible(pca)
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