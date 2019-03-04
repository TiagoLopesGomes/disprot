library(protr)
library(readr)
library(dplyr)
library(parallel)

# Notes
# - class: compositional, pseudo-amino acid composition, autocorrelational
# - feature: AAC, CTDC, CTDD, ...
# - descriptor: (e.g., for PAAC) Xc1.A, Ac1.R, Xc1.N, ...

# --- CONSTANTS --- #

# compositional features
featureFuns_comp = c(
  AAC = protr::extractAAC,
  CTDC = protr::extractCTDC,
  CTDD = protr::extractCTDD,
  CTDT = protr::extractCTDT,
  CTriad = protr::extractCTriad,
  DC = protr::extractDC,
  TC = protr::extractTC
)
# pseudo-amino acid composition features
featureFuns_paac = c(
  APAAC = protr::extractAPAAC,
  PAAC = protr::extractPAAC
)
# autocorrelational features
featureFuns_cor = c(
  Geary = protr::extractGeary,
  Moran = protr::extractMoran,
  MoreauBroto = protr::extractMoreauBroto,
  QSO = protr::extractQSO,
  SOCN = protr::extractSOCN
)
featureFuns_all = c(featureFuns_comp, featureFuns_paac, featureFuns_cor)

# --- FUNCTIONS --- #

featurizeSeq <- function(seq, features, lambda, nlag, name_prefix = TRUE) {
  # Featurize a single sequence.
  # 
  # Args
  # - seq: character
  #     sequence to featurize
  # - features: vector
  #     features to include. See names of featureFuns* vectors
  # - lambda: integer
  #     Parameter for pseudo amino acid composition features (PAAC, APAAC)
  # - nlag: integer
  #     Parameter for autocorrelational features (Geary, Moran, MoreauBroto, QSO, SOCN)
  # - name_prefix: character. default = TRUE
  #     Append feature name as prefix to each of its descriptor names to ensure unique
  #     descriptor names.
  # 
  # Returns: vector, character
  #   Named vector
  
  vec = NULL
  for (feature in features) {
    # See https://stackoverflow.com/q/13992367 for conditional inclusion of arguments (lambda, nlag)
    new_vec = do.call(
      featureFuns_all[[feature]],
      list(x = seq, lambda = lambda, nlag = nlag)[c(T, feature %in% names(featureFuns_paac), feature %in% names(featureFuns_cor))])
    if (name_prefix) {
      names(new_vec) = paste(feature, names(new_vec), sep = "_")
    }
    vec = c(vec, new_vec)
  }
  invisible(vec)
}

featurizeSeqs <- function(seqs, features = NULL, lambda = 0, nlag = lambda, format = "matrix", save = NULL, cores = NULL) {
  # Featurize a set of sequences. Exclude sequences with length <= lambda, nlag.
  # 
  # Notes: requires parallel package
  # 
  # Args
  # - seqs: vector, character
  #     named vector of sequences
  # - features: vector, character. default = NULL
  #     Features to use. See the featureFuns map.
  #     If NULL, all features are extracted.
  # - lambda: integer. default = 0
  #     Parameter for pseudo amino acid composition features (PAAC, APAAC)
  #     Must be > 0 for pseudo amino acid composition features to be used.
  # - nlag: integer. default = lambda
  #     Parameter for autocorrelational features (Geary, Moran, MoreauBroto, QSO, SOCN)
  #     Must be > 0 for autocorrelational features to be used.
  # - format: character. default = "matrix"
  #     Class and format of returned featurizations
  #     - "matrix": feature matrix. rows = sequences, columns = features
  #     - "list": each element (feature category) is a named list mapping sequence name 
  #         to a vector of values in that feature category
  # - save: character. default = NULL
  #     Path to save created features database.
  #     If `format` == "matrix": save as TSV file
  #     If `format` == "list": save as RDS file
  # - cores: integer. default = NULL
  #     Number of cores to use. By default, uses parallel::detectCores() - 1
  # 
  # Returns: see `format` argument
  
  # validate arguments
  stopifnot(format %in% c("matrix", "list"))
  stopifnot(is.null(save) || is.character(save))
  stopifnot(is.null(cores) || is.numeric(cores))
  cores = ifelse(is.numeric(cores), cores, parallel::detectCores() - 1)

  # only keep sequences consisting of the 20 canonical amino acids and with length > lambda, nlag
  validSeqs = sapply(seqs, function(x) protr::protcheck(x) & nchar(x) > lambda & nchar(x) > nlag)
  if (sum(validSeqs) < length(seqs)) {
    warning(paste("Excluding", sum(!validSeqs), "out of", length(seqs),
                  "sequences with length <= min(lambda, nlag) or non-canonical AAs from features database.\n"))
  }
  seqs = seqs[validSeqs]
  
  # choose features if no features specified (`features` is NULL)
  if (is.null(features)) {
    features = names(featureFuns_comp)
    if (lambda > 0) {
      features = c(features, names(featureFuns_paac))
    }
    if (nlag > 0) {
      features = c(features, names(featureFuns_cor))
    }
  }
  
  # validate features
  stopifnot(features %in% names(featureFuns_all))
  if (lambda <= 0 && any(features %in% names(featureFuns_paac))) {
    stop("Must supply lambda > 0 to use pseudo amino acid composition descriptors")
  }
  if (nlag <= 0 && any(features %in% names(featureFuns_cor))) {
    stop("Must supply nlag > 0 to use autocorrelational descriptors")
  }

  # compute features
  featuresDB = list()
  if (.Platform$OS.type == "unix") {
    # use multicore-based parallelisation (fork)
    if (format == "list") {
      for (feature in features) {
        featuresDB[[feature]] = parallel::mcmapply(
          FUN = featureFuns_all[[feature]], seqs,
          MoreArgs = list(lambda = lambda, nlag = nlag)[c(feature %in% names(featureFuns_paac), feature %in% names(featureFuns_cor))],
          mc.cores = cores)
      }
    } else {
      featuresDB = t(parallel::mcmapply(FUN = function(seq) featurizeSeq(seq, features, lambda, nlag), seqs, mc.cores = cores))
    }
  } else {
    # use SNOW-based parallelisation (clusters) --> need to export variables to be used by workers
    cl = parallel::makeCluster(cores)
    parallel::clusterExport(cl, varlist = c("featureFuns_all", "featureFuns_paac", "featureFuns_cor"))
    
    if (format == "list") {
      for (feature in features) {
        featuresDB[[feature]] = parallel::clusterMap(
          cl, featureFuns_all[[feature]], seqs,
          MoreArgs = list(lambda = lambda, nlag = nlag)[c(feature %in% names(featureFuns_paac), feature %in% names(featureFuns_cor))])
      }
    } else {
      featuresDB = t(parallel::clusterMap(cl, featurizeSeq, seqs,
                                          MoreArgs = list(features = features, lambda = lambda, nlag = nlag),
                                          SIMPLIFY = TRUE))
    }
    
    parallel::stopCluster(cl)
  }

  if (is.character(save)) {
    if (format == "list") {
      saveRDS(featuresDB, file = save)
    } else {
      readr::write_tsv(x = tibble::as_tibble(featuresDB, rownames = "rowname"), path = save)
    }
  }

  invisible(featuresDB)
}

matchFeatures <- function(dataset, featuresDB,
  dataset_cols = c("p1_uniprotName", "p2_uniprotName"), features = NULL,
  features_path = NULL, labels_path = NULL) {
  # lookup features for entries in dataset from list-format features database (i.e., as created
  # by featurizeSeqs(..., format="list")).
  # 
  # Example: Extract features from featuresDB for genes in data frame disorderSeqs
  #   as.matrix(matchFeatures(dataset = disorderSeqs, featuresDB = featuresDB,
  #                           dataset_cols = "names", features = features))
  # 
  # Args
  # - dataset: data.frame
  #     Dataset of peptides to lookup in `featuresDB`. Must include columns given in `dataset_cols`.
  # - featuresDB: list, list, vector
  #     Database of features of peptides, as returned by createFeaturesDB().
  # - dataset_cols: vector, character. default = c("p1_uniprotName", "p2_uniprotName")
  #     Columns of dataset to use for lookup. If multiple column names are supplied,
  #     feature vectors from each column are simply concatenated.
  # - features: vector, character. default = NULL
  #     Features to use from featuresDB. If NULL, all features are used.
  # - features_path: character. default = NULL
  #     Path to save features file. Compression will be applied automatically based on extension.
  # - labels_path: character. default = NULL
  #     Path to save labels file. Compression will be applied automatically based on extension.
  #     Assumes `dataset` additionally contains "labels" column.
  # 
  # Returns: data.frame
  #   Features data.frame
  
  if (is.null(features)) {
    features = names(featuresDB)
  }
  stopifnot(all(features %in% names(featuresDB)))
  ids = names(featuresDB[[features[1]]])
  stopifnot(all(sapply(featuresDB, function(x) identical(names(x), ids))))

  # lookup features for data from features database
  mat = NULL
  for (feature in features) {
    print(paste0("matching feature ", feature))
    for (col in dataset_cols) {
      # mat = cbind(mat, do.call(rbind, featuresDB[[feature]][match(dataset[[col]], ids)]))
      mat = cbind(mat, do.call(rbind, featuresDB[[feature]][dataset[[col]]]))
    }
    # mat = cbind(mat, do.call(rbind, featuresDB[[feature]][match(dataset[["p2_uniprotName"]], ids)]))
  }
  df = as.data.frame(mat, stringsAsFactors = FALSE)
  # print(dim(df))
  
  # save features and labels files
  if (!is.null(features_path)) {
    readr::write_tsv(df, features_path)
  }
  if (!is.null(labels_path)) {
    readr::write_lines(dataset[["label"]], labels_path)
  }

  invisible(df)
}

main <- function() {
  # setup
  if (.Platform$OS.type == "unix") {
    projectDir = file.path(Sys.getenv("HOME"), "projects", "IDP_PPI_Prediction")
  } else {
    projectDir = "D:/OneDrive/Documents/College (Stanford)/2018-19 Coterm/Q1 Autumn/CS 229/Project/IDP_PPI_Prediction"
  }
  compress_save = TRUE
  useExistingFeaturesDB = TRUE
  createFullFeaturesDB = FALSE

  # read data
  dirDataProc = file.path(projectDir, "data", "proc")
  trainset_files = sort(list.files(dirDataProc, pattern = "trainset\\d+\\.tsv.gz", full.names = FALSE))
  testset_files = sort(list.files(dirDataProc, pattern = "testset\\d+\\.tsv.gz", full.names = FALSE))
  proteome_file = file.path(dirDataProc, "UP000005640_9606.tsv.gz")
  featuresDB_full_file = file.path(dirDataProc, "featuresDB_full.RDS")
  featuresDB_small_file = file.path(dirDataProc, "featuresDB_small.RDS")

  nTrainFiles = length(trainset_files)
  nTestFiles = length(testset_files)

  if (createFullFeaturesDB) {
    featuresDB = createFeaturesDB(proteome_file, save = featuresDB_full_file)
  } else {
    if (useExistingFeaturesDB && file.exists(featuresDB_full_file)) {
      featuresDB = readRDS(featuresDB_full_file)
    } else if (useExistingFeaturesDB && file.exists(featuresDB_small_file)) {
      featuresDB = readRDS(featuresDB_small_file)
    } else {
      # create small features database consisting only of unique identifiers in data
      ids = NULL
      for (i in 1:nTrainFiles) {
        trainset = readr::read_tsv(file.path(dirDataProc, trainset_files[i]))
        ids = c(ids, trainset[["p1_uniprotName"]], trainset[["p2_uniprotName"]])
      }
      for (i in 1:nTestFiles) {
        testset = readr::read_tsv(file.path(dirDataProc, testset_files[i]))
        ids = c(ids, testset[["p1_uniprotName"]], testset[["p2_uniprotName"]])
      }
      ids = sort(unique(ids))
  
      featuresDB = createFeaturesDB(proteome_file, ids = ids, save = featuresDB_small_file)
    }
    ids = names(featuresDB[[1]])
  }

  ext_features = ifelse(compress_save, ".tsv.gz", ".tsv")
  ext_labels = ifelse(compress_save, ".txt.gz", ".txt")

  for (i in 1:nTrainFiles) {
    trainset = readr::read_tsv(file.path(dirDataProc, trainset_files[i]))
    valid = (trainset[["p1_uniprotName"]] %in% ids) & (trainset[["p2_uniprotName"]] %in% ids)
    trainset = trainset[valid, ]
    features_path = file.path(dirDataProc, paste0(sprintf("train_features%02d", i), ext_features))
    labels_path = file.path(dirDataProc, paste0(sprintf("train_labels%02d", i), ext_labels))
    matchFeatures(trainset, featuresDB, features_path = features_path, labels_path = labels_path, compress = compress_save)
  }

  for (i in 1:nTestFiles) {
    testset = readr::read_tsv(file.path(dirDataProc, testset_files[i]))
    valid = (testset[["p1_uniprotName"]] %in% ids) & (testset[["p2_uniprotName"]] %in% ids)
    testset = testset[valid, ]
    features_path = file.path(dirDataProc, paste0(sprintf("test_features%02d", i), ext_features))
    labels_path = file.path(dirDataProc, paste0(sprintf("test_labels%02d", i), ext_labels))
    matchFeatures(testset, featuresDB, features_path = features_path, labels_path = labels_path, compress = compress_save)
  }
}

# main()