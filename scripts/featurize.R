library(protr)
library(readr)
library(dplyr)
library(parallel)

createFeaturesDB <- function(proteome, seq_col = "seq", name_col = "uniprotName", ids = NULL,
                             features = NULL, lambda = 0, nlag = lambda, save = NULL, cores = NULL) {
  # Create features database. Feature extraction on each unique sequence. Only keeps sequences
  # consisting of the 20 canonical amino acids and with length > lambda
  # 
  # Notes: requires parallel package
  # 
  # Args
  # - proteome: data.frame or character
  #     Must contain `seq_col` and `name_col` columns. The `name_col` column is assumed to
  #     consist of unique values. If a character, assumed to be the path to a tab-separated
  #     proteome data frame file.
  # - seq_col: character. default = "seq"
  #     Column in `proteome` containing sequences
  # - name_col: character. default = "uniprotName"
  #     Column in `proteome` containing sequence identifiers
  # - ids: vector, character. default = NULL
  #     Sequence identifiers to include. If NULL, include all peptides in `proteome`.
  # - features: vector, character. default = NULL
  #     Features to use. See the featureFuns map in the code below.
  #     If NULL, all features are extracted.
  # - lambda: integer. default = 0
  #     Parameter for pseudo amino acid composition descriptors (PAAC, APAAC)
  #     Must be > 0 for pseudo amino acid composition descriptors to be used.
  # - nlag: integer. default = lambda
  #     Parameter for autocorrelational descriptors (Geary, Moran, MoreauBroto, QSO, SOCN)
  #     Must be > 0 for autocorrelational descriptors to be used.
  # - save: character. default = NULL
  #     Path to save created features database.
  # - cores: integer. default = NULL
  #     Number of cores to use. By default, uses parallel::detectCores() - 1
  # 
  # Returns: list, list, vector
  #   Features database. Each element (feature category) is a named list mapping peptide name
  #   to a vector of values in that feature category.

  if (is.character(proteome)) {
    proteome = readr::read_tsv(proteome)
  }
  
  # only keep sequences consisting of the 20 canonical amino acids and with length > lambda
  validSeqs = sapply(proteome$seq, function(x) protr::protcheck(x) & nchar(x) > lambda)
  if (sum(validSeqs) != nrow(proteome)) {
    warning(paste0("Excluding ", sum(!validSeqs), " sequences with length <= lambda from features database."))
  }
  proteome = proteome[validSeqs,]
  
  featureFuns_0 = c(AAC = protr::extractAAC,
                    CTDC = protr::extractCTDC,
                    CTDD = protr::extractCTDD,
                    CTDT = protr::extractCTDT,
                    CTriad = protr::extractCTriad,
                    DC = protr::extractDC,
                    TC = protr::extractTC)
  
  featureFuns_lambda = c(APAAC = protr::extractAPAAC,
                         PAAC = protr::extractPAAC)
  
  featureFuns_nlag = c(Geary = protr::extractGeary,
                       Moran = protr::extractMoran,
                       MoreauBroto = protr::extractMoreauBroto,
                       QSO = protr::extractQSO,
                       SOCN = protr::extractSOCN)
  
  # choose features if no features specified (`features` is NULL)
  if (is.null(features)) {
    features = names(featureFuns_0)
    if (lambda > 0) {
      features = c(features, names(featureFuns_lambda))
    }
    if (nlag > 0) {
      features = c(features, names(featureFuns_nlag))
    }
  }
  
  # validate features
  stopifnot(features %in% c(names(featureFuns_0), names(featureFuns_lambda), names(featureFuns_nlag)))
  if (lambda <= 0 && any(features %in% names(featureFuns_lambda))) {
    stop("Must supply lambda > 0 to use pseudo amino acid composition descriptors")
  }
  if (nlag <= 0 && any(features %in% names(featureFuns_nlag))) {
    stop("Must supply nlag > 0 to use autocorrelational descriptors")
  }
  
  # choose sequences to featurize by id
  if (is.null(ids)) {
    seqs = proteome[[seq_col]]
    names(seqs) = proteome[[name_col]]
  } else {
    idx = match(ids, proteome[[name_col]])
    valid = !is.na(idx)
    ids = ids[valid]
    seqs = proteome[["seq"]][idx[valid]]
    names(seqs) = ids
  }

  featuresDB = list()
  if (.Platform$OS.type == "unix") {
    if (!is.numeric(cores)) {
      cores = parallel::detectCores() - 1
    }
    for (feature in intersect(names(featureFuns_0), features)) {
      featuresDB[[feature]] = parallel::mclapply(X = seqs, FUN = featureFuns_0[[feature]], mc.cores = cores)
    }
    for (feature in intersect(names(featureFuns_lambda), features)) {
      featuresDB[[feature]] = parallel::mclapply(X = seqs, FUN = function(x) featureFuns_lambda[[feature]](x, lambda = lambda), mc.cores = cores)
    }
    for (feature in intersect(names(featureFuns_nlag), features)) {
      featuresDB[[feature]] = parallel::mclapply(X = seqs, FUN = function(x) featureFuns_nlag[[feature]](x, nlag = nlag), mc.cores = cores)
    }
    # featuresDB[["aac"]] = parallel::mclapply(X = seqs, FUN = protr::extractAAC, mc.cores = cores)
    # featuresDB[["dc"]] = parallel::mclapply(X = seqs, FUN = protr::extractDC, mc.cores = cores)
    # featuresDB[["tc"]] = parallel::mclapply(X = seqs, FUN = protr::extractTC, mc.cores = cores)
    # featuresDB[["ctdc"]] = parallel::mclapply(X = seqs, FUN = protr::extractCTDC, mc.cores = cores)
    # featuresDB[["ctdt"]] = parallel::mclapply(X = seqs, FUN = protr::extractCTDT, mc.cores = cores)
    # featuresDB[["ctdd"]] = parallel::mclapply(X = seqs, FUN = protr::extractCTDD, mc.cores = cores)
    # featuresDB[["ctriad"]] = parallel::mclapply(X = seqs, FUN = protr::extractCTriad, mc.cores = cores)
    # featuresDB[["paac"]] = parallel::mclapply(X = seqs, FUN = function(x) protr::extractPAAC(x, lambda = lambda), mc.cores = cores)
  } else {
    if (is.numeric(cores)) {
      cl = parallel::makeCluster(cores)
    } else {
      cl = parallel::makeCluster(parallel::detectCores() - 1)
    }
    
    for (feature in intersect(names(featureFuns_0), features)) {
      featuresDB[[feature]] = parallel::clusterApply(cl, seqs, featureFuns_0[[feature]])
      names(featuresDB[[feature]]) = names(seqs)
    }
    for (feature in intersect(names(featureFuns_lambda), features)) {
      featuresDB[[feature]] = parallel::clusterApply(cl, seqs, function(x) featureFuns_lambda[[feature]](x, lambda = lambda))
      names(featuresDB[[feature]]) = names(seqs)
    }
    for (feature in intersect(names(featureFuns_nlag), features)) {
      featuresDB[[feature]] = parallel::clusterApply(cl, seqs, function(x) featureFuns_nlag[[feature]](x, nlag = nlag))
      names(featuresDB[[feature]]) = names(seqs)
    }
    # featuresDB[["aac"]] = parallel::clusterApply(cl, seqs, protr::extractAAC)
    # names(featuresDB[["aac"]]) = names(seqs)
    # featuresDB[["dc"]] = parallel::clusterApply(cl, seqs, protr::extractDC)
    # names(featuresDB[["dc"]]) = names(seqs)
    # featuresDB[["tc"]] = parallel::clusterApply(cl, seqs, protr::extractTC)
    # names(featuresDB[["tc"]]) = names(seqs)
    # featuresDB[["ctdc"]] = parallel::clusterApply(cl, seqs, protr::extractCTDC)
    # names(featuresDB[["ctdc"]]) = names(seqs)
    # featuresDB[["ctdt"]] = parallel::clusterApply(cl, seqs, protr::extractCTDT)
    # names(featuresDB[["ctdt"]]) = names(seqs)
    # featuresDB[["ctdd"]] = parallel::clusterApply(cl, seqs, protr::extractCTDD)
    # names(featuresDB[["ctdd"]]) = names(seqs)
    # featuresDB[["ctriad"]] = parallel::clusterApply(cl, seqs, protr::extractCTriad)
    # names(featuresDB[["ctriad"]]) = names(seqs)
    # featuresDB[["paac"]] = parallel::clusterApply(cl, seqs, function(x) protr::extractPAAC(x, lambda = lambda))
    # names(featuresDB[["paac"]]) = names(seqs)
  }

  if (is.character(save)) {
    saveRDS(featuresDB, file = save)
  }

  invisible(featuresDB)
}

matchFeatures <- function(dataset, featuresDB,
  dataset_cols = c("p1_uniprotName", "p2_uniprotName"), features = NULL,
  features_path = NULL, labels_path = NULL) {
  # lookup features for entries in dataset from features database.
  # 
  # Args
  # - dataset: data.frame
  #     Dataset of peptides to lookup in `featuresDB`. Must include columns given in `dataset_cols`.
  # - featuresDB: list, list, vector
  #     Database of features of peptides, as returned by createFeaturesDB().
  # - dataset_cols: vector, character. default = c("p1_uniprotName", "p2_uniprotName")
  #     Columns of dataset to use for lookup.
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

featurize <- function(disorderSeqs, features, ...) {
  # Create feature matrix for a set of sequences
  # 
  # Args
  # - disorderSeqs: AAStringSet or data.frame
  #     Set of sequences. Assumed to already be length/width-filtered.
  #     If data.frame, must have columns "seq" and "names", and "names" should consist of
  #     unique values (i.e., no duplicates).
  # - features: vector, character
  #     Features to include. See createFeaturesDB()
  # - ...
  #     Arguments to pass onto createFeaturesDB(), especially "lambda" and "nlag"
  # 
  # Returns: matrix
  #   Rows represent sequences. Columns are different features.
  
  if (inherits(disorderSeqs, 'AAStringSet')) {
    disorderSeqs = as.data.frame(disorderSeqs) %>%
      dplyr::rename(seq = x) %>%
      dplyr::mutate(names = as.character(1:length(disorderSeqs))) %>%
      tibble::as.tibble()
  }
  featuresDB = createFeaturesDB(proteome = disorderSeqs, name_col = "names", features = features, ...)
  featuresMat = as.matrix(matchFeatures(dataset = disorderSeqs, featuresDB = featuresDB,
                                        dataset_cols = "names", features = features))
  invisible(featuresMat)
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