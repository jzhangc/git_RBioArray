#' @title rbio_randomforest_fs
#'
#' @description Supervised feature selection with recursive random forest (RFFS)
#' @param object \code{sig} object. Input \code{sig} class object.
#' @param sample_id.var String. Sample id variable name, present in \code{object$input_data$targets}.
#' @param sample_group.car String. Sample group variable name, present in \code{object$input_data$targets}.
#' @param export.name String. Prefix for export files and objects. When \code{NULL}, the function uses the input object name.
#' @param gene_symbol.only Boolean. Whether or not to remove probes without gene symbol. Default is \code{FALSE}.
#' @param fs_on_de Boolean. Whether to conduct RFFS on only the DE targets (features), according to the iput object.
#' @param rffs.center_scale Boolean. Whether to center.scale the data for RFFS. Default is \code{FALSE}.
#' @param rffs.quantile Boolean. Whether to quantile normalize the data for RFFS. Default is \code{FALSE}.
#' @param rffs.ntimes Integer. Number of random forest recursion. Default is \code{50}.
#' @param rffs.inital_fs.ntree Integer. ntree setting for the RFFS initial selection. Default is \code{501}.
#' @param rffs.inital_fs.mtry Integer. mtry setting for the RFFS initial selection. Default is \code{NULL}.
#' @param rffs.sfs Boolean. Whether to run the secondary sequential forward selection (SFS). Default is \code{FALSE}.
#' @param rffs.sfs_ntree Integer. Set when \code{rffs.sfs = TRUE}, ntree setting for the RFFS SFS step. Default is \code{501}.
#' @param rffs.randomstate Integer. Set random state for RFFS. If \code{0}, the random state is off. Default is \code{0}.
#' @param parallelComputing Whether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Set when \code{parallelComputing = TRUE}, number of CPU threads/cores for parallel computing.
#' @param cluterType clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @return Export an \code{rbio_rffs} object to the global environment.
#' @details The resulting \code{rbio_rffs} object contains the following items:
#'
#'          \code{ml_raw_dataframe} Raw data dataframe for RFFS. Data not normalized or center.scaled by the function.
#'                                    The data is usually previously normalized or processed.
#'                                    The E data is from \code{object$input_data$norm_E} from the input \code{sig} object.
#'                                    It is noted the data may be subset with only the DE features according to the \code{fs_on_de} setting.
#'                                    This data frame can be used as input data for other ML purpouses featured in the \code{RBioFS} packages.
#'
#'          \code{rffs_working_data} It is list containing the follow items:
#'              \code{rffs_working_E} Working datafarme as the \code{x} input for RFFS. Columns: features; Rows: samples. Normalized/center.scaled when set.
#'              \code{rffs_working_y}  Working vector (\code{factor}) as the \code{y} input for RFFS.#'
#'              \code{center_scale} Center.scale setting.
#'              \code{quantile} Quantile setting.
#'
#'          \code{rffs_initial_results}
#'
#'          \code{rffs_sfs_results}
#'
#'          \code{rffs_subset_dataframe}: Input data subset with only the RFFS selected features. No data or fas
#'
#'          \code{rffs_ntimes}
#'
#'          \code{rffs_inital_ntree}
#'
#'          \code{rffs_inital_mtry}
#'
#'          \code{rffs_sfs}
#'
#'          \code{rffs_sfs_ntree}
#'
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom RBioFS rbioFS_rf_initialFS rbioFS_rf_SFS center_scale
#' @export
rbio_randomforest_fs <- function(object, sample_id.var = NULL, sample_group.var = NULL,
                                 export.name = NULL,
                                 gene_symbol.only = FALSE,
                                 fs_on_de = FALSE,
                                 rffs.center_scale = FALSE, rffs.quantile = FALSE,
                                 rffs.ntimes = 50,
                                 rffs.inital_fs.ntree = 501, rffs.inital_fs.mtry = NULL,
                                 rffs.sfs =  FALSE, rffs.sfs_ntree = 501,
                                 rffs.randomstate = 0,
                                 parallelComputing = TRUE, n_cores = parallel::detectCores() - 1, clusterType = "FORK") {
  # - check arguments -
  if (!any(class(object) %in% c("sig"))) stop("The input object needs to be an \"sig\" class.")
  if (is.null(sample_id.var)) {
    stop("Please provide sample_id.var.")
  } else if (!sample_id.var %in% names(object$input_data$targets)) {
    stop("Provided sample_id.var not found in the input object (check: names(object$input_data$targets)).")
  }
  if (is.null(sample_group.var)) {
    stop("Please provide sample_group.var.")
  } else if (!sample_group.var %in% names(object$input_data$targets)) {
    stop("Provided sample_group.var not found in the input object (check: names(object$input_data$targets)).")
  }

  # - random state -
  if (rffs.randomstate) set.seed(rffs.randomstate)

  # - set up data vars -
  E <- object$input_data$norm_E
  genes <- object$input_data$genes
  sample_id <- object$input_data$targets[, sample_id.var]
  sample_group <- object$input_data$targets[, sample_group.var]
  sample_group <- factor(sample_group, levels = unique(sample_group))

  input.genes_annotation.gene_symbol.var_name <- object$input_data$input.genes_annotation.gene_symbol.var_name
  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(object$input_data$genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  } else if (is.null(input.genes_annotation.gene_symbol.var_name)) {
    cat("input.genes_annotation.gene_symbol.var_name is NULL, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  input.genes_annotation.gene_id.var_name = object$input_data$input.genes_annotation.gene_id.var_name

  if (is.null(export.name)) {
    cat("Object/file export name prefix automatically set to the input object name. \n")
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }
  input.sample_groups <- object$input_data$sample_groups
  input.genes_annotation.control_type <- object$input_data$input.genes_annotation.control_type
  thresholding_summary <- object$thresholding_summary
  row.lab.var_name <- input.genes_annotation.gene_id.var_name

  # - finalize settings -
  n_de <- as.numeric(object$significant_change_summary[object$significant_change_summary[, "comparisons"] == "F stats", "True"])

  if (fs_on_de & n_de < 2) {
    warning("fs_on_de = TRUE, however, the number of DE targts < 2, automatically reset fs_on_de = FALSE\n")
    fs_on_de <- FALSE
  }

  if (is.null(input.genes_annotation.control_type)) {
    cat("Argument input.genes_annotation.control_type is NULL, no control probes are removed.\n\n")
    rm.control <- FALSE
  } else {
    rm.control <- TRUE
  }

  # - set up working data -
  dfm <-  data.frame(genes, E, check.names = FALSE)

  # if set, set up control filter
  if (rm.control) {
    is.control <- dfm[, input.genes_annotation.control_type$control_type.var_name] == input.genes_annotation.control_type$exp_type.value
    dfm <- dfm[is.control, ]
    cat(paste0("Number of control targets removed: ", length(which(is.control))))
  } else {
    is.control <- rep(FALSE, times = nrow(dfm))  # if not set, the function assumes none of the targets are controls
  }

  # if set, set up gene symbol filter
  if (gene_symbol.only) {
    cat("gene_symbol.only = TRUE, any genes without a gene symbol are removed from clustering.\n\n")
    gene_symbol <- dfm[, input.genes_annotation.gene_symbol.var_name]
    has.gene_symbol <- complete.cases(gene_symbol) & gene_symbol != ""
  } else {
    has.gene_symbol <- rep(TRUE, times = nrow(dfm))   # if not set, the function assumes all targets have symbols
  }

  # set up the initial working data
  cat(paste0("Only running on the DE targets (based on F stats): ", fs_on_de, "\n"))
  if (fs_on_de) {
    is_test <- object$thresholding_summary$F_stats  & !is.control & has.gene_symbol
    n_target_working <- length(which(is_test))
    if (n_target_working < 2) {
      warning("Number of filtered final targetes < 2, revert to using all targets.\n")
      is_test <- rep(TRUE, times = nrow(dfm))
    }
  } else {
    is_test <- rep(TRUE, times = nrow(dfm))  # if to use all targets
  }

  dfm_working <- dfm[is_test, ]
  E_working <- E[is_test, ]

  if (gene_symbol.only) {
    rownames(E_working) <- dfm_working[, input.genes_annotation.gene_symbol.var_name]
  } else {
    rownames(E_working) <- dfm_working[, input.genes_annotation.gene_id.var_name]
  }

  if (rffs.center_scale){
    cat(paste0("Data centered with scaling prior to RF-FS.\n"))
    centered_X <- center_scale(E_working, scale = TRUE)  # center data with the option of scaling
    E_working <- centered_X$centerX
  }

  if (rffs.quantile){
    cat(paste("Quantile normalization...", sep = ""))  # initial message
    E_working <- t(fs_data)
    E_working <- RBioFS::rbioNorm(E_working, correctBG = FALSE)
    E_working <- t(fs_data)
    E_working <- data.frame(E_working, check.names = FALSE)
    cat(paste("Done!\n", sep = ""))  # final message
  }

  E_working <- t(E_working)
  # dim(E_working)
  if (any(colSums(is.na(E_working)) != 0)) {
    warning("Predictors with NA values are removed (hint: check center.scale settings to avoid this).")
    E_working <- E_working[ , colSums(is.na(E_working)) == 0]
  }

  # dim(E_working)
  E_rffs_raw <- E_working  # for output
  rffs_raw_dfm <- data.frame(sample_id = sample_id,
                             sample_group = sample_group,
                             E_rffs_raw, check.names = FALSE)  # for output

  # write.csv(E_working, file = "test.csv", row.names = FALSE)


  # - rRF-FS -
  # initial selection
  cat(paste0("RFFS initial selection..." ))  # initial message
  if (is.null(rffs.inital_fs.mtry)) {
    initial_fs_mtry = if (!is.factor(sample_group)) max(floor(ncol(E_working)/3), 1) else floor(sqrt(ncol(E_working)))
  } else {
    initial_fs_mtry = rffs.inital_fs.mtry
  }
  rbioFS_rf_initialFS(objTitle = export.name, x = E_working, y = sample_group,
                      nTimes = rffs.ntimes, nTree = rffs.inital_fs.ntree,
                      mTry = initial_fs_mtry,
                      parallelComputing = parallelComputing, n_cores = n_cores, clusterType = clusterType,
                      plot = FALSE, verbose = FALSE) # initial FS
  rffs.initial_res <- get(paste0(export.name, "_initial_FS"))
  fs <- rffs.initial_res$feature_initial_FS
  cat(paste("Done!\n", sep = ""))  # final message


  # if set, sfs
  if (rffs.sfs) {
    cat(paste0("RFFS SFS..." ))  # initial message
    rbioFS_rf_SFS(objTitle = export.name,
                  x = rffs.initial_res$training_initial_FS,
                  y = sample_group, nTimes = rffs.ntimes, nTree = rffs.sfs_ntree,
                  parallelComputing = parallelComputing, n_cores = n_cores, clusterType = clusterType,
                  plot = FALSE, verbose = FALSE) # SFS
    rffs.sfs_res <- get(paste0(export.name, "_SFS"))
    fs <- rffs.sfs_res$selected_features
    cat(paste("Done!\n", sep = ""))  # final message
  } else {
    rffs.sfs_res <- NULL
  }

  # - output -
  rffs_subset_dfm <- data.frame(sample_id = sample_id,
                                sample_group = sample_group,
                                E_rffs_raw[, fs, drop = FALSE])  # for output

  out <- list(
    ml_raw_dataframe = rffs_raw_dfm,
    rffs_working_data = list(
      rffs_working_E = E_working,
      rffs_working_y = sample_group,
      center_scale = rffs.center_scale,
      quantile = rffs.quantile
    ),
    rffs_initial_results = rffs.initial_res,
    rffs_sfs_results = rffs.sfs_res,
    rffs_subset_dataframe = rffs_subset_dfm,
    rffs_ntimes = rffs.ntimes,
    rffs_inital_ntree = rffs.inital_fs.ntree,
    rffs_inital_mtry = initial_fs_mtry,
    rffs_sfs = rffs.sfs, rffs_sfs_ntree = rffs.sfs_ntree,
    rffs_randomstate = ifelse(rffs.randomstate, rffs.randomstate, "not set")
  )

  class(out) <- "rbio_rffs"
  assign(paste0(export.name, "_rffs"), out, envir = .GlobalEnv)
}
