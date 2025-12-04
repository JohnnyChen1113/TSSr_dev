#' Kolmogorov-Smirnov Test for TSS Distribution Comparison
#'
#' @description Performs Kolmogorov-Smirnov test to compare TSS distributions
#' between samples within consensus clusters.
#'
#' @version 2.0.0
#' @date 2025-12-03


#' Perform KS test on consensus clusters between sample pairs
#'
#' @description Compares TSS distributions between two samples using the
#' Kolmogorov-Smirnov test. For each consensus cluster, the function extracts
#' TSS positions from TSSrawMatrix and performs a two-sample KS test.
#'
#' @usage ksTestCluster(object, comparePairs, min_tags = 3,
#'                      useMultiCore = FALSE, numCores = NULL)
#'
#' @param object A TSSr object with consensusClusters and TSSrawMatrix.
#' @param comparePairs A list of sample pairs for comparison.
#'        Example: list(c("treatment", "control"), c("sample1", "sample2"))
#' @param min_tags Minimum total tags required in both samples for a cluster
#'        to be tested. Default = 3.
#' @param useMultiCore Logical, whether to use parallel processing. Default = FALSE.
#' @param numCores Number of cores for parallel processing. If NULL and
#'        useMultiCore = TRUE, will use detectCores().
#' @return The TSSr object is updated with KSTestResults slot containing
#'         a list of data.tables with KS test results for each comparison.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes global unified boundaries for each cluster across ALL samples
#'   \item For each cluster, extracts TSS positions using raw counts from TSSrawMatrix
#'   \item Performs two-sample KS test comparing the TSS distributions
#'   \item Calculates FDR-adjusted p-values using Benjamini-Hochberg method
#' }
#'
#' Results include:
#' \itemize{
#'   \item cluster: Cluster identifier
#'   \item chr, start, end, strand: Genomic coordinates (unified boundaries)
#'   \item dominant_tss_1, dominant_tss_2: Dominant TSS positions in each sample
#'   \item tags_1, tags_2: Tag counts in each sample
#'   \item ks_D: KS test D statistic
#'   \item ks_pvalue: KS test p-value
#'   \item ks_padj: FDR-adjusted p-value
#'   \item ks_test_status: Test status (success, no_tss_data, insufficient_tags, error)
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(exampleTSSr)
#' # ksTestCluster(exampleTSSr, comparePairs = list(c("treat", "control")))
#' }
setGeneric("ksTestCluster", function(object, comparePairs,
                                      min_tags = 3,
                                      useMultiCore = FALSE,
                                      numCores = NULL)
  standardGeneric("ksTestCluster"))

#' @rdname ksTestCluster
#' @export
setMethod("ksTestCluster", signature(object = "TSSr"), function(object,
                                                                  comparePairs,
                                                                  min_tags,
                                                                  useMultiCore,
                                                                  numCores) {

  message("\n=== Cluster KS Test Analysis ===\n")

  objName <- deparse(substitute(object))

  # Validate input
  if (length(object@consensusClusters) == 0) {
    stop("Error: consensusClusters is empty. Please run consensusCluster() first.")
  }

  if (is.null(object@TSSrawMatrix) || nrow(object@TSSrawMatrix) == 0) {
    stop("Error: TSSrawMatrix is empty.")
  }

  consensus_clusters <- object@consensusClusters
  tss_matrix <- object@TSSrawMatrix

  # Step 1: Compute global unified boundaries across ALL samples
  message("Step 1: Computing global unified boundaries across all samples...")
  unified_boundaries <- .compute_unified_boundaries(consensus_clusters)
  message(paste("Computed unified boundaries for", nrow(unified_boundaries), "clusters\n"))

  # Step 2: Perform KS test for each comparison pair
  ks_results <- lapply(as.list(seq(comparePairs)), function(i) {

    sample1 <- comparePairs[[i]][1]
    sample2 <- comparePairs[[i]][2]

    # Validate samples exist
    if (!(sample1 %in% names(consensus_clusters))) {
      stop(paste("Error: Sample", sample1, "not found in consensusClusters"))
    }
    if (!(sample2 %in% names(consensus_clusters))) {
      stop(paste("Error: Sample", sample2, "not found in consensusClusters"))
    }

    message(paste("======", sample1, "vs", sample2, "======"))

    results <- .perform_ks_test(
      consensus_clusters = consensus_clusters,
      tss_matrix = tss_matrix,
      sample1 = sample1,
      sample2 = sample2,
      unified_boundaries = unified_boundaries,
      min_tags = min_tags,
      useMultiCore = useMultiCore,
      numCores = numCores
    )

    n_tested <- sum(results$ks_test_status == "success")
    n_sig <- sum(results$ks_padj < 0.05, na.rm = TRUE)
    message(paste("KS test completed:", n_tested, "clusters tested,",
                  n_sig, "significant (padj < 0.05)\n"))

    return(results)
  })

  # Name results by comparison pairs
  result_names <- sapply(as.list(seq(comparePairs)), function(i) {
    paste0(comparePairs[[i]][1], "_VS_", comparePairs[[i]][2])
  })
  names(ks_results) <- result_names

  # Store results in object
  object@KSTestResults <- ks_results

  message("=== All KS tests completed ===")

  assign(objName, object, envir = parent.frame())
})


#' Export KS test results to file
#'
#' @description Exports KS test results to tab-separated text files.
#'
#' @usage exportKSTestTable(object, outputDir = "./")
#'
#' @param object A TSSr object with KSTestResults.
#' @param outputDir Directory to save output files. Default = "./"
#' @return Exports tab-separated files for each comparison.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(exampleTSSr)
#' # exportKSTestTable(exampleTSSr, outputDir = "./ks_results/")
#' }
setGeneric("exportKSTestTable", function(object, outputDir = "./")
  standardGeneric("exportKSTestTable"))

#' @rdname exportKSTestTable
#' @export
setMethod("exportKSTestTable", signature(object = "TSSr"), function(object, outputDir) {

  if (length(object@KSTestResults) == 0) {
    stop("Error: No KS test results found. Please run ksTestCluster() first.")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  ks_results <- object@KSTestResults

  for (name in names(ks_results)) {
    output_file <- file.path(outputDir, paste0(name, ".ks_test.txt"))
    fwrite(ks_results[[name]], output_file, sep = "\t", quote = FALSE)
    message(paste("Exported:", output_file))
  }

  message(paste("\nAll KS test results exported to:", outputDir))
})
