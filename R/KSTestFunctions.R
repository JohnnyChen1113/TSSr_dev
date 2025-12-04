#' KS Test Helper Functions
#'
#' @description Internal helper functions for Kolmogorov-Smirnov test analysis
#' of TSS distributions between samples.
#'
#' @version 2.0.0
#' @date 2025-12-03


#' Get TSS distribution within a cluster region
#'
#' @param tss_matrix TSSrawMatrix data.table
#' @param sample_name Sample name (column name in the matrix)
#' @param target_chr Chromosome
#' @param target_strand Strand
#' @param start Start position
#' @param end End position
#' @return A vector of TSS positions (repeated by raw read counts)
#' @keywords internal
.get_tss_distribution <- function(tss_matrix, sample_name, target_chr, target_strand, start, end) {

  # Check if sample_name exists as a column
  if (!(sample_name %in% names(tss_matrix))) {
    return(NULL)
  }

  # Filter TSS positions within the cluster region
  tss_in_region <- tss_matrix[chr == target_chr & strand == target_strand &
                               pos >= start & pos <= end]

  if (nrow(tss_in_region) == 0) {
    return(NULL)
  }

  # Get positions and their raw counts
  positions <- tss_in_region$pos
  counts <- tss_in_region[[sample_name]]

  # Only keep positions with counts > 0
  valid_idx <- counts > 0
  if (sum(valid_idx) == 0) {
    return(NULL)
  }

  positions <- positions[valid_idx]
  counts <- counts[valid_idx]

  # Repeat each position by its raw count (no weighting)
  tss_positions <- rep(positions, times = counts)

  return(tss_positions)
}


#' Compute global unified boundaries across all samples
#'
#' @param consensus_clusters List of consensus clusters from TSSr object
#' @return A data.table with unified boundaries for each cluster
#' @keywords internal
.compute_unified_boundaries <- function(consensus_clusters) {

  sample_names <- names(consensus_clusters)

  # Collect all cluster data from all samples
  all_clusters_list <- lapply(sample_names, function(sample) {
    dt <- consensus_clusters[[sample]]
    if (is.null(dt) || nrow(dt) == 0) return(NULL)

    if (!inherits(dt, "data.table")) dt <- as.data.table(dt)

    data.table(
      cluster = as.character(dt$cluster),
      chr = as.character(dt$chr),
      strand = as.character(dt$strand),
      start = dt$start,
      end = dt$end
    )
  })

  all_clusters <- rbindlist(all_clusters_list, use.names = TRUE)

  # Compute unified boundaries: min(start) and max(end) for each cluster
  unified_boundaries <- all_clusters[, .(
    unified_start = min(start),
    unified_end = max(end)
  ), by = .(cluster, chr, strand)]

  unified_boundaries[, cluster_key := paste(cluster, chr, strand, sep = "_")]

  return(unified_boundaries)
}


#' Process a single cluster for KS test
#'
#' @param key Cluster key
#' @param dt1 Data table for sample 1
#' @param dt2 Data table for sample 2
#' @param unified_boundaries Unified boundaries data table
#' @param tss_matrix TSSrawMatrix
#' @param sample1 Sample 1 name
#' @param sample2 Sample 2 name
#' @param min_tags Minimum tags required
#' @return A data.table row with KS test results
#' @keywords internal
.process_ks_cluster <- function(key, dt1, dt2, unified_boundaries, tss_matrix,
                                 sample1, sample2, min_tags) {
  row1 <- dt1[cluster_key == key][1]
  row2 <- dt2[cluster_key == key][1]

  # Get the GLOBAL unified boundary (same for all samples)
  boundary <- unified_boundaries[key]
  unified_start <- boundary$unified_start
  unified_end <- boundary$unified_end

  # Get TSS distributions using global unified boundary
  dist1 <- .get_tss_distribution(
    tss_matrix, sample1,
    target_chr = row1$chr,
    target_strand = row1$strand,
    start = unified_start,
    end = unified_end
  )

  dist2 <- .get_tss_distribution(
    tss_matrix, sample2,
    target_chr = row2$chr,
    target_strand = row2$strand,
    start = unified_start,
    end = unified_end
  )

  # Initialize result
  result <- data.table(
    cluster = row1$cluster,
    chr = row1$chr,
    start = unified_start,
    end = unified_end,
    strand = row1$strand,
    dominant_tss_1 = row1$dominant_tss,
    dominant_tss_2 = row2$dominant_tss,
    tags_1 = row1$tags,
    tags_2 = row2$tags,
    ks_D = NA_real_,
    ks_pvalue = NA_real_,
    ks_test_status = "not_tested"
  )

  # Perform KS test if both distributions have enough data
  if (!is.null(dist1) && !is.null(dist2) &&
      length(dist1) >= min_tags && length(dist2) >= min_tags) {

    tryCatch({
      # Suppress "ties" warning - this is expected for discrete position data
      ks_result <- suppressWarnings(ks.test(dist1, dist2))
      result$ks_D <- ks_result$statistic
      result$ks_pvalue <- ks_result$p.value
      result$ks_test_status <- "success"
    }, error = function(e) {
      result$ks_test_status <- paste("error:", e$message)
    })

  } else if (is.null(dist1) || is.null(dist2)) {
    result$ks_test_status <- "no_tss_data"
  } else {
    result$ks_test_status <- "insufficient_tags"
  }

  return(result)
}


#' Perform KS test between two samples
#'
#' @param consensus_clusters List of consensus clusters
#' @param tss_matrix TSSrawMatrix
#' @param sample1 Sample 1 name
#' @param sample2 Sample 2 name
#' @param unified_boundaries Pre-computed unified boundaries
#' @param min_tags Minimum tags required
#' @param useMultiCore Use parallel processing
#' @param numCores Number of cores
#' @return A data.table with KS test results
#' @keywords internal
.perform_ks_test <- function(consensus_clusters, tss_matrix, sample1, sample2,
                              unified_boundaries, min_tags, useMultiCore, numCores) {

  dt1 <- consensus_clusters[[sample1]]
  dt2 <- consensus_clusters[[sample2]]

  # Convert to data.table
  if (!inherits(dt1, "data.table")) dt1 <- as.data.table(dt1)
  if (!inherits(dt2, "data.table")) dt2 <- as.data.table(dt2)

  # Create cluster keys
  dt1[, cluster_key := paste(cluster, chr, strand, sep = "_")]
  dt2[, cluster_key := paste(cluster, chr, strand, sep = "_")]

  # Find common clusters between the two samples AND in unified_boundaries
  common_keys <- intersect(intersect(dt1$cluster_key, dt2$cluster_key),
                           unified_boundaries$cluster_key)

  message(paste("Found", length(common_keys), "common clusters between", sample1, "and", sample2))

  # Set key for fast lookup
  setkey(unified_boundaries, cluster_key)

  # Process clusters: parallel or sequential
  if (useMultiCore) {
    if (is.null(numCores)) {
      numCores <- detectCores()
    }
    print(paste("process is running on", numCores, "cores..."))

    results_list <- mclapply(common_keys, function(key) {
      .process_ks_cluster(key, dt1, dt2, unified_boundaries, tss_matrix,
                          sample1, sample2, min_tags)
    }, mc.cores = numCores)

  } else {
    results_list <- lapply(common_keys, function(key) {
      .process_ks_cluster(key, dt1, dt2, unified_boundaries, tss_matrix,
                          sample1, sample2, min_tags)
    })
  }

  results <- rbindlist(results_list)

  # Calculate adjusted p-values (FDR)
  valid_pvalues <- !is.na(results$ks_pvalue)
  if (sum(valid_pvalues) > 0) {
    results[valid_pvalues, ks_padj := p.adjust(ks_pvalue, method = "BH")]
  } else {
    results[, ks_padj := NA_real_]
  }

  # Sort by cluster
  setorder(results, chr, start)

  return(results)
}
