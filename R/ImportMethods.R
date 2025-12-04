###############################################################################
#' Precisely identify TSSs from bam files, paired end bam files, bed files,
#' BigWig files, tss files, or tss tables.
#'
#' @description getTSS function is used to precisely identify TSSs from multiple
#' input file formats. The files include users' home-made alignment files (bam format)
#' or downloaded files from public databases. See inputFilesType for details on
#' the supported input file formats.
#'
#' @usage getTSS(object, sequencingQualityThreshold = 10,
#' mappingQualityThreshold = 20, softclippingAllowed = FALSE,
#' useMultiCore = FALSE, numCores = NULL, readsPerPiece = 1000000)
#'
#' @param object A TSSr object.
#' @param sequencingQualityThreshold Used only if inputFilesType == "bam" or
#' "bamPairedEnd", otherwise ignored.
#' @param mappingQualityThreshold Used only if inputFilesType == "bam" or
#' "bamPairedEnd", otherwise ignored.
#' @param softclippingAllowed Used only if inputFilesType == "bam" or
#' "bamPairedEnd". Default is FALSE.
#' @param useMultiCore Logical indicating whether multiple cores are used (TRUE)
#' or not (FALSE). Default is FALSE. Only used for BAM files.
#' @param numCores Number of cores to use. Used only if useMultiCore = TRUE.
#' Default is NULL (auto-detect).
#' @param readsPerPiece Number of reads to process per chunk. Used only for BAM
#' files. Set to 'auto' for automatic memory-based detection. Default is 1000000.
#' @return Large List of elements - one element for each sample
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(exampleTSSr)
#' #getTSS(exampleTSSr)
#' }
setGeneric("getTSS",function(object
                             ,sequencingQualityThreshold = 10
                             ,mappingQualityThreshold = 20
                             ,softclippingAllowed = FALSE
                             ,useMultiCore = FALSE
                             ,numCores = NULL
                             ,readsPerPiece = 1000000)standardGeneric("getTSS"))
#' @rdname getTSS
#' @export
setMethod("getTSS",signature(object = "TSSr"), function(object
                                    ,sequencingQualityThreshold
                                    ,mappingQualityThreshold
                                    ,softclippingAllowed
                                    ,useMultiCore
                                    ,numCores
                                    ,readsPerPiece){
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabels <- object@sampleLabels
  inputFilesType <- object@inputFilesType
  if (length(object@sampleLabelsMerged) == 0) {
    object@sampleLabelsMerged <- sampleLabels
  }
  objName <- deparse(substitute(object))

  if(inputFilesType == "bam" | inputFilesType == "bamPairedEnd"){

    # Determine readsPerPiece
    if (readsPerPiece == 'auto') {
      readsPerPieceToUse <- .get_readsPerPiece()
    } else {
      readsPerPieceToUse <- readsPerPiece
    }

    if (useMultiCore) {
      # Multicore processing for BAM files
      if (is.null(numCores)) {
        numCores <- detectCores()
      }
      message("\nProcessing BAM files using ", numCores, " cores...")

      inputFilesID <- as.list(seq_len(length(object@inputFiles)))

      results <- mclapply(inputFilesID, function(x) {
        tssMC <- .getTSS_from_bam_chunked(
          bam.files = object@inputFiles[x],
          Genome = Genome,
          sampleLabels = sampleLabels[x],
          inputFilesType = inputFilesType,
          sequencingQualityThreshold = sequencingQualityThreshold,
          mappingQualityThreshold = mappingQualityThreshold,
          softclippingAllowed = softclippingAllowed,
          readsPerPiece = readsPerPieceToUse
        )
        return(tssMC)
      }, mc.cores = numCores)

      # Merge results from separate cores
      tss <- NULL
      for (i in results) {
        if (is.null(tss)) {
          tss <- i
        } else {
          tss <- merge(tss, i, by = c("chr", "pos", "strand"), all = TRUE)
        }
      }

      # Replace NA with 0 only in sample columns (not in chr, pos, strand!)
      if (ncol(tss) > 3) {
        sample_cols <- names(tss)[4:ncol(tss)]
        for (col in sample_cols) {
          set(tss, which(is.na(tss[[col]])), col, 0L)
        }
      }

    } else {
      # Single-core processing with chunked reading
      tss <- .getTSS_from_bam_chunked(
        bam.files = object@inputFiles,
        Genome = Genome,
        sampleLabels = sampleLabels,
        inputFilesType = inputFilesType,
        sequencingQualityThreshold = sequencingQualityThreshold,
        mappingQualityThreshold = mappingQualityThreshold,
        softclippingAllowed = softclippingAllowed,
        readsPerPiece = readsPerPieceToUse
      )
    }

  } else if(inputFilesType == "bed"){
    tss <- .getTSS_from_bed(object@inputFiles, Genome, sampleLabels)
  } else if(inputFilesType == "BigWig"){
    tss <- .getTSS_from_BigWig(object@inputFiles, Genome, sampleLabels)
  } else if(inputFilesType == "tss"){
    tss <- .getTSS_from_tss(object@inputFiles, sampleLabels)
  } else if(inputFilesType == "TSStable"){
    tss <- .getTSS_from_TSStable(object@inputFiles, sampleLabels)
  } else {
    stop("Unknown inputFilesType: '", inputFilesType, "'. ",
         "Supported types: bam, bamPairedEnd, bed, BigWig, tss, TSStable")
  }

  setorder(tss, "strand", "chr", "pos")
  # get library sizes
  object@librarySizes <- colSums(tss[, 4:ncol(tss), drop = FALSE], na.rm = TRUE)

  object@TSSrawMatrix <- tss
  object@TSSprocessedMatrix <- tss
  assign(objName, object, envir = parent.frame())
})


