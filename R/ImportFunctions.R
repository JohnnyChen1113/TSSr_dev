################################################################################
.getGenome <- function(genomeName) {
  if (is.null(genomeName) || genomeName == "") {
    stop("Error: 'genomeName' must not be NULL or empty.")
  }
  
  if (!genomeName %in% rownames(installed.packages())) {
    stop(paste("Error: Requested genome package '", genomeName, "' is not installed. Please install it before running this function.", sep=""))
  }
  
  if (!requireNamespace(genomeName, quietly = TRUE)) {
    stop(paste("Error: Namespace for '", genomeName, "' could not be loaded. Check if the package is correctly installed.", sep=""))
  }
  
  getExportedValue(genomeName, genomeName)
}

################################################################################
##.getTSS_from_bam function calls TSS from bam files
.getTSS_from_bam <- function(bam.files, Genome, sampleLabels, inputFilesType
                             ,sequencingQualityThreshold
                             ,mappingQualityThreshold
                             ,softclippingAllowed){
  ## check the bam.files is empty or not
  if (length(bam.files) == 0) {
    stop("Error: 'bam.files' is empty.")
  }
  if (!exists("Genome")) {
    stop("Error: 'Genome' is not defined.")
  }
  if (length(sampleLabels) != length(bam.files)) {
    stop("Error: Length of 'sampleLabels' must match the number of 'bam.files'.")
  }
  
  ##define variable as a NULL value
  chr = pos = tag_count = strand = NULL

  what <- c("rname", "strand", "pos", "seq", "qual", "mapq","flag","cigar")
  param <- ScanBamParam( what = what
                         , flag = scanBamFlag(isUnmappedQuery = FALSE,
                                              isNotPassingQualityControls = FALSE)
                         , mapqFilter = mappingQualityThreshold)
  if (inputFilesType == "bamPairedEnd"){
    Rsamtools::bamFlag(param) <- scanBamFlag( isUnmappedQuery = FALSE
                                   , isProperPair    = TRUE
                                   , isFirstMateRead = TRUE)}
  chunksize <- 1e6
  first <- TRUE
  for(i in seq_len(length(bam.files))) {
    message("\nReading in file: ", bam.files[i], "...")
    bam <- scanBam(bam.files[i], param = param)
    message("\t-> Filtering out low quality reads...")
    qual <- bam[[1]]$qual
    start <- 1
    # chunksize <- 1e6
num_chunks <- ceiling(length(qual) / chunksize)
qa.avg <- vector(mode = "integer", length = num_chunks)
chunk_index <- 1
while (start <= length(qual)) {
  end <- min(start + chunksize - 1, length(qual))
  qa.avg[chunk_index] <- as.integer(mean(as.integer(qual[start:end])))
  start <- end + 1
  chunk_index <- chunk_index + 1
}

    cigar <- bam[[1]]$cigar
    start <- 1
    # chunksize <- 1e6
## update 20240605
# Initialize the vector for storing mapped lengths
mapped.length <- integer(0)

# Pre-process the CIGAR strings to extract all numeric values and identify soft-clipped segments
# This avoids repeated string operations inside the loop
numeric_parts <- lapply(bam[[1]]$cigar, function(cigar) {
  as.integer(str_extract_all(cigar, "[0-9]+")[[1]])
})
soft_clipped_lengths <- if (softclippingAllowed) {
  sapply(bam[[1]]$cigar, function(cigar) {
    sc_part <- str_extract(cigar, "[0-9]+S")
    if (is.na(sc_part)) {
      0
    } else {
      as.integer(sub("S", "", sc_part))
    }
  })
} else {
  integer(length(bam[[1]]$cigar))
}

# Process in chunks to handle large arrays more efficiently
start <- 1
while (start <= length(cigar)) {
  end <- min(start + chunksize - 1, length(cigar))
  
  # Calculate the total mapped length by summing numeric parts and adjusting for soft clipping if allowed
  chunk_numeric_parts <- numeric_parts[start:end]
  chunk_soft_clipped_lengths <- soft_clipped_lengths[start:end]
  
  total_lengths <- vapply(chunk_numeric_parts, sum, numeric(1))
  if (softclippingAllowed) {
    total_lengths <- total_lengths - chunk_soft_clipped_lengths
  }
  
  # Append the computed lengths for this chunk to the overall mapped lengths vector
  mapped.length <- c(mapped.length, total_lengths)
  
  if (end == length(cigar)) {
    break
  } else {
    start <- end + 1
  }
}

    
readsGR <- GRanges(
  seqnames = as.vector(bam[[1]]$rname),
  IRanges(start = bam[[1]]$pos, width = mapped.length),
  strand = bam[[1]]$strand,
  qual = qa.avg,
  mapq = bam[[1]]$mapq,
  seq = bam[[1]]$seq,
  read.length = width(bam[[1]]$seq),
  flag = bam[[1]]$flag
)

# 过滤掉在基因组中不存在的染色体
valid_chromosomes <- as.character(seqnames(Genome))
readsGR <- readsGR[seqnames(readsGR) %in% valid_chromosomes]

# 过滤掉超过基因组长度的读取
valid_seqnames <- as.character(seqnames(readsGR))
seq_lengths <- seqlengths(Genome)[valid_seqnames]
readsGR <- readsGR[end(readsGR) <= seq_lengths[valid_seqnames]]

# 将 NA 的 mapq 值替换为 Inf
elementMetadata(readsGR)$mapq[is.na(elementMetadata(readsGR)$mapq)] <- Inf

# 过滤正链读取
readsGR.p <- readsGR[
  strand(readsGR) == "+" &
  elementMetadata(readsGR)$qual >= sequencingQualityThreshold &
  elementMetadata(readsGR)$mapq >= mappingQualityThreshold
]

# 过滤负链读取
readsGR.m <- readsGR[
  strand(readsGR) == "-" &
  elementMetadata(readsGR)$qual >= sequencingQualityThreshold &
  elementMetadata(readsGR)$mapq >= mappingQualityThreshold
]

    
    if(softclippingAllowed){
      ##------------------------------------------------------------------------
      TSS.p <- data.table(chr = as.character(seqnames(readsGR.p)), 
                          pos = start(readsGR.p), strand = "+", 
                          stringsAsFactors = FALSE)
      ##------------------------------------------------------------------------
      TSS.m <- data.table(chr = as.character(seqnames(readsGR.m)), 
                          pos = end(readsGR.m), strand = "-", 
                          stringsAsFactors = FALSE)
      #-------------------------------------------------------------------------
      TSS <- rbind(TSS.p, TSS.m)
      TSS <- TSS[,c("chr", "pos", "strand")]
      TSS$tag_count <- 1
      setDT(TSS)
      TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
      
    }else{
      # remove G mismatch
      TSS <- .removeNewG(readsGR.p, readsGR.m, Genome)
    }

    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
    gc()
  }
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}

################################################################################
.removeNewG <- function(readsGR.p, readsGR.m, Genome) {
  # Define variables as NULL values
  chr = pos = tag_count = Gp = Gm = i = NULL

  message("\t-> Removing the bases of the reads if mismatched 'Gs'...")

  # Optimized processing for the plus strand
  processPlusStrand <- function(reads, Genome) {
    Gp <- which(substr(GenomicRanges::elementMetadata(reads)$seq, start = 1, stop = 1) == "G")
    i <- 1
    
    while (length(Gp) > 0) {
      # Identify mismatched G bases
      G.mismatch <- Gp[getSeq(Genome, GenomicRanges::resize(reads[Gp], width = 1, fix = "start"), as.character = TRUE) != "G"]
      # Adjust start position of reads with mismatched G
      start(reads[G.mismatch]) <- start(reads[G.mismatch]) + as.integer(1)
      i <- i + 1
      # Update Gp for the next iteration
      Gp <- G.mismatch[which(substr(GenomicRanges::elementMetadata(reads[G.mismatch])$seq, start = 1, stop = i) == paste(rep("G", i), collapse = ""))]
    }
    
    # Create data table for plus strand TSS
    data.table(chr = as.character(seqnames(reads)), pos = start(reads), strand = "+", stringsAsFactors = FALSE)
  }

  # Optimized processing for the minus strand
  processMinusStrand <- function(reads, Genome) {
    Gm <- which(substr(GenomicRanges::elementMetadata(reads)$seq, start = GenomicRanges::elementMetadata(reads)$read.length, stop = GenomicRanges::elementMetadata(reads)$read.length) == "C")
    i <- 1
    
    while (length(Gm) > 0) {
      # Identify mismatched C bases
      G.mismatch <- Gm[getSeq(Genome, GenomicRanges::resize(reads[Gm], width = 1, fix = "start"), as.character = TRUE) != "G"]
      # Adjust end position of reads with mismatched C
      end(reads[G.mismatch]) <- end(reads[G.mismatch]) - as.integer(1)
      i <- i + 1
      # Update Gm for the next iteration
      Gm <- G.mismatch[which(substr(GenomicRanges::elementMetadata(reads[G.mismatch])$seq, start = 1, stop = i) == paste(rep("C", i), collapse = ""))]
    }
    
    # Create data table for minus strand TSS
    data.table(chr = as.character(seqnames(reads)), pos = end(reads), strand = "-", stringsAsFactors = FALSE)
  }

  # Process plus and minus strands
  TSS.p <- processPlusStrand(readsGR.p, Genome)
  TSS.m <- processMinusStrand(readsGR.m, Genome)
  
  # Combine plus and minus strand TSS and summarize
  TSS <- rbind(TSS.p, TSS.m)
  TSS <- TSS[, .(tag_count = .N), by = .(chr, pos, strand)]
  
  return(TSS)
}


################################################################################
##.getTSS_from_bed function calls TSS from bed files
.getTSS_from_bed <- function(bed.files, Genome, sampleLabels){
  first <- TRUE
  ##define variable as a NULL value
  chr = pos = tag_count = NULL

  for(i in seq_len(length(bed.files))) {
    message("\nReading in file: ", bed.files[i], "...")
    readsGR <- import(bed.files[i], format = "BED")
    readsGR <- readsGR[as.character(seqnames(readsGR)) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    readsGR.p <- readsGR[strand(readsGR) == "+"]
    readsGR.m <- readsGR[strand(readsGR) == "-"]
    message("\t-> Making TSS table...")
    TSS.plus <- data.table(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), stringsAsFactors = FALSE)
    TSS.minus <- data.table(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), stringsAsFactors = FALSE)
    TSS <- rbind(TSS.plus, TSS.minus)
    TSS$tag_count <- 1
    TSS <- data.table(TSS)
    TSS <- TSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
    gc()
  }
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}
################################################################################
##.getTSS_from_BigWig function calls TSS from BigWig files

.getTSS_from_BigWig <- function(BigWig.files, Genome, sampleLabels){
  #library.sizes <- vector()
  ##define variable as a NULL value
  chr = pos = NULL
  first <- TRUE
  for(i in seq_len(length(BigWig.files))) {
    message("\nReading in file: ", BigWig.files[i], "...")
    readsGR <- import(BigWig.files[i], format = "BigWig")
    readsGR <- readsGR[as.character(seqnames(readsGR)) %in% seqnames(Genome)]
    readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
    readsGR.p <- readsGR[score(readsGR) > 0]
    readsGR.m <- readsGR[score(readsGR) < 0]
    message("\t-> Making TSS table...")
    TSS.plus <- data.table(chr = as.character(seqnames(readsGR.p)), pos = as.integer(start(readsGR.p)), strand = rep("+", times = length(readsGR.p)), score = as.numeric(abs(readsGR.p$score)), stringsAsFactors = FALSE)
    TSS.minus <- data.table(chr = as.character(seqnames(readsGR.m)), pos = as.integer(end(readsGR.m)), strand = rep("-", times = length(readsGR.m)), score = as.numeric(abs(readsGR.m$score)), stringsAsFactors = FALSE)
    TSS <- rbind(TSS.plus, TSS.minus)

    setDT(TSS)

    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)

    #library.sizes <- c(library.sizes, as.integer(sum(data.table(TSS)[,4])))
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
    gc()
  }
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}


################################################################################
##.getTSS_from_tss function calls TSS from tss files

.getTSS_from_tss <- function(tss.files, sampleLabels){
  first <- TRUE

  for(i in seq_len(length(tss.files))) {
    message("\nReading in file: ", tss.files[i], "...")
    TSS <- read.table(file = tss.files[i], header = TRUE, sep = "\t"
                      ,colClasses = c("character", "integer", "character", "integer")
                      ,col.names = c("chr", "pos", "strand", sampleLabels[i]))

    setDT(TSS)

    setkeyv(TSS, cols = c("chr", "pos", "strand"))
    if(first == TRUE) {
      TSS.all.samples <- TSS
    }else{
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
    gc()
  }
  TSS.all.samples <- data.table(TSS.all.samples)
  TSS.all.samples[,4:ncol(TSS.all.samples)][is.na(TSS.all.samples[,4:ncol(TSS.all.samples)])] =0
  return(TSS.all.samples)
}

################################################################################################
##.getTSS_from_TSStable function calls TSS from one TSStable file

.getTSS_from_TSStable <- function(TSStable.file, sampleLabels){
  # Check if only one file is provided
  if(length(TSStable.file) > 1){
    stop("Only one file should be provided when inputFilesType = 'TSStable', but provided ", length(TSStable.file), " files.")
  }
  
  # Check if file exists
  if(!file.exists(TSStable.file)){
    stop("Could not locate input file: ", normalizePath(TSStable.file, mustWork = FALSE))
  }

  # Read the file using fread for better performance and directly create a data.table
  TSS.all.samples <- fread(file = TSStable.file, header = TRUE, colClasses = c("character", "integer", "character", rep("integer", length(sampleLabels))))

  # Verify the number of columns matches the expected number
  expected_cols <- length(sampleLabels) + 3
  if(ncol(TSS.all.samples) != expected_cols){
    stop("Number of provided sample labels must match the number of samples in the TSS table! Expected columns: ", expected_cols, ", but found ", ncol(TSS.all.samples), " columns.")
  }

  # Replace NA values with 0 in sample columns
  TSS.all.samples[, (sampleLabels) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = sampleLabels]

  return(TSS.all.samples)
}

################################################################################
#' Process BAM chunk to TSS
#'
#' @description Internal function to convert a BAM chunk to TSS data.
#' Used for chunked reading of large BAM files.
#'
#' @param bam BAM data from scanBam
#' @param Genome BSgenome object
#' @param sequencingQualityThreshold Quality threshold
#' @param mappingQualityThreshold Mapping quality threshold
#' @param softclippingAllowed Whether softclipping is allowed
#' @param chunksize Chunk size for processing
#' @return A data.table with TSS positions and tag counts
#' @keywords internal
.bam_to_TSS <- function(bam, Genome,
                        sequencingQualityThreshold,
                        mappingQualityThreshold,
                        softclippingAllowed,
                        chunksize) {

  # Calculate average quality scores in chunks
  qual <- bam[[1]]$qual
  start <- 1
  qa.avg <- vector(mode = "integer")

  repeat {
    if (start + chunksize <= length(qual)) {
      end <- start + chunksize
    } else {
      end <- length(qual)
    }
    qa.avg <- c(qa.avg, as.integer(sapply(as(qual[start:end], "IntegerList"), mean)))
    if (end == length(qual)) {
      break
    } else {
      start <- end + 1
    }
  }
  qa.avg[is.na(qa.avg)] <- -1

  # Calculate mapped lengths from CIGAR strings
  cigar <- bam[[1]]$cigar
  start <- 1
  mapped.length <- vector(mode = "integer")

  repeat {
    if (start + chunksize <= length(cigar)) {
      end <- start + chunksize
    } else {
      end <- length(cigar)
    }
    if (softclippingAllowed) {
      mapped.length <- c(mapped.length,
                         as.integer(sum(as(str_extract_all(bam[[1]]$cigar[start:end], "([0-9]+)"), "IntegerList"))) -
                           ifelse(is.na(sub("S", "", str_extract(bam[[1]]$cigar[start:end], "[0-9]+S"))),
                                  0,
                                  as.integer(sub("S", "", str_extract(bam[[1]]$cigar[start:end], "[0-9]+S")))))
    } else {
      mapped.length <- c(mapped.length,
                         as.integer(sum(as(str_extract_all(bam[[1]]$cigar[start:end], "([0-9]+)"), "IntegerList"))))
    }
    if (end == length(cigar)) {
      break
    } else {
      start <- end + 1
    }
  }

  # Create GRanges object
  readsGR <- GRanges(
    seqnames = as.vector(bam[[1]]$rname),
    IRanges(start = bam[[1]]$pos, width = mapped.length),
    strand = bam[[1]]$strand,
    qual = qa.avg,
    mapq = bam[[1]]$mapq,
    seq = bam[[1]]$seq,
    read.length = width(bam[[1]]$seq),
    flag = bam[[1]]$flag
  )

  # Filter by genome
  readsGR <- readsGR[as.character(readsGR@seqnames) %in% seqnames(Genome)]
  readsGR <- readsGR[!(end(readsGR) > seqlengths(Genome)[as.character(seqnames(readsGR))])]
  GenomicRanges::elementMetadata(readsGR)$mapq[is.na(GenomicRanges::elementMetadata(readsGR)$mapq)] <- Inf

  # Filter by quality
  readsGR.p <- readsGR[(as.character(strand(readsGR)) == "+" &
                          GenomicRanges::elementMetadata(readsGR)$qual >= sequencingQualityThreshold) &
                          GenomicRanges::elementMetadata(readsGR)$mapq >= mappingQualityThreshold]
  readsGR.m <- readsGR[(as.character(strand(readsGR)) == "-" &
                          GenomicRanges::elementMetadata(readsGR)$qual >= sequencingQualityThreshold) &
                          GenomicRanges::elementMetadata(readsGR)$mapq >= mappingQualityThreshold]

  if (softclippingAllowed) {
    TSS.p <- data.table(chr = as.character(seqnames(readsGR.p)),
                        pos = start(readsGR.p),
                        strand = "+",
                        stringsAsFactors = FALSE)
    TSS.m <- data.table(chr = as.character(seqnames(readsGR.m)),
                        pos = end(readsGR.m),
                        strand = "-",
                        stringsAsFactors = FALSE)
    TSS <- rbind(TSS.p, TSS.m)
    TSS <- TSS[, c("chr", "pos", "strand")]
    TSS$tag_count <- 1
    setDT(TSS)
    TSS <- TSS[, .(tag_count = as.integer(sum(tag_count))), by = list(chr, pos, strand)]
  } else {
    # Remove G mismatch
    TSS <- .removeNewG(readsGR.p, readsGR.m, Genome)
  }

  return(TSS)
}

################################################################################
#' Auto-detect optimal readsPerPiece based on available memory
#'
#' @description Detects available system memory and calculates an optimal
#' chunk size for reading BAM files.
#'
#' @return Numeric value for readsPerPiece parameter
#' @keywords internal
.get_readsPerPiece <- function() {

  sysInfo <- Sys.info()

  if (!is.null(sysInfo)) {
    OS <- sysInfo['sysname']
  } else {
    OS <- "Unknown"
  }

  available_memory_gb <- NULL

  # Windows
  if (OS == "Windows") {
    tryCatch({
      cmd <- "wmic OS get FreePhysicalMemory /Value"
      output <- system2("cmd", args = c("/c", cmd), stdout = TRUE)
      available_memory <- as.numeric(gsub("\\D", "", output))
      available_memory <- na.omit(available_memory)
      available_memory_gb <- available_memory[1] / (1024^2)
      message(available_memory_gb, " GB of memory available on your Windows machine")
    }, error = function(e) {
      message("Could not detect memory on Windows, using default value")
      available_memory_gb <<- 4
    })
  }

  # Linux/Mac
  if (OS != "Windows") {
    tryCatch({
      cmd <- "free -b | grep 'Mem:' | awk '{print $7}'"
      available_memory <- as.numeric(system2("bash", args = c("-c", cmd), stdout = TRUE))
      available_memory <- na.omit(available_memory)
      available_memory_gb <- available_memory[1] / (1024^3)
      message(round(available_memory_gb, 2), " GB of memory available on your machine")
    }, error = function(e) {
      message("Could not detect memory, using default value")
      available_memory_gb <<- 4
    })
  }

  if (is.null(available_memory_gb) || is.na(available_memory_gb)) {
    available_memory_gb <- 4
  }

  # Calculate readsPerPiece: roughly 0.3GB per 1M reads
  readsPerPieceLocal <- as.integer(available_memory_gb * 0.8 / 0.3 * 1000000)
  readsPerPieceLocal <- max(readsPerPieceLocal, 1000000)  # Minimum 1M reads

  message("readsPerPiece set to: ", format(readsPerPieceLocal, big.mark = ",", scientific = FALSE))

  return(readsPerPieceLocal)
}

################################################################################
#' Get TSS from BAM files with chunked reading support
#'
#' @description Modified version of .getTSS_from_bam that supports chunked
#' reading for large BAM files to reduce memory usage.
#'
#' @param bam.files Vector of BAM file paths
#' @param Genome BSgenome object
#' @param sampleLabels Sample labels
#' @param inputFilesType Input file type
#' @param sequencingQualityThreshold Quality threshold
#' @param mappingQualityThreshold Mapping quality threshold
#' @param softclippingAllowed Whether softclipping is allowed
#' @param readsPerPiece Number of reads to process at a time
#' @return A data.table with TSS data for all samples
#' @keywords internal
.getTSS_from_bam_chunked <- function(bam.files, Genome, sampleLabels, inputFilesType,
                                      sequencingQualityThreshold,
                                      mappingQualityThreshold,
                                      softclippingAllowed,
                                      readsPerPiece) {
  ##define variable as a NULL value
  chr = pos = tag_count = strand = NULL

  what <- c("rname", "strand", "pos", "seq", "qual", "mapq", "flag", "cigar")
  param <- ScanBamParam(what = what,
                        flag = scanBamFlag(isUnmappedQuery = FALSE,
                                           isNotPassingQualityControls = FALSE),
                        mapqFilter = mappingQualityThreshold)
  if (inputFilesType == "bamPairedEnd") {
    Rsamtools::bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE,
                                              isProperPair = TRUE,
                                              isFirstMateRead = TRUE)
  }

  chunksize <- 1e6
  first <- TRUE

  for (i in seq_len(length(bam.files))) {
    message("\nReading in file: ", bam.files[i], "...")

    bamFile <- BamFile(bam.files[i])
    yieldSize(bamFile) <- readsPerPiece

    records <- countBam(bamFile)
    n_chunks <- as.integer(records$records / readsPerPiece) + 1
    message("\t-> Processing in ", n_chunks, " chunk(s)...")

    open(bamFile)
    TSS <- NULL

    for (j in seq_len(n_chunks)) {
      bami <- scanBam(bamFile, param = param)

      # Skip if no reads
      if (length(bami[[1]]$pos) == 0) next

      tssi <- .bam_to_TSS(bami, Genome,
                          sequencingQualityThreshold,
                          mappingQualityThreshold,
                          softclippingAllowed,
                          chunksize)

      if (is.null(TSS)) {
        TSS <- tssi
      } else {
        TSS <- rbind(TSS, tssi)
      }
    }

    close(bamFile)

    if (is.null(TSS) || nrow(TSS) == 0) {
      warning("No TSS data obtained from file: ", bam.files[i])
      next
    }

    # Aggregate tag counts
    TSS <- TSS[, .(V1 = as.integer(sum(tag_count))), by = list(chr, pos, strand)]
    setnames(TSS, c("chr", "pos", "strand", sampleLabels[i]))
    setkey(TSS, chr, pos, strand)

    if (first == TRUE) {
      TSS.all.samples <- TSS
    } else {
      TSS.all.samples <- merge(TSS.all.samples, TSS, all = TRUE)
    }
    first <- FALSE
    gc()
  }

  # Replace NA with 0 only in sample columns (not in chr, pos, strand!)
  sample_cols <- names(TSS.all.samples)[4:ncol(TSS.all.samples)]
  for (col in sample_cols) {
    set(TSS.all.samples, which(is.na(TSS.all.samples[[col]])), col, 0L)
  }

  return(TSS.all.samples)
}
