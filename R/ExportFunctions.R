###############################################################################
.plotCorrelation <- function(TSS.all.samples, maxPoints = 50000)
{
  z_full <- TSS.all.samples[,-c(1,2,3)]

  # Get sample names and calculate appropriate font size
  sample_names <- colnames(z_full)
  max_name_length <- max(nchar(sample_names))
  n_samples <- ncol(z_full)
  n_points <- nrow(z_full)

  # Pre-calculate correlation matrix using ALL data (before sampling)
  cor_matrix <- cor(z_full, use = "pairwise.complete.obs")

  # Apply sampling if needed for visualization
  if (!is.null(maxPoints) && n_points > maxPoints) {
    message(paste("Sampling", maxPoints, "points from", n_points, "for visualization..."))
    message("(Correlation coefficients are calculated using ALL data)")
    set.seed(42)  # For reproducibility
    sample_idx <- sample(n_points, maxPoints)
    z <- z_full[sample_idx, ]
  } else {
    z <- z_full
  }

  # Calculate font size based on longest name and number of samples
  # Base size decreases with more samples and longer names
  base_cex <- min(1.2, 8 / max_name_length, 6 / n_samples)
  label_cex <- max(0.4, base_cex)  # Minimum size of 0.4

  # Customize lower panel (scatter plot with sampled data)
  panel.scatter <- function(x, y)
  {
    points(x, y, pch = ".", col = "#00AFBB")
  }

  # Customize upper panel (correlation coefficient from FULL data)
  panel.cor <- function(x, y)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)

    # Get variable indices to look up pre-calculated correlation
    x_name <- deparse(substitute(x))
    y_name <- deparse(substitute(y))

    # Find matching column names
    x_idx <- which(sapply(seq_along(z), function(i) identical(x, z[[i]])))
    y_idx <- which(sapply(seq_along(z), function(i) identical(y, z[[i]])))

    if (length(x_idx) > 0 && length(y_idx) > 0) {
      # Use pre-calculated correlation from full dataset
      r <- round(cor_matrix[x_idx[1], y_idx[1]], digits = 2)
    } else {
      # Fallback: calculate from current data
      r <- round(cor(x, y), digits = 2)
    }

    txt <- paste0(r)
    cex.cor <- 0.8 / strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  }

  # Set margins: increase outer margin for long labels
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # Adjust outer margins based on label length
  oma_size <- max(1, min(4, max_name_length / 8))
  par(oma = c(oma_size, oma_size, oma_size, oma_size))

  # Create the plots with auto-scaled labels
  pairs(z,
        lower.panel = panel.scatter,
        upper.panel = panel.cor,
        diag.panel = NULL,  # Use text.panel instead
        log = "xy",
        cex.labels = label_cex,
        font.labels = 2,
        gap = 0.3)
}

###############################################################################
###############################################################################
##.plotTSS function plots TSS graph
##.plotTSS function takes three input files, tss.tpm, cs, and ref
##tss.tpm table has at least 4 columns (chr, pos, strand, tpm). tpm value is negative in the minus strand
##cs table has 11 columns (cluster,chr,start,end,strand,dominant_tss,tpm,tpm.dominant_tss,q_0.1,q_0.9,interquantile_width)
##ref table has at least 5 columns (gene,chr, start, end, strand)
##run script with the following example command:
##.plotTSS(tss.tpm,cs.cl, ref, up.dis = 500, down.dis=100)
.plotTSS <- function(tss, cs,df, samples, Bidirection, up.dis, down.dis,yFixed)
{
  setnames(df, colnames(df)[c(1,6)], c("chr","gene"))

  if(df$strand == "+")
  {
    p <- df$start - up.dis
    q <- df$end + down.dis
  }
  else
  {
    p <- df$start - down.dis
    q <- df$end + up.dis
  }
  ##Genome range (genomic coordinates) track
  range.gr <- GRanges(seqnames = as.character(df$chr),
                      ranges = IRanges(start = p,end = q),
                      strand = as.character(df$strand))
  gtrack <- GenomeAxisTrack(range = range.gr,fontcolor = "black")
  ##Gene region track
  gene.gr <- makeGRangesFromDataFrame(df,keep.extra.columns=FALSE,ignore.strand=FALSE,
                                      seqinfo=NULL,seqnames.field=c("chr"),
                                      start.field="start",end.field=c("end"),
                                      strand.field="strand",starts.in.df.are.0based=FALSE)
  atrack.gene <- GeneRegionTrack(gene.gr, name = "gene",col.title="black",
                                 transcriptAnnotation = "gene")
  ##TSS cluster track and TSS Data track
  if(Bidirection == TRUE)
  {
    tss_sub <- tss[tss$chr == as.character(df$chr) & tss$pos >= p & tss$pos <= q,]
    tss_sub[, strand:= df$strand]
  }
  else
  {
    tss_sub <- tss[tss$chr == as.character(df$chr) & tss$strand == as.character(df$strand) & tss$pos >= p & tss$pos <= q,]
  }

  # setup y_range
  # To hide TSS signals of the opposite strand, change "min(tss_sub1)" to "0" for positive strand and change "max(tss_sub1)" to "0" for negative strand
  tss_sub1 <- tss_sub[,.SD, .SDcols = c(samples)]
  if(yFixed==TRUE)
  {
    if(tss_sub$strand[1] == "+")
    {
    y_range=c(min(tss_sub1),max(tss_sub1))
    }
    else
    {
      y_range=c(min(tss_sub1),max(tss_sub1))
     }
  }
  else
  {
    y_range=NULL
  }

  data_tss_track <- list()
  dtrack <- list()
  s <- c()
  for (my.sample in seq(samples))
  {
    ##TSS cluster track
    cs.temp <- cs[[samples[my.sample]]]
    cs_sub <- cs.temp[cs.temp$chr == as.character(df$chr) & cs.temp$strand == as.character(df$strand) & cs.temp$q_0.1 >= p & cs.temp$q_0.9 <= q,]
    cs.gr <- makeGRangesFromDataFrame(cs_sub,keep.extra.columns=FALSE,
                                      ignore.strand=FALSE,seqinfo=NULL,
                                      seqnames.field=c("chr"),
                                      start.field="q_0.1",end.field=c("q_0.9"),
                                      strand.field="strand",
                                      starts.in.df.are.0based=FALSE)
    # atrack.cs <- AnnotationTrack(cs.gr, name = paste(samples[my.sample],"clusters",sep = " "),col.title="black")
    atrack.cs <- AnnotationTrack(cs.gr, name = paste(samples[my.sample],"clusters",sep = " "),
                                 id = cs_sub$cluster, col.title="black")
    ##TSS track
    temp <- tss_sub[,.SD, .SDcols = c("chr","pos","strand",samples[my.sample])]
    data_tss_track <- makeGRangesFromDataFrame(temp,keep.extra.columns=TRUE,
                                               ignore.strand=FALSE,seqinfo=NULL,
                                               seqnames.field=c("chr"),
                                               start.field="pos",end.field=c("pos"),
                                               strand.field="strand",starts.in.df.are.0based=FALSE)
    dtrack <- DataTrack(data_tss_track, name = paste(samples[my.sample],"TSS (TPM)", sep = "\n"),
                        type = "h",col = rainbow(length(samples))[my.sample],
                        baseline = 0,col.baseline = "grey",
                        col.title="black",col.axis = "black",ylim = y_range)
    s <- c(s, atrack.cs, dtrack)
  }
  ##plot Genome range track, gene track, clusters track, TSS track
  plotTracks(c(list(gtrack, atrack.gene),s),
             main = df$gene,
             featureAnnotation="id",
             fontcolor.feature="black",
             cex.feature=0.7, shape = "arrow")
}
###############################################################################
.getBed <- function(p){
  ##define variable as a NULL value
  chr = q_0.1 = q_0.9 = dominant_tss = NULL
  pbed <- lapply(as.list(seq_len(p[,.N])), function(i){
    if(p[i,q_0.1] == p[i, q_0.9]){
      nrBlocks = length(unique(c(p[i,start],p[i,end],p[i,q_0.1],p[i,q_0.9])))
    }else{nrBlocks = length(unique(c(p[i,start],p[i,end],p[i,q_0.1],p[i,q_0.9]))) -1}

    if(nrBlocks == 3){
      block.sizes = paste(1, p[i,q_0.9]-p[i,q_0.1]+1, 1, sep = ",")
      block.starts = paste(0, p[i,q_0.1]-p[i,start], p[i,end]-p[i,start], sep = ",")
    }else if(nrBlocks ==2){
      if(p[i,q_0.1] == p[i, start]){
        block.sizes = paste(p[i,q_0.9]-p[i,q_0.1]+1,1, sep = ",")
        block.starts = paste(0, p[i,end]-p[i,start], sep = ",")
      }else if(p[i,end] == p[i, q_0.9]){
        block.sizes = paste(1, p[i,q_0.9]-p[i,q_0.1]+1, sep = ",")
        block.starts = paste(0, p[i,q_0.1]-p[i,start], sep = ",")
      }else{
        message("\nWhat is the additional condition which has two blocks...")
        print(p[i,])
      }
    }else{
      block.sizes = paste(p[i,end]-p[i,start])
      block.starts = paste(0)
    }
    list(p[i,chr]
         ,p[i,start]-1
         ,p[i,end]
         ,"."
         ,0
         ,p[i,strand]
         ,p[i,dominant_tss]-1
         ,p[i,dominant_tss]
         ,0
         ,nrBlocks
         ,block.sizes
         ,block.starts)
  })
  pbed <- rbindlist(pbed)
  return(pbed)
}
