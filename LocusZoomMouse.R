#locuszoomr is hard-coded to only consider human data.  These two functions
#fix that to work in mice (only, at least for now)
locus_mouse=function (gene = NULL, data = NULL, xrange = NULL, seqname = NULL, 
                      flank = NULL, fix_window = NULL, ens_db, chrom = NULL, pos = NULL, 
                      p = NULL, yvar = NULL, labs = NULL, index_snp = NULL, LD = NULL) 
{
  if (is.character(ens_db)) {
    if (!ens_db %in% (.packages())) {
      stop("Ensembl database not loaded. Try: library(", 
           ens_db, ")", call. = FALSE)
    }
    edb <- get(ens_db)
  }
  else edb <- ens_db
  if (!is.null(flank) & !is.null(fix_window)) 
    stop("both `flank` and `fix_window` cannot be specified at the same time")
  if (is.null(flank)) 
    flank <- 1e+05
  flank <- rep_len(flank, 2)
  if (!is.null(gene)) {
    locus <- genes(edb, filter = AnnotationFilterList(GeneNameFilter(gene), 
                                                      SeqNameFilter(c(1:22, "X", "Y"))))
    if (length(locus) > 1) {
      message(sprintf("Identified %d genes matching name '%s', taking first", 
                      length(locus), gene))
      locus <- locus[1]
    }
    seqname <- names(seqlengths(locus))
    if (is.null(fix_window)) {
      xrange <- c(start(locus) - flank[1], end(locus) + 
                    flank[2])
    }
    else {
      m <- mean(c(start(locus), end(locus)))
      xrange <- as.integer(c(m - fix_window/2, m + fix_window/2))
    }
    xrange[xrange < 0] <- 0
  }
  if (!is.null(data)) {
    dc <- detect_cols(data, chrom, pos, p, labs, yvar)
    chrom <- dc$chrom
    pos <- dc$pos
    p <- dc$p
    labs <- dc$labs
    if (!is.null(index_snp) & is.null(gene) & is.null(seqname) & 
        is.null(xrange)) {
      if (!index_snp %in% data[, labs]) 
        stop("SNP specified by `index_snp` not found")
      ind <- which(data[, labs] == index_snp)
      if (length(ind) > 1) 
        message("SNP found more than once")
      seqname <- data[ind[1], chrom]
      snp_pos <- data[ind[1], pos]
      xrange <- if (is.null(fix_window)) {
        c(snp_pos - flank[1], snp_pos + flank[2])
      }
      else {
        as.integer(c(snp_pos - fix_window/2, snp_pos + 
                       fix_window/2))
      }
      xrange[xrange < 0] <- 0
    }
  }
  if (is.null(xrange) | is.null(seqname)) 
    stop("No locus specified")
  msg <- paste0("chromosome ", seqname, ", position ", xrange[1], 
                " to ", xrange[2])
  if (!is.null(gene)) 
    msg <- paste(gene, msg, sep = ", ")
  if (!is.null(index_snp)) 
    msg <- paste(index_snp, msg, sep = ", ")
  message(msg)
  if (!is.null(data)) {
    data <- data[which(data[, chrom] == seqname), ]
    data <- data[which(data[, pos] > xrange[1] & data[, pos] < 
                         xrange[2]), ]
    data[data[, p] < 4.94065645841247e-324, p] <- 4.94065645841247e-324
    if (is.null(yvar)) {
      data$logP <- -log10(data[, p])
      yvar <- "logP"
    }
    data <- as.data.frame(data)
    if (nrow(data) == 0) {
      message("Locus contains no SNPs/datapoints")
      data <- NULL
    }
    else {
      message(nrow(data), " SNPs/datapoints")
      if (is.null(index_snp)) 
        index_snp <- data[which.max(data[, yvar]), labs]
      if (is.character(LD)) {
        colnames(data)[which(colnames(data) == LD)] <- "ld"
      }
    }
  }
  seqname <- gsub("chr|[[:punct:]]", "", seqname, ignore.case = TRUE)
  if (!seqname %in% c(1:22, "X", "Y")) 
    warning("`seqname` refers to a non-conventional chromosome")
  TX <- ensembldb::genes(edb, filter = AnnotationFilterList(SeqNameFilter(seqname), 
                                                            TxStartFilter(xrange[2], condition = "<"), TxEndFilter(xrange[1], 
                                                                                                                   condition = ">"), GeneIdFilter("ENSMUSG", "startsWith")))
  TX <- data.frame(TX)
  TX <- TX[!is.na(TX$start), ]
  TX <- TX[!duplicated(TX$gene_id), ]
  if (nrow(TX) == 0) {
    message("No gene transcripts")
    EX <- ensembldb::exons(edb, filter = AnnotationFilterList(SeqNameFilter(seqname), 
                                                              ExonStartFilter(xrange[2], condition = "<"), ExonEndFilter(xrange[1], 
                                                                                                                         condition = ">"), GeneIdFilter("ENSMUSG", "startsWith")))
  }
  else {
    EX <- ensembldb::exons(edb, filter = GeneIdFilter(TX$gene_id))
  }
  loc <- list(seqname = seqname, xrange = xrange, gene = gene, 
              ens_db = ens_db, chrom = chrom, pos = pos, p = p, yvar = yvar, 
              labs = labs, index_snp = index_snp, data = data, TX = TX, 
              EX = EX)
  class(loc) <- "locus"
  loc
}
detect_cols <- function(data, chrom, pos, p, labs = NULL, yvar = NULL) {
  # autodetect headings
  if (is.null(chrom)) {
    w <- grep("chr", colnames(data), ignore.case = TRUE)
    if (length(w) == 1) {
      chrom <- colnames(data)[w]
    } else stop("unable to autodetect chromosome column")
  }
  if (is.null(pos)) {
    w <- grep("pos", colnames(data), ignore.case = TRUE)
    if (length(w) == 1) {
      pos <- colnames(data)[w]
    } else stop("unable to autodetect SNP position column")
  }
  if (!is.null(p) && !is.null(yvar)) stop("cannot specify both `p` and `yvar`")
  if (is.null(p) && is.null(yvar)) {
    if ("p" %in% colnames(data)) {
      p <- "p"
    } else {
      w <- grep("^p?val", colnames(data), ignore.case = TRUE)
      if (length(w) == 1) {
        p <- colnames(data)[w]
      } else stop("unable to autodetect p-value column")
    }
  }
  if (is.null(labs)) {
    w <- grep("rs?id", colnames(data), ignore.case = TRUE)
    if (length(w) > 1) stop("unable to autodetect SNP id column")
    if (length(w) == 0) {
      w <- grep("SNP", colnames(data), ignore.case = TRUE)
    }
    if (length(w) == 1) {
      labs <- colnames(data)[w]
    } else stop("unable to autodetect SNP id column")
  }
  
  # check headings
  if (!chrom %in% colnames(data)) {
    stop("Column specified by `chrom` not found in `data`")}
  if (!pos %in% colnames(data)) {
    stop("Column specified by `pos` not found in `data`")}
  if (is.null(yvar) && !p %in% colnames(data)) {
    stop("Column specified by `p` not found in `data`")}
  if (!labs %in% colnames(data)) {
    stop("Column specified by `labs` not found in `data`")}
  if (!is.null(yvar)) {
    if (!yvar %in% colnames(data)) {
      stop("Column specified by `yvar` not found in `data`")
    }
  }
  
  list(chrom = chrom, pos = pos, p = p, labs = labs)
}