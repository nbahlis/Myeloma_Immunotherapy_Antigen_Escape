###############
#### Utils ####
###############

get.cnv.data <- function(file, n, barcodes){
  require(BSgenome.Hsapiens.UCSC.hg38)
  
  chr.lengths = seqlengths(Hsapiens)[1:22]
  chr.names = names(chr.lengths)
  
  cnv.data <- rhdf5::H5Fopen(file)
  
  windows  <- GenomicRanges::tileGenome(chr.lengths, tilewidth = cnv.data$constants$bin_size, cut.last.tile.in.chrom = TRUE)
  windows <- split(windows, GenomicRanges::seqnames(windows))
  
  mtx <- do.call(rbind, sapply(chr.names, function(x){ print(x); as.matrix( get.chr.cnv.data(file, x, chr.lengths=chr.lengths, windows=windows, n=n, barcodes=barcodes) ) }))
  return(mtx)
}

get.chr.cnv.data <- function(file, chr, chr.lengths, windows, n, barcodes, max.n.copies=13){
  
  cnv.data <- rhdf5::H5Fopen(file)
  
  windows  <- GenomicRanges::tileGenome(chr.lengths, tilewidth = cnv.data$constants$bin_size, cut.last.tile.in.chrom = TRUE)
  windows <- split(windows, GenomicRanges::seqnames(windows))
  
  cnv.matrix.chr <- cnv.data$cnvs[[chr]][, 1:cnv.data$constants$num_cells]
  
  rownames(cnv.matrix.chr) <- as.character( windows[[chr]] )
  colnames(cnv.matrix.chr) <- cnv.data$cell_barcodes
  # see https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/output/hdf5
  cnv.matrix.chr[cnv.matrix.chr == -128] <- NA
  cnv.matrix.chr[cnv.matrix.chr == -127] <- 0
  cnv.matrix.chr <- abs(cnv.matrix.chr)
  cnv.matrix.chr[ cnv.matrix.chr > max.n.copies ] <- max.n.copies
  
  cnv.matrix.chr <- cnv.matrix.chr[, barcodes]
  
  # thanks to https://stackoverflow.com/questions/30359427/calculate-the-mean-of-every-13-rows-in-data-frame
  aggregated.cnv.matrix <- aggregate(cnv.matrix.chr, list( rep(1:(nrow(cnv.matrix.chr)%/%n+1), each=n, len=nrow(cnv.matrix.chr))), mean, na.rm=T)[-1];
  
  windows.2  <- GenomicRanges::tileGenome(chr.lengths[chr], tilewidth = cnv.data$constants$bin_size*n, cut.last.tile.in.chrom = TRUE)
  windows.2  <- split(windows.2, GenomicRanges::seqnames(windows.2))
  
  rownames(aggregated.cnv.matrix) <- as.character(windows.2[[chr]])
  return(aggregated.cnv.matrix)
}
