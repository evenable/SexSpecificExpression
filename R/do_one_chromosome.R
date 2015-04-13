#' Annotate Genome with SNPs for a Single Chromosome
#' 
#' This function annotates the genome of one chromosome with SNP count and produces data frames containing SNPs in exons, introns
#' and intergenic regions of the chromosome.
#' @param chr_gene A data frame containing the genome of the species for one chromosome.
#' Each row of data frame represents a specific section of the genome.  The data 
#' should at least contain the listed columns:
#' \describe{
#'  \item{Region}{Labeled either "exon", "CDS", "5UTR", or "3UTR".}
#'  \item{Start}{The begining of the genomic region.}
#'  \item{Stop}{The end of the genomic region.}
#'  \item{GeneID}{The name of the gene to which the region belongs}
#'  \item{direction}{A "+" or "-" indicating the strand on which the gene segment is located}
#' }
#' @param chr_SNP A data frame containing the start and stop of each SNP in the specific 
#' chromosome of the genome. The start codon should be labled "Start," and the stop codon should be labeled "Stop."
#' 
#' @return A list with five components:
#' \describe{
#'  \item{SNP_introns}{A data frame containing all the SNPs in introns}
#'  \item{SNP_exons}{A data frame containing all the SNPs in exons. The data frame also contains information on the gene name and genomic region in 
#'  which the SNP is found.}
#'  \item{SNP_stream}{A data frame containing all the SNPs in the intergenic regions}
#'  \item{stream_withSNP}{Contains the start and stop codon of each intergenic region and the number of SNPs found int the region.}
#'  \item{chr_withSNP}{The orginal "chr_gene" data frame with an additional column containing the number of SNPs in each gene region.}
#'}
#' @export
do_one_chromosome <- function(chr_gene, chr_SNP, w = NULL, do_stream = TRUE){
  
  #### Getting Rid of repeat UTRs ####
  repeat_UTR <- chr_gene[grep("UTR", chr_gene[,3]),] # matrix of only UTRs
  take_out = c()
  for (i in 1:(length(repeat_UTR$Region)-1)){
    if(repeat_UTR$Region[i] == repeat_UTR$Region[i+1]){  # finds a repeat
      if(repeat_UTR$Region[i] == "5UTR"){  # if 5UTR it will take out the first of the repeats
        take_out = append(take_out,i) 
      }
      else{   # if 3UTR it will take out the last of the repeats
        take_out = append(take_out,(i+1))
      }
    }
  }
  index_out <- as.numeric(row.names(repeat_UTR)[take_out]) # gets index of UTRs to remove
  new_chr <- chr_gene[-(match(index_out, as.numeric(row.names(chr_gene)))),] # removes UTRs
  
  #### Adding up and down stream ####
  repeat_UTR2 <- repeat_UTR[-(take_out),] # all of UTRs with no repeats

  if(repeat_UTR2$Region[nrow(repeat_UTR2)] == "5UTR"){
    repeat_UTR2 <- repeat_UTR2[-nrow(repeat_UTR2),]
  }
  if(repeat_UTR2$Region[1] == "3UTR"){
    repeat_UTR2 <- repeat_UTR2[-1,]
  }
  pos <- which(repeat_UTR2$direction == "+")
  UTR_ordertmp <- rep(0,nrow(repeat_UTR2))
  UTR_ordertmp[pos] <- repeat_UTR2$Start[pos]
  UTR_ordertmp[-pos] <- repeat_UTR2$Stop[-pos]
  order_index <- order(UTR_ordertmp)
  UTR_order <- sort(UTR_ordertmp)
  
  index_start <- seq(2,nrow(repeat_UTR2)-1,2)
  index_stop <- seq(3,nrow(repeat_UTR2)-1,2)
  
  stream <- data.frame(Start = UTR_order[index_start],Stop=UTR_order[index_stop])
  
  #### Run IRanges on non-stream data sets #### 
  if(length(w) == 0){
    gr0 <- with(new_chr, IRanges(start=Start, end=Stop))
    gr1 <- with(chr_SNP, IRanges(start=Start, end=Stop))
    hits = findOverlaps(gr0, gr1) # queryHits are genes and subjectHits are SNPs
  }
  else{
    gr0 <- with(new_chr, IRanges(start=Start, end=Stop))
    gr1 <- with(chr_SNP, IRanges(start=Start, width=1))
    hits = findOverlaps(gr0, gr1) # queryHits are genes and subjectHits are SNPs
  }
  
  #### Add to chromosome table number of snps at each exon, UTR, TE ####
  SNP_tab <- as.numeric(table(queryHits(hits)))
  snp_count <- rep(0,nrow(new_chr))
  snp_count[as.numeric(row.names(table(queryHits(hits))))] = SNP_tab
  chr_withSNP <- cbind(new_chr,snp_count)

  # find SNPs that were in coding sequencuences
  SNP_gene_match <- rep(0,nrow(chr_SNP))
  SNP_Region <- rep(0,nrow(chr_SNP))
  SNP_gene_match[subjectHits(hits)] = new_chr$GeneID[queryHits(hits)]
  SNP_Region[subjectHits(hits)] = as.character(new_chr$Region[queryHits(hits)])
  
  SNP_with_data <- cbind(chr_SNP,SNP_Region,SNP_gene_match) # final table containing data
  SNP_exons <- SNP_with_data[subjectHits(hits),]
  
  if(do_stream==TRUE){
    SNP_no_CDS <- chr_SNP[which(SNP_with_data$SNP_gene_match=="0"),] # removes genes in CDS, UTR, TE

    #### Run IRANGES on down and up stream data and collect information ####
    # Run Iranges
    gr2 <- with(stream, IRanges(start=Start, end=Stop))
    gr3 <- with(SNP_no_CDS, IRanges(start=Start, end=Stop))
    stream_hits = findOverlaps(gr2, gr3)
  
    # get counts
    SNP_tab_stream <- as.numeric(table(queryHits(stream_hits)))
    snp_count_stream <- rep(0,nrow(stream))
    snp_count_stream[as.numeric(row.names(table(queryHits(stream_hits))))] = SNP_tab_stream
    stream_withSNP <- cbind(stream,snp_count_stream) # new final matrix of up and down stream info
  
    SNP_stream <- cbind(SNP_no_CDS[subjectHits(stream_hits),]) # final table with SNP info in up and down stream

    # table with introns
    SNP_introns <- SNP_no_CDS[-subjectHits(stream_hits),]
  
    list(SNP_introns = SNP_introns, SNP_exons = SNP_exons, SNP_stream = SNP_stream, stream_withSNP = stream_withSNP, chr_withSNP = chr_withSNP)
  }
  else{
    list(chr_withSNP = chr_withSNP, SNP_exons = SNP_exons )
  }
  }

