#' Annotate Genome with SNPs
#' 
#' This function annotates the genome with SNP count and produces data frames containing SNPs in exons, introns
#' and intergenic regions of the genome.
#' @param allgenes A data frame containing the entire genome of the species.
#' Each row of data frame represents a specific section of the genome.  The data 
#' should at least contain the listed columns:
#' \describe{
#'  \item{Header}{Contains a number representing the chromosome in which the gene region is located. All chromosomes with numbers should come first in the list.}
#'  \item{Region}{Labeled either "exon", "CDS", "5UTR", or "3UTR".}
#'  \item{Start}{The begining of the genomic region.}
#'  \item{Stop}{The end of the genomic region.}
#'  \item{GeneID}{The name of the gene to which the region belongs}
#'  \item{direction}{A "+" or "-" indicating the strand on which the gene segment is located}
#' }
#' @param SNPs A data frame containing the chromosome location, start, and stop of each SNP. The start 
#' codon should be labled "Start," the stop codon should be labeled "Stop," and the chromosome should be
#' labeled "Chromosome."
#' @param chr_num The number of chromosomes in the genome.
#' 
#' @return A list with five components:
#' \describe{
#'  \item{introns}{A data frame containing all the stop and start codons of all the SNPs located in introns.}
#'  \item{exons}{A data frame containing all the stop and start codons, gene location, and genomic region location of all the SNPs in exons}
#'  \item{stream}{A data frame containing all the start and stop codons for the SNPs in the intergenic regions.}
#'  \item{genome_with_SNP}{The orginal "allgenes" genome data frame with an additional column containing the SNP count for each region.}
#'  \item{stream_with_SNP}{A data frame containing the start and stop codon for each intergenic region and the SNP count for each region.}
#'}
#' @export
categorize_SNPs <- function(allgenes, SNPs, chr_num){
  
  gene_list <- split(x = allgenes, f = allgenes$Header)
  snp_list <- split(x=SNPs, f = SNPs$Chromosome)
  
  types <- list()
  for(i in 1:chr_num){
    types[[i]] <- i
  }
  SNP_categorized <- lapply(types, function(x) 
    do_one_chromosome(chr_gene = gene_list[[x]], chr_SNP = snp_list[[x]]))
  all_compiled_data <- all_chromosomes(compiled_SNP_data = SNP_categorized)
  all_compiled_data
}