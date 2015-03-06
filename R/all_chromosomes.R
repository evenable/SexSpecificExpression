#' Compile all chromosome lists into data frames.
#' 
#' This function takes the annotated SNP data for each chromosome and compiles all the data into five data frames containing
#' the entire genome, a table with the intergenic data, and tables containing SNPs in exons, introns, and intergenic regions.
#' @param compiled_SNP_data some A list containing the annotated data of each chromosome.
#' 
#' @return A list of five components.
#' \describe{
#'  \item{introns}{A data frame containing all the stop and start codons of all the SNPs located in introns.}
#'  \item{exons}{A data frame containing all the stop and start codons, gene location, and genomic region location of all the SNPs in exons}
#'  \item{stream}{A data frame containing all the start and stop codons for the SNPs in the intergenic regions.}
#'  \item{genome_with_SNP}{The orginal "allgenes" genome data frame with an additional column containing the SNP count for each region.}
#'  \item{stream_with_SNP}{A data frame containing the start and stop codon for each intergenic region and the SNP count for each region.}
#'}
#' @export
all_chromosomes <- function(compiled_SNP_data){
  introns <- do.call(what = rbind, args = lapply(compiled_SNP_data, function(x) x$SNP_introns))
  exons <- do.call(what = rbind, args = lapply(compiled_SNP_data, function(x) x$SNP_exons))
  stream <- do.call(what = rbind, args = lapply(compiled_SNP_data, function(x) x$SNP_stream))
  genome_with_SNP <- do.call(what = rbind, args = lapply(compiled_SNP_data, function(x) x$chr_withSNP))
  up_and_down_stream <- do.call(what = rbind, args = lapply(compiled_SNP_data, function(x) x$stream_withSNP))
  num_in_exons <- nrow(exons)
  num_in_introns <- nrow(introns)
  num_in_stream <- nrow(stream)
  list(introns = introns, exons = exons, stream = stream, genome_with_SNP = genome_with_SNP, stream_with_SNP = up_and_down_stream)
}