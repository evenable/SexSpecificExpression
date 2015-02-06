#' Title
#' 
#' Description
#' @param param some stuff
#' @param param some stuff
#' 
#' @return list
#' \describe{
#'  \item{stuff}{stuff }
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
  list(introns = introns, exons = exons, stream = stream, genome_with_SNP = genome_with_SNP, up_and_down_stream = up_and_down_stream, SNP_numbers = list(num_in_exons = num_in_exons, num_in_introns = num_in_introns, num_in_stream = num_in_stream))
}