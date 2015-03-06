#' Title
#' 
#' Description
#' @param gene_list some stuff
#' @param snp_list some stuff
#' @param chr_num
#' @param bpairs
#' 
#' @return list
#' \describe{
#'  \item{SNP_introns}{stuff }
#'  \item{SNP_exons}{}
#'  \item{SNP_stream}{}
#'  \item{stream_withSNP}{}
#'  \item{chr_withSNP}{}
#'}
#' @export
categorize_SNPs <- function(gene_list, snp_list, chr_num){
  types <- list()
  for(i in 1:chr_num){
    types[[i]] <- i
  }
  SNP_categorized <- lapply(types, function(x) 
    do_one_chromosome(chr_gene = gene_list[[x]], chr_SNP = snp_list[[x]]))
  SNP_categorized
}