#' Title
#' 
#' Description
#' @param GeneIDs some stuff
#' @param genome
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
calculate_snp_rate <- function(GeneIDs, prepared_genome){
  index_genes <- which(prepared_genome$GeneID %in% GeneIDs)
  snp_rate <- sum(prepared_genome$snp_count[index_genes])/sum(prepared_genome$gene_length[index_genes])
  snp_rate
}