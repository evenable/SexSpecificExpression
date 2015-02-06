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
calculate_snp_rate <- function(GeneIDs, genome){
  gene_names_temp <- gsub("^.*A","A",as.character(genome$GeneID))
  genome$GeneID <- gsub("\\..*","",gene_names_temp)
  ID_and_count_temp <- data.frame(GeneID = genome$GeneID, snp_count = (genome$snp_count), gene_length = abs(genome$Start - genome$Stop))
  ID_and_count <- aggregate(. ~ GeneID, data = ID_and_count_temp, FUN = sum)
  index_genes <- which(ID_and_count$GeneID %in% GeneIDs)
  snp_rate <- sum(ID_and_count$snp_count[index_genes])/sum(ID_and_count$gene_length)
  snp_rate
}