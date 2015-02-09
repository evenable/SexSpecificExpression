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
prepare_genome <- function(genome){
gene_names_temp <- gsub("^.*A","A",as.character(genome$GeneID))
genome$GeneID <- gsub("\\..*","",gene_names_temp)
index <- grep("TE",genome$GeneID)
genome$GeneID[index] <- gsub("[\\].*","",genome$GeneID[index])
ID_and_count_temp <- data.frame(GeneID = genome$GeneID, snp_count = (genome$snp_count), gene_length = abs(genome$Start - genome$Stop))
ID_and_count <- aggregate(. ~ GeneID, data = ID_and_count_temp, FUN = sum)
ID_and_count
}