#' Prepare Genome for the Functions
#' 
#' Description
#' @param genome
#' 
#' @return The edited genome with with the GeneID column containing only
#' the names of the gene with no extra characters.
#'
#' @export
prepare_genome <- function(genome){
gene_names_temp <- gsub("^.*A","A",as.character(genome$GeneID))
genome$GeneID <- gsub("\\..*","",gene_names_temp)
index <- grep("TE",genome$GeneID)
genome$GeneID[index] <- gsub("[\\].*","",genome$GeneID[index])
UTR <- grep(pattern="UTR",genome$Region)
new_genome <- genome[-UTR,]
ID_and_count_temp <- data.frame(GeneID = new_genome$GeneID, snp_count = (new_genome$snp_count), gene_length = abs(new_genome$Start - new_genome$Stop))
ID_and_count <- aggregate(. ~ GeneID, data = ID_and_count_temp, FUN = sum)
ID_and_count
}