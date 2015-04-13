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
background_impact_rate <- function(all_data, impact_list,chr_num){
  prepared_genome <- prepare_genome(genome = all_data$genome_with_SNP)
  TU_index <- which(all_data$exons$SNP_Region == "3UTR")
  FU_index <- which(all_data$exons$SNP_Region == "5UTR")
  exon_table <- all_data$exons[-c(TU_index,FU_index),]
  impact <- split(x = impact_list,f=impact_list$X.CHROM)
  exons <- split(x = exon_table, f = exon_table$Chromosome)

  gene_names_list <- lapply(1:chr_num, function(x) {gene_names <-exons[[x]]$SNP_gene_match[which(exons[[x]]$Stop %in% impact[[x]]$POS)]})
  gene_names <- unlist(gene_names_list)
  freq_table <- as.data.frame(ftable(gene_names))
  freq_table$gene_names <- as.character(factor(all_data$genome_with_SNP$GeneID)[freq_table$gene_names])
  gene_names_temp <- gsub("^.*A","A",as.character(freq_table$gene_names))
  freq_table$gene_names <- gsub("\\..*","",gene_names_temp)
  TE_index <- grep("TE",freq_table$gene_names)
  if(length(TE_index) != 0){
    freq_table$gene_names[TE_index] <- gsub("[\\].*","",freq_table$gene_names[TE_index])
  }
  freq_table_norep <- aggregate(. ~ gene_names, data = freq_table, FUN = sum)
  
  impact_snp_count <- rep(0,nrow(prepared_genome))
  gene_index <- match(freq_table_norep$gene_names, prepared_genome$GeneID)
  impact_snp_count[gene_index] <- freq_table_norep$Freq
  genome_data <- cbind(prepared_genome,impact_snp_count)
  rates <- (genome_data$impact_snp_count)/(genome_data$gene_length)
  ave_rate <- sum(rates)/length(rates)
  list(average_rate = ave_rate, prepared_impact_list = genome_data)
}