#background rate of all genes
#background rate for specific gene lists
background_impact_rate <- function(all_data,gene_list){
  prepared_genome <- prepare_genome(genome = all_data$genome_with_SNP)
  
  gene_names <- all_data$exons$SNP_gene_match[which(all_data$exons$Start %in% gene_list$POS)]
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
  ave_rate
}