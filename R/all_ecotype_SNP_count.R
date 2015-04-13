#' Compile all chromosome lists into data frames.
#' 
#' Description
#' @param compiled_SNP_data some A list containing the annotated data of each chromosome.
#' 
#' @return A data frame.
#' @export
all_ecotype_SNP_count <- function(genome_file, gene_ids, gene_directory, ecotype_list,chr_num){
  gene_list <- split(x = genome_file, f = genome_file$Header)
  file_list <- list.files(gene_directory)
  ecotype_list <- list()
  index <- 0
  for (file in file_list){
    index <- index + 1
    f_name <- load(paste(gene_directory,file,sep=""))
    f_data <- get(f_name)
    ecotype_list[[index]] <- list()
    ecotype_list[[index]]$name <- as.character(f_data$type[1])
    ecotype_list[[index]]$snp_list <- split(x=f_data, f = f_data$Chromosome)
  }
  
  
  SNP_categorized <- lapply(ecotype_list, function(z) {list(name = z$name, data = lapply(1:chr_num, function(x) 
    do_one_chromosome(chr_gene = gene_list[[x]], chr_SNP = z$snp_list[[x]],w = 1,do_stream = FALSE)))})
  
  all_compiled_data <- lapply(SNP_categorized, function(y) {list(name = y$name, data = genome_with_SNP <- do.call(what = rbind, args = lapply(y$data, function(x) x$chr_withSNP)))})
  
  snp_count_tables <- lapply(all_compiled_data, function(x) {list(name = x$name, genome = prepare_genome(genome = x$data))})
  
  snp_count <- lapply(snp_count_tables, function(y) {list(name=y$name, count = y$genome$snp_count[match(gene_ids, y$genome$GeneID)])})
  
  table <- data.frame(Gene = gene_ids)
  tmp_table <- do.call(what = cbind, args = lapply(snp_count, function(x) x$count))
  colnames(tmp_table) = lapply(snp_count, function(x) unlist(x$name))
  ecotype_table <- cbind(table,tmp_table)
}
