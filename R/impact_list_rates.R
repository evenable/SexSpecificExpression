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
impact_list_rates <- function(all_data,impact_list,gene_list, chr_num){
  gene_list <- as.matrix(gene_list)
  list_of_geneID <- list()
  for(i in 1:ncol(gene_list)){
    list_of_geneID[[i]] <- list()
    remove_indices <- which(gene_list[,i] == "")
    list_of_geneID[[i]]$name <- colnames(gene_list)[i]
    if(length(remove_indices) == 0){
      list_of_geneID[[i]]$genes <- gene_list[,i]
    }
    else{
      list_of_geneID[[i]]$genes <- gene_list[-(remove_indices),i]
    }
  }
  background_data <- background_impact_rate(all_data = all_data,impact_list = impact_list, chr_num = chr_num)
  
  calc_rate <- function(gene_names, impact_list, all_data, background_data){
    index_genes <- which(background_data$prepared_impact_list$GeneID %in% gene_names)
    rate <- sum(background_data$prepared_impact_list$impact_snp_count[index_genes])/sum(background_data$prepared_impact_list$gene_length[index_genes])
    rate
  }
  
  rates <- lapply(list_of_geneID, function(z) list(name = z$name, rate = calc_rate(gene_names = z$genes, impact_list = impact_list, all_data = all_data,background_data = background_data)))
  rate_table <- matrix(unlist(rates),ncol=2,byrow=TRUE)
  list(rate_table = rate_table, background_rate = background_data)
}