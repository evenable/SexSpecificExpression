#' Compare SNP Rates in Gene Lists to Background Rates
#' 
#' This function both calculates that SNP rate per nucleotide in specific gene lists and 
#' simulates background rates for comparision to each gene list.  The background rates are simulated
#' using random lists of the same size as the gene lists and are simulated using three sets of
#' genes: all genes, only TEs, and all genes but TEs.
#' @param gene_list A data frame containing the a list of genes to be tested in each column.
#' @param all_IDs A list of all genes in the genome.
#' @param genome The genome of the species with the snp count of each gene region.
#' @param rep_num The number of simulations replications used to estimate the background rateß.
#' 
#' @return A list with two elements:
#' \describe{
#'  \item{all_trails}{A data frame containing all of the simulated background rates for each gene list.}
#'  \item{rates}{A table containing the true rates and simulated rates for each gene list.  There are 
#'  three simulated rates: rates including all genes, rates including only TEs, and rates containing
#'  only genes without TEs.}
#'}
#' @export
compare_background_rate <- function(gene_list, all_IDs, genome, rep_num){
  #### Prepare Genome ####
  prepared_genome <- prepare_genome(genome)
  
  #### Function for Generating and Running Random Sample ####
  simulate_background_rate <- function(N, all_IDs, prepared_genome){
  TEs <- grep(pattern = "TE", prepared_genome$GeneID)
  index_TEs <- grep(pattern = "TE", all_IDs)
  different_genomes <- list()
  different_genomes$all$data <- prepared_genome
  different_genomes$TE$data <- prepared_genome[TEs,]
  different_genomes$noTE$data <- prepared_genome[-TEs,]
  different_genomes$all$gene_index <- all_IDs[sample(1:length(all_IDs), N, replace=FALSE)]
  different_genomes$TE$gene_index <- all_IDs[TEs][sample(1:length(all_IDs[TEs]), N, replace=FALSE )]
  different_genomes$noTE$gene_index <- all_IDs[-TEs][sample(1:length(all_IDs[-TEs]), N, replace=FALSE )]
  snp_rate <- lapply(different_genomes, function(x) calculate_snp_rate(GeneIDs = x$gene_index,prepared_genome = x$data))
  snp_rate
  }
  
  #### Make the Data Table Lists and Remove NAs ####
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
  
  #### Run the simulations for rep_num times ####
  simulated_rates <- lapply(list_of_geneID, function(z) {
    lapply(1:rep_num, function(y) {
    list( rates = simulate_background_rate(N = length(z$genes), all_IDs = all_IDs,prepared_genome = prepared_genome), name = z$name)})})
  
  slurp_rates <- function(y){
    data.frame(all_rate = y$rates$all, TE_rate = y$rates$TE, noTE_rate = y$rates$noTE, name = y$name)
  }
  rates_tmp1 <- lapply(1:length(simulated_rates), function(rep) {lapply(simulated_rates[[rep]], function(x) {df <- slurp_rates(x); df})})
  rates_tmp2 <- unlist(rates_tmp1, recursive = FALSE)
  rates_df <- do.call(what = rbind, args = rates_tmp2)
  rates <- aggregate(. ~ name, data = rates_df, FUN = sum)
  rates[,-1] <- rates[,-1]/rep_num
  
  #### Calculate the True Rates ####
  real_rates <- lapply(list_of_geneID, function(z) calculate_snp_rate(GeneIDs = z$genes, prepared_genome = prepared_genome))
  
  rates$true_rates <- unlist(real_rates)
  
  list(all_trials = rates_df, rates = rates)
}