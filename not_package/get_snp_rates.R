library(SexSpecificExpression)
load("~/SexSpecificExpression/R/sysdata/snp_rate_table.rda")
load("~/SexSpecificExpression/data/gene_lists_three_sets.rda")
load("~/SexSpecificExpression/data/all_compiled_data.rda")

gene_list <- gene_lists_three_sets
data_list <- split(x=snp_rate$all_trials, f=snp_rate$all_trials$name)
rate_table <- snp_rate_table
prepared_genome <- prepare_genome(genome = all_compiled_data$genome_with_SNP)

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

length_tmp <- lapply(list_of_geneID, function(z) {sum(prepared_genome$gene_length[(prepared_genome$GeneID %in% z$genes)])})
gene_length <- data.frame(name = colnames(gene_lists_three_sets), length = unlist(length_tmp))
#### Poisson #### 
p_values <- matrix(0,nrow=nrow(rate_table),ncol=(ncol(rate_table)-2))
rownames(p_values) <- rep("",nrow(rate_table))
for(i in 1:nrow(rate_table)){
  rownames(p_values)[i] <- as.character(rate_table$name[i])
  true_rate <- rate_table$true_rates[i]
  true_count <- round(true_rate*gene_length$length[i])
  length <- gene_length$length[i]
  for(j in 1:(ncol(rate_table)-2)){
    expected_rate <- rate_table[i,j+1]
    if( true_rate > rate_table[i,j+1]){
      p_values[i,j] = poisson.test(x = true_count, T = length,r = expected_rate, alternative = "greater")$p.value
    }
    else{
      p_values[i,j] = poisson.test(x = true_count, T = length,r = expected_rate, alternative = "less")$p.value
    }
  }
}
colnames(p_values) <- colnames(rate_table)[2:(ncol(rate_table)-1)]

#### Test Normality ####
p_values <- matrix(0,nrow=length(data_list),ncol=(ncol(data_list[[1]])-1))
rownames(p_values) <- rep("",length(data_list))
for(i in 1:length(data_list)){
  rownames(p_values)[i] <- as.character(data_list[[i]]$name[1])
  for(j in 1:(ncol(data_list[[i]])-1)){
    p_values[i,j] = shapiro.test(data_list[[i]][,j])$p.value
  }
}
colnames(p_values) <- colnames(data_list[[1]])[1:3]

#### T-test or Wilcox Test ####
t_test <- matrix(0,nrow=length(data_list),ncol=(ncol(data_list[[1]])-1))
rownames(t_test) <- rep("",length(data_list))
for(i in 1:length(data_list)){
  rownames(t_test)[i] <- as.character(data_list[[i]]$name[1])
  for(j in 1:(ncol(data_list[[i]])-1)){
    if(p_values[i,j]>0.05){
      t_test[i,j] = t.test(x = data_list[[i]][,j],alternative = "two.sided",mu = comparison_table$true_rates[i])$p.value
    }
    else{
      t_test[i,j] = wilcox.test(x = data_list[[i]][,j],alternative = "two.sided",mu = comparison_table$true_rates[i])$p.value
    }
  }
}
colnames(t_test) <- colnames(p_values)

