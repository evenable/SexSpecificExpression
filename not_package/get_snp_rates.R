library(SexSpecificExpression)
load("~/SexSpecificExpression/data/compiled_data.rda")
load("~/SexSpecificExpression/data/genelist_FirstSet.rda")
load("~/SexSpecificExpression/data/gene_IDs.rda")
all_info <- all_chromosomes(compiled_SNP_data = compiled_data)
comparisons <- compare_background_rate(gene_list =  genelists_FirstSet,all_IDs = gene_IDs,genome =all_info$genome_with_SNP,rep_num = 100)

data_list <- split(x=comparison_data, f=comparison_data$name)

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

