index_genome <- which(all_data$genome_with_SNP$snp_count == 0)
index2 <- grep("UTR",all_data$genome_with_SNP$Region)
all_genes <- all_data$genome_with_SNP[-index2,]
index_genome <- which(all_genes$snp_count == 0)
all_genes <- all_genes$GeneID[-index_genome]
all_gene_names_temp <- gsub("^.*A","A",as.character(all_genes))
all_gene_names <- gsub("\\..*","",all_gene_names_temp)
TE <- grep("TE",all_gene_names)
all_gene_names[TE] <- gsub("[\\].*","",all_gene_names[TE])
genes_with_snps <- unique(all_gene_names)


all_genes2 <- all_data$genome_with_SNP
index_genome2 <- which(all_genes2$snp_count == 0)
all_genes22 <- all_genes2$GeneID[-index_genome2]
all_gene_names_temp2 <- gsub("^.*A","A",as.character(all_genes22))
all_gene_names2 <- gsub("\\..*","",all_gene_names_temp2)
TE2 <- grep("TE",all_gene_names2)
all_gene_names2[TE2] <- gsub("[\\].*","",all_gene_names2[TE2])
Ugenes_with_snps <- unique(all_gene_names2)



