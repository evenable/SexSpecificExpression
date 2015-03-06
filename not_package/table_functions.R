total_bargraph <- function(all_data){
  index_TE <- grep(pattern="TE",all_data$genome_with_SNP$GeneID)
  exon_data = all_data$genome_with_SNP[-index_TE,]
  FUTR_index <- grep(pattern="5UTR",exon_data$Region)
  TUTR_index <- grep(pattern="3UTR",exon_data$Region)
  exon = sum(exon_data$snp_count[-c(FUTR_index,TUTR_index)])
  FUTR <- sum(exon_data$snp_count[FUTR_index])
  TUTR <- sum(exon_data$snp_count[TUTR_index])
  TE = sum(all_data$genome_with_SNP$snp_count[index_TE])
  intron = nrow(all_data$introns)
  total = nrow(SNP_list)
  intergene = total - (TE + intron + exon + FUTR + TUTR)
  snp_percent = 100*(c(intergene,FUTR, TUTR, exon, intron, TE)/total)
  df <- data.frame(genome_loc = c("Intergenic","5UTR","3UTR","Exon","Intron","TE"),
                   snp_percent = snp_percent)
  a <- ggplot(data=df, aes(x=genome_loc,y=snp_percent, fill=genome_loc)) + 
    geom_bar(color="black",stat="identity") + guides(fill=FALSE) +
    ylab("Percent of SNPs") + xlab("Genomic Region") +
    scale_x_discrete(limits = c("Intergenic","5UTR","3UTR","Exon","Intron","TE")) + 
    geom_text(aes(label=round(snp_percent,digits=1), y=snp_percent+1),size=4) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
         axis.line = element_line(color = "black"), panel.background = element_blank())
  return(a)
}

total_table <- function(all_data){
  index_TE <- grep(pattern="TE",all_data$genome_with_SNP$GeneID)
  exon_data = all_data$genome_with_SNP[-index_TE,]
  FUTR_index <- grep(pattern="5UTR",exon_data$Region)
  TUTR_index <- grep(pattern="3UTR",exon_data$Region)
  exon = sum(exon_data$snp_count[-c(FUTR_index,TUTR_index)])
  FUTR <- sum(exon_data$snp_count[FUTR_index])
  TUTR <- sum(exon_data$snp_count[TUTR_index])
  TE = sum(all_data$genome_with_SNP$snp_count[index_TE])
  intron = nrow(all_data$introns)
  intergene <- nrow(SNP_list) - (intron + TE + exon + FUTR + TUTR)
  SNPs <- c(intergene,FUTR,TUTR, exon, intron, TE)
  GenomicRegion <- c("Intergenic","5UTR","3UTR", "Exons", "Introns","TE")
  data <- cbind(GenomicRegion,SNPs)
  data.table <- xtable(data)
  colnames(data.table)[1] <- "Genomic Region" 
  print(data.table, include.rownames=FALSE)
}