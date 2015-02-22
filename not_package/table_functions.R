total_bargraph <- function(all_data){
  index_TE <- grep(pattern="TE",all_data$genome_with_SNP$GeneID)
  exon = sum(all_data$genome_with_SNP$snp_count[-index_TE])
  TE = sum(all_data$genome_with_SNP$snp_count[index_TE])
  intron = nrow(all_data$introns)
  total = nrow(SNP_list)
  intergene = total - (TE + intron + exon)
  snp_percent = 100*(c(intergene,exon, intron, TE)/total)
  df <- data.frame(genome_loc = c("Intergenic","Exon","Intron","TE"),
                   snp_percent = snp_percent)
  a <- ggplot(data=df, aes(x=genome_loc,y=snp_percent, fill=genome_loc)) + 
    geom_bar(color="black",stat="identity") + guides(fill=FALSE) +
    ylab("Percent of SNPs") + xlab("Genomic Region") +
    scale_x_discrete(limits = c("Intergenic","Exon","Intron","TE")) + 
    geom_text(aes(label=round(snp_percent,digits=1), y=snp_percent+1),size=4)
  return(a)
}

total_table <- function(all_data){
  index_TE <- grep(pattern="TE",all_data$genome_with_SNP$GeneID)
  exon = sum(all_data$genome_with_SNP$snp_count[-index_TE])
  TE = sum(all_data$genome_with_SNP$snp_count[index_TE])
  intron = nrow(all_data$introns)
  intergene <- nrow(SNP_list) - (intron + TE + exon)
  SNPs <- c(intergene, exon, intron, TE)
  GenomicRegion <- c("Intergenic", "Exons", "Introns","TE")
  data <- cbind(GenomicRegion,SNPs)
  data.table <- xtable(data)
  colnames(data.table)[1] <- "Genomic Region" 
  print(data.table, include.rownames=FALSE)
}