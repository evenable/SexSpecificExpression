library(GenomicRanges)
# Separating into lists
allgenes_list <- split(x = allgenes_and_TE3, f = allgenes_and_TE3$Header)
chr1 <- allgenes_list[[1]]

#### Getting Rid of repeat UTRs ####
test <- chr1[grep("UTR", chr1[,3]),]
take_out = c()
for (i in 1:(length(test$Region)-1)){
  if(test$Region[i] == test$Region[i+1]){
    if(test$Region[i] == "5UTR"){
      take_out = append(take_out,i)  
    }
    else{
      take_out = append(take_out,(i+1))
    }
  }
}
index_out <- as.numeric(row.names(test)[take_out])
new_chr1 <- chr1[-(match(index_out, as.numeric(row.names(chr1)))),]

#### Adding up and down stream ####
direction <- c()
starts <- 0
test2 <- test[-(take_out),]
up_down1 <- test2
up_down1$Region = "down"
for(i in which(test2$Region == "5UTR")){
    if(test2$Start[i] < test2$Start[i+1]){
      starts <- test2$Start[i] - 1
      up_down1$Stop[i] <- starts
      up_down1$Start[i] <- starts - 1000
      direction[i] <- "+"
    }
    else{
      starts <- test2$Stop[i] + 1
      up_down1$Start[i] <- starts
      up_down1$Stop[i] <- starts + 1000
      direction[i] <- "-"
    }
  }
for(i in which(test2$Region == "3UTR")){
  up_down1$Region[i] = "up"
  if(test2$Start[i] < test2$Start[i-1]){
    starts <- test2$Start[i] - 1
    up_down1$Stop[i] <- starts
    up_down1$Start[i] <- starts - 1000
    direction[i] <- "-"
  }
  else{
    starts <- test2$Stop[i] + 1
    up_down1$Start[i] <- starts
    up_down1$Stop[i] <- starts + 1000
    direction[i] <- "+"
  }
}  
up_down <- cbind(up_down1,direction)

#### Run on non down stream ####
chr_SNP <- split(x = SNP_list, f = SNP_list$Chromosome)
SNP1 <- chr_SNP[[1]]
gr0 <- with(new_chr1, IRanges(start=Start, end=Stop))
gr1 <- with(SNP1, IRanges(start=Start, end=Stop))
hits = findOverlaps(gr0, gr1) # queryHits are genes and subjectHits are SNPs

 get count for each gene
SNP_tab <- as.numeric(table(queryHits(hits)))
snp_count <- rep(0,nrow(new_chr1))
snp_count[as.numeric(row.names(table(queryHits(hits))))] = SNP_tab
chr1_withSNP <- cbind(new_chr1,snp_count)
total_SNPs_in_coding <- sum(snp_count)


# find way to mark the distance from the start codon
fiveU <- up_down[(which(up_down$Region == "down")),]
fiveU_start <- rep(0, nrow(fiveU))
fiveU_start[which(fiveU$direction == "+")] = fiveU$Stop[which(fiveU$direction == "+")]
fiveU_start[which(fiveU$direction =="-")] = fiveU$Start[which(fiveU$direction == "-")]
fiveU_tot <- cbind(fiveU_start, fiveU$GeneID)

threeU <- up_down[(which(up_down$Region == "up")),]
threeU_start <- rep(0, nrow(threeU))
threeU_start[which(threeU$direction == "+")] = threeU$Start[which(threeU$direction == "+")]
threeU_start[which(threeU$direction =="-")] = threeU$Stop[which(threeU$direction == "-")]
threeU_tot <- cbind(threeU_start, threeU$GeneID)

SNP_gene_match <- rep(0,nrow(SNP1))
SNP_distance_from_start <- rep(0,nrow(SNP1))
SNP_distance_from_end <- rep(0,nrow(SNP1))
SNP_gene_match[subjectHits(hits)] = new_chr1$GeneID[queryHits(hits)]

for (i in 1:nrow(SNP1)){
  name <- SNP_gene_match[i]
  if(name != 0 & name %in% fiveU_tot[,2] & name %in% threeU_tot[,2]){
    SNP_distance_from_start[i] <- SNP1$Start[i] - fiveU_tot[which(fiveU_tot[,2] == name),1]
    SNP_distance_from_end[i] <- threeU_tot[which(threeU_tot[,2]== name),1] - SNP1$Stop[i] 
  }
}

SNP_with_data <- cbind(SNP1,SNP_gene_match,SNP_distance_from_start,SNP_distance_from_end)

# Take out hits
SNP_no_CDS <- SNP1[which(SNP_with_data$SNP_gene_match=="0"),]

#### Run on down stream ####
gr2 <- with(up_down, IRanges(start=Start, end=Stop))
gr3 <- with(SNP_no_CDS, IRanges(start=Start, end=Stop))
hits_down = findOverlaps(gr2, gr3)

# get counts and distance from start stop codon
SNP_tab_stream <- as.numeric(table(queryHits(hits_down)))
snp_count_stream <- rep(0,nrow(up_down))
snp_count_stream[as.numeric(row.names(table(queryHits(hits_down))))] = SNP_tab_stream
stream_withSNP <- cbind(up_down,snp_count_stream)

#distances
plus_down <- up_down$Region[queryHits(hits_down)] == "down" & up_down$direction[queryHits(hits_down)] == "+"
minus_down <- up_down$Region[queryHits(hits_down)] == "down" & up_down$direction[queryHits(hits_down)] == "-"
plus_up <- up_down$Region[queryHits(hits_down)] == "up" & up_down$direction[queryHits(hits_down)] == "+"
minus_up <- up_down$Region[queryHits(hits_down)] == "up" & up_down$direction[queryHits(hits_down)] == "-"
SNP_distance <- plus_down*(up_down$Stop[queryHits(hits_down)] - SNP_no_CDS$Start[subjectHits(hits_down)])+ minus_down*(up_down$Start[queryHits(hits_down)] - SNP_no_CDS$Start[subjectHits(hits_down)]) +plus_up*(up_down$Stop[queryHits(hits_down)] - SNP_no_CDS$Start[subjectHits(hits_down)]) + minus_up*(up_down$Start[queryHits(hits_down)] - SNP_no_CDS$Start[subjectHits(hits_down)])

stream_with_data <- cbind(SNP_no_CDS[subjectHits(hits_down),],queryHits(hits_down),SNP_distance)

# Take out hits
SNP_introns <- SNP_no_CDS[-subjectHits(hits_down),]

#### SNP content in expressed versus not expressed ####
subset_genes <- gene_IDs[1:100]
SNPs_and_genes <- all_chromosomes(compiled_SNP_data = compiled_SNP_data2)
genes_with_SNP_count <- SNPs_and_genes$genome_with_SNP
gene_names_temp <- gsub("^.*A","A",as.character(genes_with_SNP_count$GeneID))
genes_with_SNP_count$GeneID <- gsub("\\..*","",gene_names_temp)
ID_and_count_temp <- data.frame(GeneID = genes_with_SNP_count$GeneID, snp_rate = (genes_with_SNP_count$snp_count)/abs(genes_with_SNP_count$Start - genes_with_SNP_count$Stop))
ID_and_count <- aggregate(. ~ GeneID, data = ID_and_count_temp, FUN = sum)
index_genes <- which(ID_and_count$GeneID %in% subset_genes)

#### Background Rate function - making lists and removing NAs ####
list_of_geneID <- list()
for(i in 1:ncol(genelists_FirstSet)){
  list_of_geneID[[i]] <- list()
  remove_indices <- which(genelists_FirstSet[,i] == "")
  list_of_geneID[[i]]$name <- colnames(genelists_FirstSet)[i]
  list_of_geneID[[i]]$genes <- genelists_FirstSet[-(remove_indices),i]
}
