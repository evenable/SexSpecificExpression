library(SexSpecificExpression)
library(ggplot2)
library(xtable)
load("not_package/table_functions.R")
load("~/SexSpecificExpression/data/compiled_data.rda")
load("~/SexSpecificExpression/data/SNP_list.rda")

all_data <- all_chromosomes(compiled_SNP_data = compiled_data)
total_bargraph <- total_bargraph(all_data = all_data)

pdf("tables_and_charts/total_bargraph",width=4,height=5)
total_bargraph
dev.off()

total_table(all_data = all_data)