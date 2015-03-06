library(SexSpecificExpression)
library(ggplot2)
library(xtable)
source("~/SexSpecificExpression/not_package/table_functions.R")
load("~/SexSpecificExpression/data/all_compiled_data.rda")
load("~/SexSpecificExpression/data/SNP_list.rda")

total_bargraph <- total_bargraph(all_data = all_compiled_data)

pdf("tables_and_charts/total_bargraph2",width=4,height=5)
total_bargraph
dev.off()

total_table(all_data = all_data)
