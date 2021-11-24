setwd("D:/Uppsala/Period_6/Applied_Bioinformatics/Monarchs/Github/Monarchs")
table = read.csv("migrant_vs_resident_Final_stats.csv")

options(scipen = 999) #changing R labels to reals instead of exponential

library(ggplot2)
library(dplyr)
library(plotly)
library(hrbrthemes)
library(lattice)

library(httr)
library(jsonlite)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(purrr)
library(Sushi)
library(egg)
library(cowplot)
library(lubridate)


test1 = table[order(-table$West_vs_resident_Fst), ]
westvresident = test1[1:236, ]

test2 = table[order(-table$East_vs_resident_Fst), ]
eastvresident = test2[1:236, ]

test3 = table[order(-table$migrant_vs_resident_Fst), ]
migrantvresident = test3[1:236, ]

top2 = test3[1:472, ]
top2_migrantvresident = write.csv(top2, "top2migrantvresident.csv")

top5 = test3[1:1180,]
top5_migrantvresident = write.csv(top5, "top5migrantvresident.csv")

florida_west = table[order(-table$West_vs_Florida_FST), ]
florida_east = table[order(-table$East_vs_Florida_FST), ]


#write.csv(westvresident, "West_vs_Resident.csv")
#write.csv(eastvresident, "East_vs_Resident.csv")
#write.csv(migrantvresident, "migrant_vs_Resident.csv")
#write.csv(florida_west, "westvflorida.csv")
#write.csv(florida_east, "eastvflorida.csv")



table2 = subset(table, table$Chr_number==16)
sorted_table2 = table2[order(-table2$migrant_vs_resident_Fst), ]
top5_chr16 = sorted_table2[1:43,]

plot(table2$BIN_END, table2$migrant_vs_resident_Fst)
lines(sorted_table2$BIN_START, sorted_table2$migrant_vs_resident_Fst)

#p3 <- plot_ly(x = top5_chr16$BIN_END, y = top5_chr16$migrant_vs_resident_Fst, type="scatter", mode="markers")
#p3


#plot(top5_chr16$BIN_END, top5_chr16$migrant_vs_resident_Fst)

#plot(top5_chr16$BIN_START, top5_chr16$migrant_vs_resident_Fst)


#trend.line(top5_chr16$BIN_END, top5_chr16$migrant_vs_resident_Fst, plot=TRUE, main="Linear")

newtop5 = top5_chr16[order(top5_chr16$BIN_START), ]
plot(newtop5$BIN_START, newtop5$migrant_vs_resident_Fst)
lines(newtop5$BIN_START, newtop5$migrant_vs_resident_Fst)

ggplot(newtop5, mapping=aes(BIN_START,migrant_vs_resident_Fst)) + geom_line() + labs(y="Migrant vs Resident Fst",x="Position") + xlim(0,1000000)
plot1 = ggplot(newtop5, mapping=aes(BIN_START,migrant_vs_resident_Fst)) + geom_line() + labs(y="Migrant vs Resident Fst",x="Position") + xlim(0,1000000)

#Gene plot
data = read.table("danaus_plexippus_v3_core_32_85_1.gff")
t=c('DPSCF300098','DPSCF300147','DPSCF300267','DPSCF300067','DPSCF300240','DPSCF300322','DPSCF300095','DPSCF300353','DPSCF300170','DPSCF300127','DPSCF300331','DPSCF300114','DPSCF300105','DPSCF300047','DPSCF300236')

data2 = subset(data,data$V1==t)
genes = subset(data2,data2$V3=="gene")
genes = genes[order(genes$V4),]

table_sort = subset(table,table$CHROM==t)
#plot(genes$V4, factor(genes$V9))
genes$V9 = sub(";.*", "", genes$V9)  

dotplot(genes$V9~genes$V4)

ggplot(genes, mapping=aes(V4,V9)) +geom_line() + geom_point() + labs(y="Gene ID",x="Gene Start Position") + xlim(0,1000000)
plot2 = ggplot(genes, mapping=aes(V4,V9)) +geom_line() + geom_point() + labs(y="Gene ID",x="Gene Start Position") + xlim(0,1000000)
# + geom_line(aes(x=V4,y=V9)) + geom_line(aes(x=V5,y=V9))


plotGenes(genes, chrom=16)

cowplot::plot_grid(plot1, plot2, align = "v", ncol = 1, rel_heights = c(0.25, 0.75))
egg::ggarrange(plot1, plot2, heights = c(0.25, 0.75))

unique(sorted_table2$CHROM)
