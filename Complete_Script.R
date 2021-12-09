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

A=intersect(westvresident,eastvresident,migrantvresident)

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

#ggplot(newtop5, mapping=aes(BIN_START,migrant_vs_resident_Fst)) + geom_line() + labs(y="Migrant vs Resident Fst",x="Position") + xlim(0,1000000)
#plot1 = ggplot(newtop5, mapping=aes(BIN_START,migrant_vs_resident_Fst)) + geom_line() + labs(y="Migrant vs Resident Fst",x="Position") + xlim(0,1000000)

#Gene plot
data = read.table("danaus_plexippus_v3_core_32_85_1.gff")
t=c('DPSCF300098','DPSCF300147','DPSCF300267','DPSCF300067','DPSCF300240','DPSCF300322','DPSCF300095','DPSCF300353','DPSCF300170','DPSCF300127','DPSCF300331','DPSCF300114','DPSCF300105','DPSCF300047','DPSCF300236')
t

data2 = subset(data,data$V1==t)
genes = subset(data2,data2$V3=="gene")
genes = genes[order(genes$V4),]

#write.csv(genes,"Genes.csv")

table_sort = subset(table,table$CHROM==t)
#plot(genes$V4, factor(genes$V9))
genes$V9 = sub(";.*", "", genes$V9)  

Set_genes=read.csv("Genes.csv")
Set_genes = Set_genes[order(Set_genes$V1), ]
write.csv(Set_genes,"Genes.csv")

genes = read.csv("Genes.csv")
genes = genes[order(genes$chrom_pos), ]

#plot(genes$chrom_pos,genes$V9)

dotplot(genes$V9~genes$chrom_pos)

#ggplot(genes, mapping=aes(V4,V9)) +geom_line() + geom_point() + labs(y="Gene ID",x="Gene Start Position")
#plot2 = ggplot(genes, mapping=aes(V4,V9)) +geom_line() + geom_point() + labs(y="Gene ID",x="Gene Start Position")
# + geom_line(aes(x=V4,y=V9)) + geom_line(aes(x=V5,y=V9))


#unique(sorted_table2$CHROM)

#table2["chrom_pos"] <- "0"

#write.csv(table2, "Chr_16_table.csv")

#if(table2$CHROM=="DPSCF300047") {table2$chrom_pos == table2$BIN_START + 6619107}


chr_16=read.csv("Chr_16_table.csv")


chr_16 = chr_16[order(chr_16$chrom_pos), ]

plot(chr_16$chrom_pos ,chr_16$migrant_vs_resident_Fst)
lines(chr_16$chrom_pos ,chr_16$migrant_vs_resident_Fst)

mean1=mean(chr_16[["migrant_vs_resident_Fst"]])
mean(chr_16$migrant_vs_resident_Fst, na.rm = TRUE)  

ggplot(data=chr_16, aes(x=chrom_pos,y=migrant_vs_resident_Fst)) + geom_line() + labs(y="Migrant vs Resident Fst",x="Position on chromosome 16") + xlim(0,9000000) + theme_gray(base_size = 9)#+ geom_hline(yintercept = mean1,linetype="dashed", color = "blue")
plot1 = ggplot(data=chr_16, aes(x=chrom_pos,y=migrant_vs_resident_Fst)) + geom_line() + labs(y="Migrant vs Resident Fst",x="Position on chromosome 16") + xlim(0,9000000) + theme_gray(base_size = 9)#+ geom_hline(yintercept = mean1,linetype="dashed", color = "blue")

ggplot(data = genes, aes(x=chrom_pos, y= V9)) + geom_line() + geom_point() + labs(y="Gene ID",x="Gene Position on Chr 16") + xlim(0,9000000) + theme_gray(base_size = 9)
plot2 = ggplot(data = genes, aes(x=chrom_pos, y= V9)) + geom_line() + geom_point() + labs(y="Gene ID",x="Gene Position on Chr 16") + xlim(0,9000000) + theme_gray(base_size = 9)

cowplot::plot_grid(plot1, plot2, align = "v", ncol = 1, rel_heights = c(0.2, 0.8))
egg::ggarrange(plot1, plot2, heights = c(0.2, 0.8))


#v4=read.table("v4_annotation.txt")

#V4 ANNOTATION COMPARISON

v4=read.delim("v4_annotation.txt")

v4_chr16 = subset(v4,v4$Scaffold=="NC_045822.1")
v4_chr16_genes = subset(v4_chr16,v4_chr16$Function=="gene")
chr_16 = chr_16[order(chr_16$chrom_pos), ]

v4_chr16_genes$Gene_ID.s = sub(";.*", "", v4_chr16_genes$Gene_ID.s)
#ggplot(data = v4_chr16_genes, aes(x=Start_pos, y= Gene_ID.s)) + geom_line() + geom_point() + labs(y="Gene ID",x="Gene Position on Chr 16") + xlim(0,9000000) + theme_gray(base_size = 6)
plot3 = ggplot(data = v4_chr16_genes, aes(x=Start_pos, y= Gene_ID.s)) + geom_line() + geom_point() + labs(y="Gene ID",x="Gene Position on Chr 16") + xlim(0,9000000) + theme_gray(base_size = 6)

cowplot::plot_grid(plot1, plot3, align = "v", ncol = 1, rel_heights = c(0.2, 0.8))
egg::ggarrange(plot1, plot3, heights = c(0.2, 0.8))

filtered_v4 = subset(v4_chr16_genes,v4_chr16_genes$Start_pos>=7100000)
ggplot(data = filtered_v4, aes(x=Start_pos, y= Gene_ID.s)) + geom_line() + geom_point() + labs(y="Gene ID",x="Gene Position on Chr 16") + xlim(0,9000000) + theme_gray(base_size = 6)
plot4 = ggplot(data = filtered_v4, aes(x=Start_pos, y= Gene_ID.s)) + geom_line() + geom_point() + labs(y="Gene ID",x="Gene Position on Chr 16") + xlim(0,9000000) + theme_gray(base_size = 6)

cowplot::plot_grid(plot1, plot4, align = "v", ncol = 1, rel_heights = c(0.2, 0.8))
egg::ggarrange(plot1, plot4, heights = c(0.2, 0.8))

write.csv(filtered_v4, "Gene_List.csv")

#TEST FOR MEAN AND MEDIAN PLOTS
#unique(table$CHROM)
#mean_table=c("scaffold","mean")
#x = unique(table$CHROM)

#scaffolds = table$CHROM
#fst = table$migrant_vs_resident_Fst

#total = data.frame(scaffolds,fst)
#mean=total[ ,list(mean=mean(fst)), by=scaffolds]
#mean=aggregate(total$fst, list(total$scaffolds), FUN=mean)
#median=aggregate(total$fst, list(total$scaffolds), FUN=median)
#plot(mean$Group.1,mean$x)
#ggplot(data = mean, aes(x=Group.1, y=x)) + geom_point() + labs(y="Mean Fst",x="Scaffold_name") + theme_gray(base_size = 5)
#ggplot(data = median, aes(x=Group.1, y=x)) + geom_point() + labs(y="Median Fst",x="Scaffold_name") + theme_gray(base_size = 5)
plot1

length=read.delim("scaffolds_v3_lengths.txt")
plot(length$Scaffold,length$Length)
length2 = subset(length,length$Length>=100000)
plot99=ggplot(data = length, aes(x=Scaffold, y=Length)) + geom_point() + labs(y="Length",x="Scaffold") + geom_point(data = length2, aes(x='DPSCF300114', y=516792), color='red',size=5)#+ geom_line()
plot99

plot100 = ggplot(data = length2, aes(x=Scaffold, y=Length)) + geom_point() + geom_point(data = length2, aes(x='DPSCF300114', y=516792), color='red',size=5) + labs(y="Length",x="Scaffold") #+ geom_line()
plot100
