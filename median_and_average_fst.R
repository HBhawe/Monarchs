# Load packages
library(ggplot2)
library(forcats)

# Load csv table
fst_values <- read.csv2(file = "migrant_vs_resident_Final_stats.csv", head = TRUE, sep = ",")

# Creating an empty data frame for the fst values as doubles
fst_as_doubles <- data.frame(migrant_vs_resident_Fst_double = double(), chr_number = integer(), scaffold_name = character())

# Adding the fst values as doubles to the new vector with chromosome numbers
for (i in 1:length(fst_values$migrant_vs_resident_Fst)) {
  number <- as.double(fst_values$migrant_vs_resident_Fst[i])
  fst_as_doubles[i,"migrant_vs_resident_Fst_double"] <- number
  chr_number <- fst_values$Chr_number[i]
  fst_as_doubles[i,"chr_number"] <- chr_number
  scaffold_name <- fst_values$CHROM[i]
  fst_as_doubles[i, "scaffold_name"] <- scaffold_name
}

# CHROMOSOMES

# Creating a data frame for mean and median over chromosomes
chr_mean_median <- data.frame(matrix(ncol = 4, nrow = 30))
colnames(chr_mean_median) <- c('Chr_num', 'Mean_Fst', 'Median_Fst', 'Max_Fst')


# Calculating the mean and median for each chromosome and adding the values into the chr_mean_median data frame
for (i in 2:31) {
  chr_mean <- mean(fst_as_doubles$migrant_vs_resident_Fst_double[fst_as_doubles$chr_number == i], na.rm = TRUE)
  chr_median <- median(fst_as_doubles$migrant_vs_resident_Fst_double[fst_as_doubles$chr_number == i], na.rm = TRUE)
  chr_max <- max(fst_as_doubles$migrant_vs_resident_Fst_double[fst_as_doubles$chr_number == i], na.rm = TRUE)
  chr_mean_median$Mean_Fst[i-1] <- chr_mean
  chr_mean_median$Median_Fst[i-1] <- chr_median
  chr_mean_median$Max_Fst[i-1] <- chr_max
  chr_mean_median$Chr_num[i-1] <- i
}

# Plotting mean, median and max against chromosome number
pdf("chr_mean_fst.pdf")
plot(chr_mean_median$Chr_num, chr_mean_median$Mean_Fst, type = 'b', main = "Mean Fst values", xlab = "Chromosome number", ylab = "Mean fst")
dev.off()
pdf("chr_median_fst.pdf")
plot(chr_mean_median$Chr_num, chr_mean_median$Median_Fst, type ='b', main = "Median Fst values", xlab = "Chromosome number", ylab = "Median fst")
dev.off()
pdf("chr_max_fst.pdf")
plot(chr_mean_median$Chr_num, chr_mean_median$Max_Fst, type = 'b', main = "Max Fst values", xlab = "Chromosome number", ylab = "Max fst")
dev.off()


# SCAFFOLDS

# Creating a data frame for mean and median over scaffolds
scaffold_mean_median <- data.frame(matrix(ncol = 5, nrow = length(unique(fst_as_doubles$scaffold_name))))
colnames(scaffold_mean_median) <- c('Chromosome_number', 'Scaffold_name', 'Mean_Fst', 'Median_Fst', 'Max_Fst')

unique_scaffolds <- c(unique(fst_as_doubles$scaffold_name))

# Calculating mean and median over each scaffold
for (i in 1:length(unique_scaffolds)) {
  scaffold_fst_list <- c()
  for (j in 1:length(fst_as_doubles$scaffold_name)) {
    if (identical(unique_scaffolds[i], fst_as_doubles$scaffold_name[j])) {
      scaffold_fst_list <- append(scaffold_fst_list, fst_as_doubles$migrant_vs_resident_Fst_double[j])
    }
    else {
      next
    }
  }
  scaffold_mean_fst <- mean(scaffold_fst_list, na.rm = TRUE)
  scaffold_median_fst <- median(scaffold_fst_list, na.rm = TRUE)
  scaffold_max_fst <- max(scaffold_fst_list, na.rm = TRUE)
  scaffold_mean_median$Mean_Fst[i] <- scaffold_mean_fst
  scaffold_mean_median$Median_Fst[i] <- scaffold_median_fst
  scaffold_mean_median$Max_Fst[i] <- scaffold_max_fst
  scaffold_mean_median$Scaffold_name[i] <- unique_scaffolds[i]
  scaffold_chr_number <- fst_as_doubles$chr_number[fst_as_doubles$scaffold_name == unique_scaffolds[i]]
  scaffold_mean_median$Chromosome_number[i] <- unique(scaffold_chr_number)
}

# Ordered by chromosomes
scaffold_mean_median_chr_ordered <- scaffold_mean_median[order(scaffold_mean_median$Chromosome_number),]

# Plotting mean, median and max against scaffold name
pdf(file = "scaffold_mean_fst.pdf")
plot(factor(scaffold_mean_median_chr_ordered$Scaffold_name), scaffold_mean_median_chr_ordered$Mean_Fst, type = 'b', main = "Mean Fst values", ylab = "Mean fst", xlab = "", na.rm = TRUE, las = 2, cex.axis = 0.3, cex.lab =0.8)
title(xlab = "Scaffold name", line = 4, cex.lab = 0.8)
dev.off()
pdf(file = "scaffold_median_fst.pdf")
plot(factor(scaffold_mean_median_chr_ordered$Scaffold_name), scaffold_mean_median_chr_ordered$Median_Fst, type ='b', main = "Median Fst values", ylab = "Median fst", xlab = "", las = 2, cex.axis = 0.3, cex.lab = 0.8)
title(xlab = "Scaffold name", line = 4, cex.lab = 0.8)
dev.off()
pdf(file = "scaffold_max_fst.pdf")
plot(factor(scaffold_mean_median_chr_ordered$Scaffold_name), scaffold_mean_median_chr_ordered$Max_Fst, type = 'b', main = "Max Fst values", ylab = "Max fst", xlab = "", las = 2, cex.axis = 0.3, cex.lab = 0.8)
title(xlab = "Scaffold name", line = 4, cex.lab = 0.8)
dev.off()

pdf(file = "scaffold_mean_ordered_chr.pdf")
ggplot(scaffold_mean_median_chr_ordered, aes(x = fct_inorder(factor(Scaffold_name)), Mean_Fst))
ggsave(filename = "scaffold_mean_ordered_chr.pdf", width = 20, height = 20)


# SCAFFOLDS ON CHROMOSOME 16

# Creating necessary data frames and subsets for chromosome 16
chr_16 <- subset(fst_as_doubles, fst_as_doubles$chr_number == 16)
unique_scaffolds_chr_16 <- c(unique(chr_16$scaffold_name))
chr_16_scaffolds_mean_median <- data.frame(matrix(ncol = 4, nrow = length(unique_scaffolds_chr_16)))
colnames(chr_16_scaffolds_mean_median) <- c('Scaffold_name', 'Mean_Fst', 'Median_Fst', 'Max_Fst')

# Calculating mean, median and max Fst on all scaffolds on chromosome 16
for (i in 1:length(unique_scaffolds_chr_16)) {
  chr_16_scaffold_fst_list <- c()
  for (j in 1:length(chr_16$scaffold_name)) {
    if (identical(unique_scaffolds_chr_16[i], chr_16$scaffold_name[j])) {
      chr_16_scaffold_fst_list <- append(chr_16_scaffold_fst_list, chr_16$migrant_vs_resident_Fst_double[j])
    }
    else {
      next
    }
  }
  chr_16_scaffold_mean_fst <- mean(chr_16_scaffold_fst_list, na.rm = TRUE)
  chr_16_scaffold_median_fst <- median(chr_16_scaffold_fst_list, na.rm = TRUE)
  chr_16_scaffold_max_fst <- max(chr_16_scaffold_fst_list, na.rm = TRUE)
  chr_16_scaffolds_mean_median$Mean_Fst[i] <- chr_16_scaffold_mean_fst
  chr_16_scaffolds_mean_median$Median_Fst[i] <- chr_16_scaffold_median_fst
  chr_16_scaffolds_mean_median$Max_Fst[i] <- chr_16_scaffold_max_fst
  chr_16_scaffolds_mean_median$Scaffold_name[i] <- unique_scaffolds_chr_16[i]
}

# Plotting mean, median and max on chromosome 16 against scaffold name
pdf(file = "chr_16_scaffold_mean_fst.pdf")
plot(factor(chr_16_scaffolds_mean_median$Scaffold_name), chr_16_scaffolds_mean_median$Mean_Fst, type = 'b', main = "Mean Fst values on chromosome 16", ylab = "Mean fst", xlab = "", na.rm = TRUE, las = 2, cex.axis = 0.6, cex.lab = 0.75)
title(xlab = "Scaffold name", line = 4, cex.lab = 0.75)
dev.off()
pdf(file = "chr_16_scaffold_median_fst.pdf")
plot(factor(chr_16_scaffolds_mean_median$Scaffold_name), chr_16_scaffolds_mean_median$Median_Fst, type ='b', main = "Median Fst values on chromosome 16", ylab = "Median fst", xlab = "", las = 2, cex.axis = 0.6, cex.lab = 0.75)
title(xlab = "Scaffold name", line = 4, cex.lab = 0.75)
dev.off()
pdf(file = "chr_16_scaffold_max_fst.pdf")
plot(factor(chr_16_scaffolds_mean_median$Scaffold_name), chr_16_scaffolds_mean_median$Max_Fst, type = 'b', main = "Max Fst values on chromosome 16", ylab = "Max fst", xlab = "", las = 2, cex.axis = 0.6, cex.lab = 0.75)
title(xlab = "Scaffold name", line = 4, cex.lab = 0.75)
dev.off()

