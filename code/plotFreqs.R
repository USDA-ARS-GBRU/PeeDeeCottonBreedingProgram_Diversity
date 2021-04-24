library(tidyverse)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/group_freq/")

freq <- read.table("by_group.frq.strat", 
           sep = "",
           header = TRUE)

freq <- freq %>%
  pivot_wider(., id_cols = c(SNP, A1, A2),
              names_from = CLST,
              values_from = MAF)

names(freq) <- c("SNP", paste0("BreedingGroup_", 1:8), "WORLD")

# get the interesting blocks
intBlocks <- read.table("C:/Local/sweeps/interestingBlocks.txt")$V1

# get the SNPs
intSNPs <- read.table("C:/Local/sweeps/Windows.txt",
                      sep = "\t",
                      header = TRUE) %>%
  filter(., name %in% intBlocks)

# repair the chromosome names
chrSZ <- read.table("C:/Local/0GWAS/chrSZ.txt",
                    sep = "\t",
                    header = TRUE)
intSNPs$chr <- chrSZ[match(intSNPs$chr, chrSZ$CHR_OLD), "CHR"]

# now need to make the plots for this many SNPs
nrow(intSNPs)

freq_sub <- subset(freq, subset = SNP %in% intSNPs$snp)

# import SNP info
map <- read.table("unimputed.bim",
                 sep = "\t",
                 header = FALSE) %>%
  mutate(., V1 = vapply(V1, function(x) chrSZ[chrSZ$CHR_OLD == x, "CHR"],
                        FUN.VALUE = character(1))) %>%
  filter(., V2 %in% intSNPs$snp) %>%
  select(., 2, 1, 4) %>%
  rename(., SNP = 1, CHR = 2, POS = 3)

# go look in the other file to figure out how to plot these
# we use matplot

freq_sub <- map %>%
  left_join(., freq_sub)

# extract the table for this row
par(mfrow = c(2,3),
    mar = c(2,1,1,1),
    oma = c(5, 5, 1, 0) + 0.1,
    new = TRUE)

for (row in 1:6) {
  print(row)
  tab <- t(freq_sub[row, 6:14])
  tab <- cbind(1 - tab[,1], tab[,1])
  
  # make the name of the plot
  mname <- paste0("SNP:",
                  freq_sub[row, "SNP"],
                  " (chr ",
                  freq_sub[row, "CHR"],
                  " : ",
                  format(freq_sub[row, "POS"], big.mark = ","),
                  ")")
  
  matplot(tab,
          type = "b",
          pch = flatten_chr(freq_sub[row, c("A1", "A2")]),
          xlab = "",
          ylab = "",
          ylim = c(0,1),
          xaxt = "n",
          yaxt = "n",
          main =  mname)
  
  axis(side = 1, at = 1:9, c(1:8, "W"), xpd = NA)
  axis(side = 2, at = c(0, 0.5, 1), c(0, 0.5, 1), xpd = NA)
}

# set up the labels on outside
title(xlab = "Breeding Group",
      ylab = "Allele Frequency",
      outer = TRUE, 
      line = 3,
      cex.lab = 3)
