library(qqman)
library(calibrate)
library(radiator)
library(adegenet)
library(RColorBrewer)
library(tidyverse)

SNPs <- read.delim("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/unimputed.bim", header = FALSE, stringsAsFactors = FALSE)

SNPs <- SNPs[, c(2, 1, 4)]

names(SNPs) <- c("SNP", "chr", "pos")

# read in the bayesfactors 
BFs <- read.table("C:/Local/sweeps/filesmerged.txt",
                  header = FALSE,
                  sep = "\t",
                  colClasses = c("character", "character", "double"),
                  col.names = c("SNP", "x", "BF"))[, c("SNP", "BF")]

SNPs %>%
  left_join(., BFs) -> SNPs

# set the chromosome lengths
chr <- vector(mode = "numeric")
for (val in 1:26) {
  chr <- append(chr, nrow(SNPs[SNPs$chr == val,]))
}

# make a better plot

flatten_chr(sapply(1:26, function(x) rep(x, chr[x]))) %>%
  factor(., levels = 1:26) -> lvls 

# get the lengths
sizes <- read.table("C:/Local/0GWAS/chrSZ.txt",
                    sep = "\t",
                    header = TRUE)

# sizes[order(sizes$CHR_OLD), ] -> sizes
# As then Ds
sizes[order(sizes$CHR), ] -> sizes

sizes$ADDER <- c(0, cumsum(as.numeric(sizes$SZ))[1:(nrow(sizes) - 1)])
sizes$ADDER <- sizes$ADDER + cumsum(c(0, rep(12500000, 25)))

# reords SNPs
sizes[match(SNPs$chr, sizes$CHR_OLD), "CHR"] -> SNPs$chr
SNPs <- SNPs[order(SNPs$chr, SNPs$pos), ]

SNPs$dummypos <- SNPs$pos + sizes[match(SNPs$chr, sizes$CHR), 
                                        c("ADDER")]

temp <- sizes$ADDER
temp <- append(temp, sizes$SZ[length(sizes$SZ)] + temp[length(temp)])
sizes$ticksat <- (temp[1:length(temp) - 1] + temp[2:length(temp)])/2
rm(temp)

color <- sapply(SNPs$chr, function(x) ifelse(match(x, sizes$CHR) %% 2 == 0, "grey", "black")) %>% unname(.)

plot(log10(SNPs$BF) ~ SNPs$dummypos,
     col = color,
     # col = c(as.matrix(unlist(coldf))),
     cex = 0.75,
     pch = 16,
     # ylim = c(0, 12),
     xaxt = "n",
     ylab = expression(log[10](italic(BF))),
     xlab = "Chromosome",
     main = "Putative SNP Loci Under Selection, PD vs ALL World Improved Hirsutum")

abline(h = 1, col = "red", lty = "dashed", lwd = 1)

axis(side = 1, 
     at = sizes$ticksat, 
     labels = sizes$CHR,
     cex.axis = 0.66,
     xpd = NA,
     gap.axis = 0)

