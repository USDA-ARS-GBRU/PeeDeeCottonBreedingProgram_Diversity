library(RColorBrewer)
library(sm)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/")

set1 <- read.table("unimputed.PD.bim")[, 1:4]
set1$V1 <- as.integer(set1$V1)
names(set1) <- c("chr", "name", "x", "pos")

set2 <- read.table("unimputed.PD.thinned.bim")[, 1:4]
set2$V1 <- as.integer(set2$V1)
names(set2) <- c("chr", "name", "x", "pos")

set1 <- subset(set1, select = c("chr", "pos"))
set2 <- subset(set2, select = c("chr", "pos"))

set1$pos <- set1$pos/1000000
set2$pos <- set2$pos/1000000

par(mfrow = c(5, 6),       # 5 rows 6 columns
    oma = c(5, 5, 0, 0) + 0.1,   # outer margin 5 left 5 bottom
    mar = c(1,1,1,1) + 0.1)      # margin between cells 

# lookup the chromosome names and REORDER
chrSZ <- read.table("C:/Local/0GWAS/chrSZ.txt",
                    header = TRUE)

chrSZ[match(set1$chr, chrSZ$CHR_OLD), "CHR"] -> set1$chr
chrSZ[match(set2$chr, chrSZ$CHR_OLD), "CHR"] -> set2$chr

plot.new()

for (i in chrSZ$CHR) {
  # here we need to subset out two sets
  set1_sub <- subset(set1, set1$chr == i)$pos
  set2_sub <- subset(set2, set2$chr == i)$pos
  
  all <- as.numeric(c(set1_sub, set2_sub))
  wherefrom <- as.factor(c(rep("Dataset 1", length(set1_sub)),
                           rep("Dataset 2", length(set2_sub))))
  
  dens <- sm.density.compare(x = all, 
                             group = wherefrom, 
                             # h = 5,   # smooth param
                             lwd = 3, # line thickness
                             xlim = c(0, max(set1_sub)),  # set x min, max
                             model = "equal",
                             xlab = "", ylab = "",
                             col = brewer.pal(4, "PuOr")[2:3])
  if (dens$p < 0.05) {
    title(paste0("chr ",i," *"), line = -1, adj = 0.05)
  }
  else {
    title(paste0("chr ",i), line = -1, adj = 0.05)
  }
  
  print(dens$h)
}

title(xlab = "Chromosome Position (Mbp)",
      ylab = "Marker Density",
      outer = TRUE, 
      line = 3,
      cex.lab = 3)

for (i in 1:4) { plot(1, type = "n", axes=FALSE, xlab="", ylab="") }

legend(locator(1), levels(wherefrom),
       fill = brewer.pal(4, "PuOr")[2:3],
       cex = 2,
       xpd = NA)
