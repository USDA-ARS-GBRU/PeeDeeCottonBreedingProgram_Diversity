library(SNPRelate)
setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/")

snpgdsBED2GDS(bed.fn = "unimputed.thinned.bed", 
              fam.fn = "unimputed.thinned.fam", 
              bim.fn = "unimputed.thinned.bim", 
              out.gdsfn = "unimputed.thinned.gds",
              cvt.chr = "char",
              verbose = TRUE)

snpgdsSummary("unimputed.thinned.gds")

genofile <- openfn.gds("unimputed.thinned.gds")

datset1 <- snpgdsSNPRateFreq(genofile)$MinorFreq

closefn.gds(genofile)

groups <- read.table("./unimputed.thinned.fam")
groups <- subset(groups, select = c("V1", "V2"))
names(groups) <- c("group", "indiv")


ibsMatrix <- snpgdsIBS(genofile, sample.id = NULL, snp.id = NULL, autosome.only=FALSE, remove.monosnp=FALSE, maf = NaN, missing.rate = NaN, num.thread = 1, verbose = TRUE)

hc <- snpgdsHCluster(ibsMatrix, sample.id=NULL, need.mat=TRUE, hang=-1)

cutTree <- snpgdsCutTree(hc, z.threshold=5, outlier.n=2, n.perm = 1000000, samp.group=factor(groups$group), col.outlier="red", col.list=c("red", "black"), pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, verbose=TRUE)

par(mar = c(2,5,6,2) + 0.1, lwd = 1, xpd = TRUE)
plot.new()

plot(cutTree$dendrogram, 
     #leaflab="perpendicular", 
     leaflab = "none",
     #main="Dendrogram of 96 Pee Dee Genotypes",
     main = "Dendrogram of Pee Genotypes vs Other Improved Hirsutum",
     ylab="1-IBS Genetic Distance",
     cex=12)

legend("topright", 
       legend = c("PD", "Other World"), 
       title = "Classificiation",
       col = c("red", "black"),
       box.col = "black", box.lty=1, box.lwd=1,
       pch = 20,  pt.cex = 1.5, cex = 1 , 
       text.col = "black", horiz = FALSE, inset = c(-0.01, -.1))

