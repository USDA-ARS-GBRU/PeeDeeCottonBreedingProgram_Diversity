# now need to do DAPC
# with only dataset2
library(radiator)
library(adegenet)
library(RColorBrewer)
library(tidyverse)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/adegenet/")

# generate a genpop object
genpop_obj <- genomic_converter(data = "unimputed.PD.thinned.bed",
                                output = c("genind"),
                                verbose = TRUE)

dat <- genpop_obj$genind

rawfile <- "unimputed.PD.thinned.raw"
mapfile <- "unimputed.PD.thinned.map"

# need to manually reset group names
genlight_obj <- read.PLINK(file=rawfile, map.file = mapfile, parallel = FALSE,
                           quiet = FALSE)

# view the optimal number of PCs
dapc1 <- dapc.genlight(genlight_obj, n.pca = 75, n.da = 5)
optim.a.score(dapc1)

dapc1 <- dapc(genlight_obj)
summary.dapc(dapc1)
m <- cbind(genlight_obj$pop, genlight_obj$ind.names, dapc1$ind.coord)
m

inds <- list(genlight_obj$ind.names)
inds[[2]] <- c(.25)
inds[[3]] <- unlist(lapply(genlight_obj$pop, function(x) brewer.pal(8, "Spectral")[x]))
inds[[4]] <- c(1)

names(inds) <- c("labels", "air", "col", "cex")
dapc1$ind.coord[, 1] * -1 -> dapc1$ind.coord[, 1]

dapc1$grp.coord[, 1] * -1 -> dapc1$grp.coord[, 1] 

scatter.dapc(dapc1, 
             col=brewer.pal(8, "Spectral"), 
             label.inds=NULL,
             cex=3, 
             posi.da = "topright",
             scree.pca = TRUE, 
             posi.pca = "topleft",
             grid = TRUE,
             mstree= TRUE,
             # label.inds = inds
             )

# find some clusters!

outclust <- find.clusters.genlight(genlight_obj,
                                   max.n.clust=10)
outclust
plot(outclust$Kstat)

names(outclust$grp)
outclust_table <- data.frame("Line" = names(outclust$grp),
                             "DAPC_group" = paste("DAPC", c(2, 3, 1)[as.numeric(outclust$grp)]))
outclust_table
# rename groups to make more sense
renamer = data.frame("V1" = c(1, 2, 3), "V2" = c(2, 3, 1))
renamer <- c(2, 3, 1)
c(3, 1, 2)[as.numeric(outclust$grp)]

write.table(data.frame(pop(genlight_obj),
                       names(outclust$grp),
                       c(3, 1, 2)[as.numeric(outclust$grp)],
                       names(outclust$grp)),
            file = "C:/Users/grant/OneDrive/Documents/1workingdir/grant_trimfilter/paper1_reanalysis/results/dapcgrps.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote  = FALSE)

# now need to make the sansky table!!








