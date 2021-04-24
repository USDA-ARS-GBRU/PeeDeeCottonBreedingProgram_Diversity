library(tidyverse)
library(adegenet)
library(RColorBrewer)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/results/unimputed_fS_out-linear")

names <- read.table("../../data/Groups2.txt",
                    sep = "\t",
                    header = TRUE)

# read in the Q-matrix
# it appears that k=6 is appropropriate here

Q <- read.table("unimputed.PD.thinned.6.meanQ",
                header = FALSE,
                sep = "") %>%
  data.matrix(.) %>%
  round(., 3)

# reorder to match the predicted groups
Q <- Q[, c(4, 6, 5, 2, 3, 1)]

dimnames(Q) <- list(names$Line,
                    paste0("Pop ", 1:6))

# assign to grps
max_each <- colnames(Q)[max.col(Q, ties.method = "first")]
max_value <- Q[cbind(1:nrow(Q),
                     max.col(Q, ties.method = "first"))]

compo_df <- data.frame(Q, max_each, max_value,
                       check.names = FALSE)

compo_df <- compo_df[order(compo_df$max_each, -compo_df$max_value), ]

write.table(compo_df,
            "~/1workingdir/grant_trimfilter/paper1_reanalysis/results/fS_allInfo.txt", 
            sep = "\t", 
            quote = FALSE,
            col.names = NA)
# increase the plotter size to fit outliers
par(mar = c(6,5,5,2) + 0.1,
    oma = c(3,1,1,1), 
    lwd = 1)

compoplot(as.matrix(subset(compo_df, select = 1:(dim(compo_df)[2]-2))),
          col.pal = brewer.pal(dim(compo_df)[2]-2, "Spectral"),
          lab = row.names(compo_df),
          show.lab = TRUE,
          cleg = 1,
          posi = list(x = 0, y = 1.175),
          cex.names = .7,
          n.col=3,
          border = "Black")

write.table(cbind(as.vector(rownames(compo_df)),
                  compo_df$max_each,
                  names[match(rownames(compo_df), names$Line), "Group"]), 
           "LinesOrdered.txt",
           sep = "\t",
           quote = FALSE,
           col.names = FALSE,
           row.names = FALSE)
