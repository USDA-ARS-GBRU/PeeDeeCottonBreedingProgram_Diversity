library(tidyverse)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/old_datasets/no_SNP_repeats_ID_repeats/")

line_names <- read.table("PD_dist.mibs.id", header = FALSE)$V2
mibs <- data.matrix(read.table("PD_dist.mibs", header = FALSE))

dimnames(mibs) <- list(line_names, line_names)

mibs_tab <- mibs %>% as.data.frame(.) %>%
  rownames_to_column(., "ID1") %>%
  pivot_longer(., 2:ncol(.), names_to = "ID2", values_to = "IBS")

# make a note of which ones to remove...
mibs_tab[(mibs_tab$IBS > 0.97) & (mibs_tab$ID1 != mibs_tab$ID2), ] %>%
  View(.)

write.table(cbind(c("PD", "PD"), c("FJA_AHK", "PD-113_AHK")),
            file = "PD_Dups.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# repeats the process for the world germplasm

line_names <- read.table("WORLD_dist.mibs.id", header = FALSE)$V2
mibs <- data.matrix(read.table("WORLD_dist.mibs", header = FALSE))

dimnames(mibs) <- list(line_names, line_names)

mibs_tab <- mibs %>% as.data.frame(.) %>%
  rownames_to_column(., "ID1") %>%
  pivot_longer(., 2:ncol(.), names_to = "ID2", values_to = "IBS")

# make a note of which ones to remove...
mibs_tab[(mibs_tab$IBS > 0.97) & (mibs_tab$ID1 != mibs_tab$ID2), ] %>%
  View(.)

write.table(cbind(c("PD", "PD"), c("FJA_AHK", "PD-113_AHK")),
            file = "PD_Dups.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)