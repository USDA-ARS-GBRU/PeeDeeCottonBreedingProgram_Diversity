library(tidyverse)

freq <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/group_freq/by_PDWORLD.frq.strat",
                   sep = "",
                   header = TRUE) %>%
  pivot_wider(., id_cols = "SNP", names_from = "CLST", values_from = "MAF")

freq %>% filter(., (PD == 0) | (PD == 1)) %>% nrow(.)/(nrow(freq)*2)*100

freq %>% filter(., (WORLD == 0) | (WORLD == 1)) %>% nrow(.)/(nrow(freq)*2)*100
