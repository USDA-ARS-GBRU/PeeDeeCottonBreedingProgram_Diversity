library(tidyverse)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed")

# get the names
IDnames <- read.table("unimputed.PD.mibs.id",
                       sep = "\t",
                       header = FALSE)$V2

# and the groups
groups <- read.table("../../Groups.txt")
names(groups) <- c("Group", "ID")

dataset1 <- read.table("unimputed.PD.mibs",
                       sep = "\t",
                       header = FALSE) %>%
  data.matrix(.)

dataset2 <- read.table("unimputed.PD.thinned.mibs",
                       sep = "\t",
                       header = FALSE) %>%
  data.matrix(.)

dimnames(dataset1) <- list(IDnames, IDnames)
dimnames(dataset2) <- list(IDnames, IDnames)

dataset1 <- dataset1 %>%
  as.data.frame(.) %>%
  rownames_to_column(., "ID") %>%
  pivot_longer(., 2:ncol(.), names_to = "ID2", values_to = "ibs")

dataset1 <- dataset1 %>%
  add_column(., Group = groups[match(dataset1$ID, groups$ID), "Group"]) %>%
  add_column(., Group2 = groups[match(dataset1$ID2, groups$ID), "Group"]) %>%
  filter(ID != ID2)

dataset2 <- dataset2 %>%
  as.data.frame(.) %>%
  rownames_to_column(., "ID") %>%
  pivot_longer(., 2:ncol(.), names_to = "ID2", values_to = "ibs")

dataset2 <- dataset2 %>%
  add_column(., Group = groups[match(dataset2$ID, groups$ID), "Group"]) %>%
  add_column(., Group2 = groups[match(dataset2$ID2, groups$ID), "Group"]) %>%
  filter(ID != ID2)

group_sort <- function(x, y) {
  return(paste0(sort(c(x,y)), collapse = ""))
}

dataset1 %>%
  mutate(., Unq = map2_chr(Group, Group2, group_sort)) %>%
  mutate(., GroupA = map_chr(Unq, function(x) substr(x,1,1))) %>%
  mutate(., GroupB = map_chr(Unq, function(x) substr(x,2,2))) -> dataset1

dataset2 %>%
  mutate(., Unq = map2_chr(Group, Group2, group_sort)) %>%
  mutate(., GroupA = map_chr(Unq, function(x) substr(x,1,1))) %>%
  mutate(., GroupB = map_chr(Unq, function(x) substr(x,2,2))) -> dataset2

temp_tab <- dataset1[, c("GroupA", "GroupB", "ibs")]

min(temp_tab$ibs)
max(temp_tab$ibs)

# aggregate to get means
ag <- aggregate(ibs ~ GroupA + GroupB, temp_tab, FUN = mean)

# turn into table with xtabs
xtabs(ibs~., data=ag)

temp_tab <- dataset2[, c("GroupA", "GroupB", "ibs")]

min(temp_tab$ibs)
max(temp_tab$ibs)

ag <- aggregate(ibs ~ GroupA + GroupB, temp_tab, FUN = mean)

xtabs(ibs~., data=ag)

# plot the correlation between datset1 and 2

# read in the GRM
dataset1 %>%
  add_column(., sorted = map2_chr(dataset1$ID, dataset1$ID2, 
                                  function(x, y) paste0(sort(c(x, y)), collapse = ""))) %>%
  distinct(., sorted, .keep_all = TRUE) -> k

dataset2 %>%
  add_column(., sorted = map2_chr(dataset2$ID, dataset2$ID2, 
                                  function(x, y) paste0(sort(c(x, y)), collapse = ""))) %>%
  distinct(., sorted, .keep_all = TRUE) -> k2

left_join(k, k2, by = c("ID" = "ID",
                        "ID2" = "ID2",
                        "Group" = "Group",
                        "Group2" = "Group2",
                        "GroupA" = "GroupA",
                        "GroupB" = "GroupB",
                        "sorted" = "sorted"),
          suffix = c(".dataset1", ".dataset2")) %>%
  select(., ID, ID2, Group, Group2, ibs.dataset1, ibs.dataset2) %>%
  write.table(.,
              "~/1workingdir/grant_trimfilter/paper1_reanalysis/results/ibs.txt",
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

s <- ggplot(data=dat, aes(x=dat1, y=dat2)) +
  xlab("Dataset1 - Observed Genetic Distance (IBS, plink 1.9)") +
  ylab("Dataset2 - Observed Genetic Distance (IBS, plink 1.9)") +
  geom_point(size=3) + 
  scale_color_brewer(palette="Spectral") + 
  theme_grey(base_size=60) +
  theme(axis.title = element_text(size = rel(.75))) +
  geom_smooth(method="lm", aes(color="red"), show.legend=FALSE)

s