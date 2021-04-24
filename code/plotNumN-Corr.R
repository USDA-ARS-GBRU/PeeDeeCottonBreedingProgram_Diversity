library(tidyverse)
library(ggrepel)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/")

unthinned <- read.table("NumN-unthinned.txt",
                        sep = "\t",
                        header = TRUE)

names(unthinned)[3] <- "GRM_unthinned"

thinned <- read.table("NumN-thinned.txt",
                        sep = "\t",
                        header = TRUE)

names(thinned)[3] <- "GRM_thinned"

dat <- full_join(unthinned, thinned)
# get the the group names

groups <- read.table("Groups.txt",
                     sep = "\t",
                     header = FALSE)[, c(2,1)]

names(groups) <- c("ID_OLD", "GROUP")
# groups$ID_OLD <- toupper(groups$ID_OLD)

renamer <- read.table("NameUpdates.txt",
                      sep = "\t",
                      header = TRUE)
names(renamer) <- c("ID_OLD", "ID")
# renamer$ID_OLD <- toupper(renamer$ID_OLD)

renamer <- renamer %>% left_join(., groups)

# fix the names in dat
dat$ID1 <- renamer[match(dat$ID1, renamer$ID_OLD), "ID"]
dat$ID2 <- renamer[match(dat$ID2, renamer$ID_OLD), "ID"]

dat <- dat %>%
  add_column(., Group1 = renamer[match(dat$ID1, renamer$ID), "GROUP"],
             .after = "ID2") %>%
  add_column(., Group2 = renamer[match(dat$ID2, renamer$ID), "GROUP"],
             .after = "Group1")


# get distincts

vapply(1:nrow(dat), function(x) paste0(sort(dat[x, 1:2]), collapse = ""), 
       FUN.VALUE = character(1)) %>%
  add_column(dat, combo = .) %>%
  distinct(., combo, .keep_all = TRUE) -> dat

# plot the thinned vs unthinned

s <- ggplot(data=dat, aes(x=GRM_unthinned, y=GRM_thinned)) +
  xlab("Dataset1 - Additive vanRaden Relationship Matrix") +
  ylab("Dataset2 - Additive vanRaden Relationship Matrix") +
  geom_point(size=3) + 
  scale_color_brewer(palette="Spectral") + 
  theme_grey(base_size=60) +
  theme(axis.title = element_text(size = rel(.75))) +
  geom_smooth(method="lm", aes(color="red"), show.legend=FALSE)

lm0 <- lm(GRM_thinned ~ GRM_unthinned, data = dat)
lm0
summary(lm0)
# plot the unthinned vs NumN

q <- ggplot(data=dat, aes(x=GRNM, y=GRM_unthinned)) +
  xlab("Coancestry, rxy, Calculated from Pedigrees (IBD, Numericware N)") +
  ylab("Dataset1 - Additive vanRaden Relationship Matrix") +
  geom_point(size=3) + 
  scale_color_brewer(palette="Spectral") + 
  theme_grey(base_size=60) +
  theme(axis.title = element_text(size = rel(.75))) +
  geom_smooth(method="lm", aes(color="red"), show.legend=FALSE)  +
  geom_label_repel(aes(label=ifelse(GRM_unthinned>1.5,
                                    as.character(paste0(ID1, " (", Group1, ") ",
                                                        ID2, " (", Group2, ")")),'')),
                   box.padding = .35, point.padding = .5,
                   segment.color = "black", size = 12, show.legend=FALSE) +
  geom_label_repel(aes(label=ifelse(GRNM > 1.25 & GRM_unthinned < 0.25,
                                    as.character(paste0(ID1, " (", Group1, ") ",
                                                        ID2, " (", Group2, ")")),'')),
                   box.padding = .35, point.padding = .5,
                   segment.color = "black", size = 12, show.legend=FALSE)
  # geom_label_repel(aes(label=ifelse(dat1<.53, 
  #                                   as.character(paste0(LineA, " (", GroupA, ") ",
  #                                                       LineB, " (", GroupB, ")")),'')), 
  #                  box.padding = .35, point.padding = .5, 
  #                  segment.color = "black", size = 12, show.legend=FALSE)

q

r <- ggplot(data=dat, aes(x=GRNM, y=GRM_thinned)) +
  xlab("Coancestry, rxy, Calculated from Pedigrees (IBD, Numericware N)") +
  ylab("Dataset2 - Additive vanRaden Relationship Matrix") +
  geom_point(size=3) + 
  scale_color_brewer(palette="Spectral") + 
  theme_grey(base_size=60) +
  theme(axis.title = element_text(size = rel(.75))) +
  geom_smooth(method="lm", aes(color="red"), show.legend=FALSE)  +
  geom_label_repel(aes(label=ifelse(GRM_thinned>1.5,
                                    as.character(paste0(ID1, " (", Group1, ") ",
                                                        ID2, " (", Group2, ")")),'')),
                   box.padding = .35, point.padding = .5,
                   segment.color = "black", size = 12, show.legend=FALSE) +
  geom_label_repel(aes(label=ifelse(GRNM > 1.25 & GRM_thinned < 0.25,
                                    as.character(paste0(ID1, " (", Group1, ") ",
                                                        ID2, " (", Group2, ")")),'')),
                   box.padding = .35, point.padding = .5,
                   segment.color = "black", size = 12, show.legend=FALSE)

lm1 <- lm(GRNM ~ GRM_unthinned, data = dat)
lm1
summary(lm1)

lm2 <- lm(GRNM ~ GRM_thinned, data = dat)
lm2
summary(lm2)
