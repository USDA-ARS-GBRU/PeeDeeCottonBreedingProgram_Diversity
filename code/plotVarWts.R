library(tidyverse)

var1 <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/unimputed.PD.eigenvec.var",
                   sep = "",
                   header = FALSE)
names(var1)[2] <- "SNP"

var2 <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/unimputed.PD.thinned.eigenvec.var",
                   sep = "",
                   header = FALSE)
names(var2)[2] <- "SNP"

# read in a bim file
bimFile <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/unimputed.PD.bim",
                      sep = "\t")[, c(2, 1, 4)]

names(bimFile) <- c("SNP", "CHR", "POS")

# read in chrSZ
sizes <- read.table("C:/Local/0GWAS/chrSZ.txt",
                    sep = "\t",
                    header = TRUE)

# resort the block map and give them the cumsums
sizes[order(sizes$CHR), ] -> sizes

sizes$ADDER <- c(0, cumsum(as.numeric(sizes$SZ))[1:(nrow(sizes) - 1)])
sizes$ADDER <- sizes$ADDER + cumsum(c(0, rep(12500000, 25)))

bimFile$CHR <- sizes[match(bimFile$CHR, sizes$CHR_OLD), "CHR"]
bimFile <- bimFile[order(bimFile$CHR, bimFile$POS), ]

bimFile$dummypos <- bimFile$POS + sizes[match(bimFile$CHR, sizes$CHR), 
                                        c("ADDER")]

temp <- sizes$ADDER
temp <- append(temp, sizes$SZ[length(sizes$SZ)] + temp[length(temp)])
sizes$ticksat <- (temp[1:length(temp) - 1] + temp[2:length(temp)])/2
rm(temp)

var1 %>%
  right_join(bimFile, .) -> var1

var2 %>%
  right_join(bimFile, .) -> var2



color <- sapply(var1$CHR, function(x) ifelse(match(x, sizes$CHR) %% 2 == 0, "grey", "black")) %>% unname(.)

var1 <- cbind(var1, "color" = color)

color <- sapply(var2$CHR, function(x) ifelse(match(x, sizes$CHR) %% 2 == 0, "grey", "black")) %>% unname(.)

var2 <- cbind(var2, "color" = color)

joined <- rbind(var1, var2)

joined <- joined %>% add_column(., Dataset = c(rep("Dataset1", nrow(var1)), 
                                               rep("Dataset2", nrow(var2))))

joined$V6 <- abs(joined$V6)

g <- ggplot(joined, aes(x = dummypos, y = V6, fill = Dataset, alpha = Dataset)) +
  geom_bar(stat = "identity", width = 5, color = NA, size = 0) +
  scale_fill_manual(values = c("purple", "#F56600")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_x_continuous(breaks = sizes$ticksat,
                     labels = sizes$CHR,
                     limits = c(0, max(joined$dummypos + 1000000))) +
  ylab("Contibution to PCA (Dimension 2)") +
  xlab("Chromosome") +
  coord_cartesian(ylim = c(0, 4)) +
  theme_dark() +
  theme(panel.background = element_rect(fill = "black", colour = "black",
                                        linetype = "solid", size =  0.5))

ggsave("C:/Local/temp.pdf", plot = g,
       width = 14, height = 7)

g <- ggplot(blockDF, aes(x = AVG, y = PIC)) + 
  geom_point(size = 1,
             show.legend = FALSE) +
  geom_segment(aes(x = MIN, y = PIC, xend = MAX, yend = PIC),
               lineend = "round",
               size = 1,
               show.legend = FALSE) +
  scale_color_manual(values = rep(c("darkgrey", "black"), 13)) +
  scale_x_continuous(breaks = sizes$ticksat,
                     labels = sizes$CHR,
                     limits = c(0, max(blockDF$MAX + 1000000))) +
  ylim(0, 1) + 
  ylab("Polymorphic Information Content (PIC)") +
  xlab("Chromosome")


# joined is what we're plotting now