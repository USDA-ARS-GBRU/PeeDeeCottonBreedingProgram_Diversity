library(tidyverse)

chrSZ <- read.table("C:/Local/0GWAS/chrSZ.txt",
                    header = TRUE,
                    sep = "\t")

bimfile <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/unimputed.PD.thinned.bim",
                      header = FALSE,
                      sep = "\t")

chrvec <- bimfile$V1

chrvec <- chrSZ[match(chrvec, chrSZ$CHR_OLD), "CHR"]

table(chrvec) -> a

bimfile <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/unimputed.PD.bim",
                      header = FALSE,
                      sep = "\t")

chrvec <- bimfile$V1

chrvec <- chrSZ[match(chrvec, chrSZ$CHR_OLD), "CHR"]


table(chrvec) -> b

tab <- cbind("chr_old" = chrSZ$CHR_OLD, "ds1" = b, "ds2" = a,
      "dif" = as.numeric(b) - as.numeric(a), perc = round(100*a/b,0))
# %>% 
  # write.table(.,
  #             "~/1workingdir/grant_trimfilter/paper1_reanalysis/results/markdens.txt",
  #             sep = "\t",
  #             quote = FALSE,
  #             col.names = NA)

tab %>% as.data.frame(.) %>% rownames_to_column(., "CHR") -> tab

tab$CHR %>% as.factor(.) -> tab$CHR

tab <- tab %>%
  pivot_longer(., ds1:ds2, names_to = "Dataset", values_to = "Number_SNPs")

ggplot(tab, aes(x = CHR, y = Number_SNPs, fill = Dataset)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.3)

# now normalize by MBp

tab_norm <- tab
tab_norm$SNP_dens <- tab_norm$Number_SNPs / (chrSZ[match(tab_norm$CHR, chrSZ$CHR), "SZ"]/1000000)                  

ggplot(tab_norm, aes(x = CHR, y = SNP_dens, fill = Dataset)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.75) +
           # fill = c(rep("#f1a340", 26), rep("#998ec3",26))) +
  scale_fill_manual(values = c("#f1a340", "#998ec3"), 
                    name = "Dataset", labels = c("Dataset 1", "Dataset2")) + 
  xlab("Chromosome") +
  ylab("SNP Density (Number of SNPs/Mb)") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size = 20))
  
