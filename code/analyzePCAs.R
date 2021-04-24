# load all libraries
library("matrixStats") # for extracting matrix characteristics
library("snpStats")    # for dealing with genotype matrix

# needed for the DC-PCA scripts
library("ca")
library("schoolmath")
library("assertive")

# for graphing PCA
library("tidyverse")
library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library("SNPRelate")   # for phylogenetic analysis

source("C:/Users/grant/OneDrive/Documents/1workingdir/r/DC-PCA2.r")

setwd("C:/Users/grant/OneDrive/Documents/1workingdir/grant_trimfilter/paper1_reanalysis/data/")

folder_in <- "/full_datasets/unimputed/"

file_in <- "unimputed.PD"

fullrun_smt <- function(foldername, fileroot) {
  foldername <- folder_in
  fileroot <- file_in
  print(paste("fullrun function running for",foldername))
  
  setwd(paste0(getwd(), folder_in))
  # read in the single traw file
  geno <- read.table(paste0(getwd(), "/", fileroot, ".traw"),
                     header = TRUE,
                     check.names = FALSE)
  
  markers <- geno[, c("CHR", "SNP", "(C)M", "POS")]
  
  markers[markers$CHR == "*", "CHR"] <- NA
  markers[is.na(markers$CHR), "POS"] <- NA
  markers[is.na(markers$CHR), "(C)M"] <- NA
  
  # apply marker names to row names
  rownames(markers) <- markers$SNP
  
  # initialize an empty numeric vector
  markers_per_chr <- vector(mode = "integer", length = 26)
  
  # make the chr vector, count number of markers for each chr
  for (val in 1:26) {
    markers_per_chr[val] <- nrow(markers[markers$CHR == val 
                                        & !(is.na(markers$CHR)),])
  }
  
  # read in genotype matrix to geno
  
  # -----
  # load geno matrix
  
  # orderby sorted markers, linear order
  geno <- geno[, 7:ncol(geno)]
  rownames(geno) <- rownames(markers)
  
  # convert to a data matrix, add 1 to each cell to be of form 0/1/2
  #   +1 for minor allele count
  geno_matrix <- data.matrix(geno)
  
  # replace NA values with row averages
  indx <- which(is.na(geno_matrix), arr.ind = TRUE)
  
  geno_matrix[indx] <- as.integer(matrixStats::rowMedians(geno_matrix, na.rm = TRUE)[indx[,1]])
  
  ###### SKIP FOR JUST ASSOC
  # # write to file for DC-PCA, columns are indivs and rows are SNPS
  # # header is 96 GENO 16679 SNPS GRAD
  header <- paste(ncol(geno_matrix), "GENO",
                  nrow(geno_matrix), "SNPS",
                  "GRANT", sep = " ")
  
  write(header, file = "dcpca_in.txt")
  
  write.table(geno_matrix, file = "dcpca_in.txt", quote = FALSE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "", append = TRUE)
  #
  # # interactive
  PCA7(names(geno))
  
  print("")
  print("PCA done")
  
  # # how about imputation??
  # library(snpStats)
  # browseVignettes("snpStats")
  
  # transpose matrix
  # a <- t(a)
  
  # deal with phenotypes
  # -----
  
  # function for performing single marker analysis
  # we can use this to repeat the analysis across phenotypes
  # and these results can be fed into the manhattan plotter
  
  ## exit here
  return (NULL)
}

makegraph <- function(eigenvec, eigenval, program, foldername,
                      x_flip=1, y_flip=1) {
  
  # now pass to PCA plotter!
  # need to have ability to:
  #   choose outlier sections
  #   flip axises if necessary to match others
  
  # 
  # eigenvec <- as.data.frame(lapply(eigenvec, unlist))
  # 
  # eigenvec$Group <- lapply(eigenvec$Group, as.character)
  # eigenvec <- as.data.frame(lapply(eigenvec, unlist))
  # eigenvec$Group <- as.factor(eigenvec$Group)
  # eigenvec
  
  xlab <- paste("PC1 (", format(eigenval$ev[1]) , "%)", sep = "")
  
  ylab <- paste("PC2 (", format(eigenval$ev[2]), "%)", sep = "")
  
  # # change the outliers preferences here
  # 
  ifel <-  ifelse(eigenvec$PC2 > 100, as.character(eigenvec$line),'')
  
  
  #
  # ifel <- ifelse(rep(TRUE, length(eigenval$line)), as.character(eigenval$line), '')
  
  p <- ggplot(data = eigenvec, 
              aes(x = (x_flip)*PC1, y = (y_flip)*PC2, 
                  color = group, label)) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(paste(program, foldername, sep=" - ")) +
    geom_point(size=12) +
    scale_color_brewer(palette="Spectral") +
    theme_grey (base_size = 60) +
    theme(legend.position="right", aspect.ratio=1) +
    geom_label_repel(aes(label=ifel),   # label=ifel
                     box.padding = .35,
                     point.padding = .5,
                     segment.color = "black",
                     size = 18,
                     show.legend=FALSE)
  
  return(p)
  # 
  # #+ geom_text(aes(label=ifelse(eig2>.15,as.character(Line),'')), hjust=0, vjust=0)
  # 
  # p
  # 
}

# double check DC-PCA file root name
pca_plot <- function(foldername, fileroot) {
  foldername <- folder_in
  fileroot <- file_in
  # read in data for line names and group numbers
  groups <- read.table("Groups.txt", header = FALSE)
  
  groups <- subset(groups, select = c("V2", "V1"))
  names(groups) = c("line", "group")
  
  groups$group <- as.factor(groups$group)
  
  setwd(paste0(getwd(), folder_in))
  
  # read in data for plink eigenvec and eigenvals
  plink.eigenvec <- read.table(paste0(getwd(), "/PCA_res/", file_in, ".eigenvec"), 
                               header = FALSE)
  
  names(plink.eigenvec) <- c("group", "line", 
                             paste(mapply(paste, rep("PC", 40), 
                                          1:40, sep = "")))
  
  plink.eigenvec$group <- groups[match(plink.eigenvec$line, groups$line), "group"]
  
  plink.eigenval <- read.table(paste0(getwd(), "/PCA_res/", file_in, ".eigenval"),
                               header = FALSE)
  
  names(plink.eigenval) <- c("ev")
  
  plink.total <- sum(plink.eigenval$ev)
  plink.eigenval$ev <- round(sapply(plink.eigenval$ev, 
                                    function(i) i/plink.total*100), digits = 2)
  
  # read in data from dcpca
  dcpca.eigenvec <- read.table(paste0(getwd(), "/PCA_res/", "dcpca.thinned.eigenvec"),
                               header = TRUE)
  dcpca.eigenvec
  dcpca.eigenvec <- cbind(line = rownames(dcpca.eigenvec), dcpca.eigenvec)
  
  dcpca.eigenvec <- merge(dcpca.eigenvec, groups, by = "line")
  
  dcpca.eigenvec$group <- as.factor(dcpca.eigenvec$group)
  
  dcpca.eigenval <- read.table(paste0(getwd(), "/PCA_res/", "dcpca.thinned.eigenval"),
                               header = TRUE)
  
  dcpca.eigenval <- dcpca.eigenval[sapply(dcpca.eigenval$source,
                                          function(i) startsWith(as.character(i),
                                                                 "PC")),]
  dcpca.eigenval <- as.data.frame(dcpca.eigenval$ss)
  names(dcpca.eigenval) <- c("ev")
  dcpca.total <- sum(dcpca.eigenval$ev)
  dcpca.eigenval$ev <- round(sapply(dcpca.eigenval$ev, 
                                    function(i) i/dcpca.total*100), digits = 2)
  
  
  p <- makegraph(plink.eigenvec, plink.eigenval, "PLINK", foldername)
  
  print(p)
  
  while (TRUE) {
    flipx <- as.numeric(readline("flip x (1/-1) "))
    flipy <- as.numeric(readline("flip y (1/-1) "))
    
    p <- makegraph(plink.eigenvec, plink.eigenval, "PLINK", foldername, 
                   flipx, flipy)
    
    print(p)
    
    break
  }
  
  wait <- readline("enter to continue....")
  
  p <- makegraph(dcpca.eigenvec, dcpca.eigenval, "DC-PCA", foldername)
  
  print(p)
  
  while (TRUE) {
    flipx <- as.numeric(readline("flip x (1/-1) "))
    flipy <- as.numeric(readline("flip y (1/-1) "))
    
    p <- makegraph(dcpca.eigenvec, dcpca.eigenval, "DC-PCA", foldername,
                   flipx, flipy)
    
    print(p)
    
    break
  }
  
  return(NULL)
  
}

pca_plot("a", "b")
