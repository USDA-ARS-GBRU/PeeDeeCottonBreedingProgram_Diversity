# need to make a GRNM
library("snpReady")
setwd("C:/Users/grant/OneDrive/Documents/1workingdir/grant_trimfilter/paper1_reanalysis/data/")

NumN <- read.table("NumericwareN.txt",
                   check.names = FALSE)


geno <- read.table(paste0("./full_datasets/imputed/imputed.PD.traw"),
                   header = TRUE,
                   check.names = FALSE)

markers <- geno[, c("CHR", "SNP", "(C)M", "POS")]

markers[markers$CHR == "*", "CHR"] <- NA
markers[is.na(markers$CHR), "POS"] <- NA
markers[is.na(markers$CHR), "(C)M"] <- NA

# apply marker names to row names
rownames(markers) <- markers$SNP

# orderby sorted markers, linear order
geno <- geno[, 7:ncol(geno)]
rownames(geno) <- rownames(markers)

toupper(names(geno))[which(!(toupper(names(geno)) %in% NumN$ID1))] %>% 
  write.table("plink_to_nn.txt",
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)

map_over <- read.table("plink_to_nn.txt",
                       header = TRUE)

for (row in 1:nrow(map_over)) {
  # first, add the column
  new_name <- map_over[row, "NEW"]
  old_name <- map_over[row, "OLD"]
  new_column <- NumN[, old_name]
  add_column(NumN, !!as.symbol(new_name) := new_column) -> NumN
  
  # make the new row vector
  new_row <- vector(mode = "numeric", length = ncol(NumN))
  new_row[1:(ncol(NumN)-1)] <- NumN[, new_name]
  
  # add the self distance
  new_row[ncol(NumN)] <- NumN[old_name, old_name]
  
  # grab the names
  names(new_row) <- colnames(NumN)
  
  # add it to the matrix
  rbind(NumN, new_row) -> NumN
  
  # set the new row name
  rownames(NumN)[nrow(NumN)] <- new_name
}

NumN[toupper(names(geno)), toupper(names(geno))] -> NumN

write.table(NumN,
            file = "NumericwareN_new.txt",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE
)

NumN <- data.matrix(read.table("NumericwareN_new.txt",
                               sep = "",
                               header = TRUE,
                               check.names = FALSE))

# make the vanRaden matrix
kk_vR_method1 <- snpReady::G.matrix(t(geno), method = c("VanRaden"),
                          format = "wide")

G.additive <- kk_vR_method1$Ga

cor(c(kk_vR_method1$Ga), c(NumN))

# adjust G.additive

A <- NumN
G <- G.additive

Beta = (mean(diag(A)) - mean(A))/(mean(diag(G)) - mean(G))
Alpha = mean(A) - mean(G)*Beta
# Alpha = mean(diag(A)) - mean(diag(G)) * Beta

Gs <- Beta * G + Alpha

G.additive <- Gs

G.additive.nas <- G.additive
diag(G.additive.nas) <- NA

NumN.nas <- NumN
diag(NumN.nas) <- NA

# cor(c(G.additive.nas), c(NumN.nas), use = "complete.obs")

# make tables foreach
G.additive.nas.long <- G.additive.nas %>%
  as.data.frame(.) %>%
  rownames_to_column(., "ID1") %>%
  pivot_longer(., 2:ncol(.), names_to = "ID2", values_to = "val")

NumN.nas.long  <- NumN.nas %>%
  as.data.frame(.) %>%
  rownames_to_column(., "ID1") %>%
  pivot_longer(., 2:ncol(.), names_to = "ID2", values_to = "val")

G.additive.nas.long <- G.additive.nas.long %>% drop_na(.)
NumN.nas.long <- NumN.nas.long %>% drop_na(.)

names(G.additive.nas.long)[3] <- "GRM"

names(NumN.nas.long)[3] <- "GRNM"

joined <- G.additive.nas.long %>%
  add_column(., GRNM = NumN.nas.long$GRNM)

plot(joined$GRNM ~ joined$GRM)

summary(lm(joined$GRNM ~ joined$GRM))

# try again for the unthinned dataset....
write.table(joined, 
            file = "NumN-unthinned.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
