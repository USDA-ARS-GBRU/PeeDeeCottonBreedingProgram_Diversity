library(tidyverse)

# read in the traw file
dat <- read.table("unimputed.PD.thinned.traw", sep = "\t",
                  header = TRUE,
                  check.names = FALSE)[, -c(1, 3:6)] %>%
  column_to_rownames(., var = "SNP") %>%
  data.matrix(.)

# count the number in each class
length(which(dat == 0))
length(which(dat == 1))
length(which(dat == 2))
length(which(is.na(dat)))

# set heterozygotes (1's) to NAs
dat[which(dat == 1)] <- NA

# make a isREF matrix
isREF <- matrix(data = FALSE, nrow = nrow(dat), ncol = ncol(dat),
                dimnames = dimnames(dat))
isREF[dat == 0] <- TRUE

# make an isALT matrix
isALT <- matrix(data = FALSE, nrow = nrow(dat), ncol = ncol(dat),
                dimnames = dimnames(dat))
isALT[dat == 2] <- TRUE

# read in the imputed data in order to calculate similarities;
# this might help in ordering
dat2 <- read.table("imputed.PD.thinned.traw", sep = "\t",
                  header = TRUE,
                  check.names = FALSE)[, -c(1, 3:6)] %>%
  column_to_rownames(., var = "SNP") %>%
  data.matrix(.)

dist(t(dat2)) %>%
  as.matrix(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(., "ID1") %>%
  pivot_longer(., cols = 2:ncol(.), names_to = "ID2",
               values_to = "Dist") %>%
  filter(., ID1 != ID2) %>%
  rowwise() %>%
  mutate(., concat = paste(sort(c(ID1, ID2)),  collapse = "_")) %>%
  ungroup() %>%
  distinct(., concat, .keep_all = TRUE) %>%
  select(., - concat) %>%
  arrange(., Dist) %>%
  as.data.frame(.) -> d2

head(d2, n = 100)

# need to find first instance of a true
# test <- isREF[1:nrow(isREF), 1:10]
# 
# apply(test, 1, FUN = function(x) Position(isTRUE, x))
# 
# rm(test)

# get the probability of discovering that allele
probALT <- apply(isALT, 1, mean)
probREF <- apply(isREF, 1, mean)

ceiling(1/probREF) %>% head

isALTREF <- rbind(isALT, isREF)
rownames(isALTREF)[1:nrow(isALT)] <- paste0(rownames(isALT), "_ALT")
rownames(isALTREF)[(nrow(isALT) + 1):nrow(isALTREF)] <- paste0(rownames(isREF), "_REF")

# probabilities
data.frame(probALT, probREF) %>%
  write.table(., file = "probabilities.txt",
              sep = "\t", quote = FALSE,
              col.names = NA)

# write the 0s matrix
write.table(isREF, file = "isREF.txt",
            sep = "\t", quote = FALSE,
            col.names = NA)

# write the 2s matrix
write.table(isALT, file = "isALT.txt",
            sep = "\t", quote = FALSE,
            col.names = NA)

# write the concat matrix
write.table(isALTREF, file = "isALTREF.txt",
            sep = "\t", quote = FALSE,
            col.names = NA)

probs <- c(probALT, probREF)
#calculate the poisson-binomial
library(poisbinom)

poisbinom::dpoisbinom(1:length(probs), probs) -> n
plot(x = 1:length(probs), y = n)
poisbinom::qpoisbinom(0.05, probs)

# do the prob that each one will not be covered:
num_covered <- vector(length = ncol(isALTREF), mode = "integer")
names(num_covered) <- 1:length(num_covered)

# loop thru each value of c and compute the distribution
for (c in 1:length(num_covered)) {
  probs_ <- (1 - probs)^c
  poisbinom::qpoisbinom(0.50, probs_) -> a
  num_covered[c] <- length(probs) - a
}

# plot the collectors curve
plot(x = 1:length(num_covered), 
     y = num_covered, 
     ylim = c(1, length(probs)),
     xlab = "number of individuals examined",
     ylab = "number of SNP alleles discovered")

abline(h = 0.8 * length(probs), col = "blue")
abline(h = 0.95 * length(probs), col = "red")

as.data.frame("percent alleles discovered", num_covered/max(num_covered))

# perform a similar calculation, but sampling within each breeding group first

# get groups
groups <- read.table("C:/Users/grant/OneDrive/Documents/1workingdir/grant_trimfilter/paper1_reanalysis/data/Groups.txt",
                     sep = "\t",
                     header = FALSE)[,c(2,1)]

# reorder groups
groups <- groups[match(colnames(isALTREF), groups$V2), "V1"] %>% as.factor

# calculate the groupwise allele probabilities
vapply(1:nrow(isALTREF), function(x) tapply(isALTREF[x, ], groups, mean),
       FUN.VALUE = vector(mode = "numeric", length = 8)) %>%
  t(.) %>%
  data.frame(., row.names = rownames(isALTREF)) %>%
  setNames(., 1:8) -> probs2

# estimate the number of alleles discovered if you alternate taking
# one individual from each until you run out

# build vector
takeOrder <- vector(mode = "integer", length = length(groups))
groups2 <- as.integer(groups)
ptr <- 1
while (ptr <= length(groups)) {
  # get the unique numbers and sort them
  insert <- sort(unique(na.omit(groups2)))
  
  # add the number numbers to the ordering vector
  takeOrder[ptr:(length(insert) + ptr - 1)] <- insert
  
  # move the pointer
  ptr <- ptr + length(insert)
  
  # remove these numbers from the collection if theyve been used
  groups2 <- groups2[-match(insert, groups2)]
}

rm(groups2)
rm(insert)
rm(ptr)

# now need to go through and sample in order

# this part is the same as above
num_covered2 <- vector(length = ncol(isALTREF), mode = "integer")
names(num_covered2) <- 1:length(num_covered2)

# loop thru each value of c and compute the distribution
# length(num_covered)
# loop through each individual
for (c in 1:114) {
  # select current and all previous columns
  probs_temp <- 1 - probs2[, takeOrder[1:c], drop = FALSE]
  
  # reduce each row to a single value, the product of 
  # all of the probabilities of not finding the allele
  probs_ <- apply(probs_temp, 1, function(x) Reduce("*", x))

  # improvement: save the probabilities from the previous iteration
  rm(probs_temp)
  poisbinom::qpoisbinom(0.50, probs_) -> a2
  num_covered2[c] <- length(probs) - a2
}

# get the best case ordering
# to get 80% cv the fastest
core80 <- read.table("C:/Users/grant/OneDrive/Documents/1workingdir/collectors_curve/Genocore-master/PDOut_80_Coreset.csv",
                     sep = ",",
                     header = TRUE,
                     check.names = FALSE) %>% colnames(.)

# to get 95% cv the fastest
core95 <- read.table("C:/Users/grant/OneDrive/Documents/1workingdir/collectors_curve/Genocore-master/PDOut_95_Coreset.csv",
                     sep = ",",
                     header = TRUE,
                     check.names = FALSE) %>% colnames(.)

# to get 99% cv the fastest
core99 <- read.table("C:/Users/grant/OneDrive/Documents/1workingdir/collectors_curve/Genocore-master/PDOut_99_Coreset.csv",
                     sep = ",",
                     header = TRUE,
                     check.names = FALSE) %>% colnames(.)

# calculate percent coverage OUR way for plotting
# for each coreset, calculate the number of SNP alleles discovered for each "round"
core80_disc <- apply(isALTREF[, core80], 1, 
                     function(x) Position(isTRUE, as.logical(x), nomatch = 10000))

core95_disc <- apply(isALTREF[, core95], 1, 
                     function(x) Position(isTRUE, as.logical(x), nomatch = 10000))

core99_disc <- apply(isALTREF[, core99], 1, 
                     function(x) Position(isTRUE, as.logical(x), nomatch = 10000))

# calculate the number of alleles discovered at each step
core80_nsnp <- sapply(1:length(core80), function(x) sum(core80_disc < (x + 1)))
core95_nsnp <- sapply(1:length(core95), function(x) sum(core95_disc < (x + 1)))
core99_nsnp <- sapply(1:length(core99), function(x) sum(core99_disc < (x + 1)))

# check how close these values are to the ones we expect
tail(core80_nsnp, n = 1) / nrow(isALTREF)
tail(core95_nsnp, n = 1) / nrow(isALTREF)
tail(core99_nsnp, n = 1) / nrow(isALTREF)

# check where num_covered(2) tails off at 99%
nc_99 <- Position(isTRUE, num_covered/nrow(isALTREF) > .99)
nc2_99 <- Position(isTRUE, num_covered2/nrow(isALTREF) > .99)
# calculate the exponential for num_covered
# library(drc)
# library(aomisc)
# devtools::install_github("OnofriAndreaPG/aomisc")
# mod <- drm(y ~ x, fct = DRC.asymReg(),
#            data = data.frame("x" = 1:length(num_covered),
#                              "y" = num_covered))
# plot(mod, log="")

# add the new points to the old graph
plot(x = 1:length(num_covered), 
     y = num_covered, 
     # ylim = c(1, length(probs)),
     xlim = c(0, nc_99 + 1),
     ylim = c(min(num_covered), max(num_covered)),
     xlab = "number of individuals examined",
     ylab = "number of SNP alleles discovered",
     col = "white",
     pch = 2)

# the line is plotted in random order
lines(x = 1:length(num_covered), 
      y = num_covered)

# plot the points for the breeding group discovery order
# filter to STOP plotting at 99%
points(x = 1:nc2_99,
       y = num_covered2[1:nc2_99],
       col = "blue",
       pch = 3)

# 99% cv
points(x = 1:length(core99_nsnp),
       y = core99_nsnp,
       col = "purple",
       pch = 4)

# 95% cv
points(x = 1:length(core95_nsnp),
       y = core95_nsnp,
       col = "red",
       pch = 4)

# plot the points for the 80% cv group
points(x = 1:length(core80_nsnp),
       y = core80_nsnp,
       col = "darkgreen",
       pch = 4)

abline(h = 0.8 * length(probs), col = "black", lty = "dashed")
abline(h = 0.95 * length(probs), col = "black", lty = "dotted")
abline(h = 0.99 * length(probs), col = "black", lty = "twodash")

# find the best combination that maximizes diversity
combn(x = 1:ncol(isALTREF), m = 5) -> combos_

# iterate thru the columns in combn
com <- combos_[, 2]

out_vec <- vector(length = ncol(combos_),
                  mode =  "integer")

z <- ncol(combos_)
print(z)
z <- 5000

# iterate thru the columns in combos_
for (col in 1:z) {
  ij <- combos_[, col]
  out_vec[col] <- sum(rowSums(isALTREF[, ij]) > 0L)
}

# write the file with just the genotypes from GWas
# read in the list of genotypes...
pheno_list <- read.table("C:/Local/0GWAS/plink/plink.fam", 
                         sep = "", header = FALSE)[, 2]

pheno_list[14] <- "PD2164_AHK"
pheno_list

#####
# read in another file with the thinned, imputed genotypes
# read in the imputed data in order to calculate similarities;
# this might help in ordering
dat3 <- read.table("C:/Local/0GWAS/plink/thinned.traw", sep = "\t",
                   header = TRUE,
                   check.names = FALSE)[, -c(1, 3:6)] %>%
  column_to_rownames(., var = "SNP") %>%
  data.matrix(.)

dat3 %>%
  write.table(., 
              file = "phenoPD_imput.csv", 
              sep = ",", col.names = NA, quote = FALSE)

dat3_nohet <- dat3
dat3_nohet[dat3_nohet == 1] <- NA

dat3_nohet %>%
  write.table(., 
              file = "phenoPD_imput_nohet.csv", 
              sep = ",", col.names = NA, quote = FALSE)


# now try and find a coreset of SNPs
# replace NAs
dat2_nohet_nona <- dat2_nohet
dat2_nohet_nona[is.na(dat2_nohet_nona)] <- 0

### ALGORITHM
### find most informative SNP (highest MAF)
### find the most different SNP -> most pairwise differences
### add that row
### remove unique individuals 
### repeat

# transpose the mat
working_mat <- t(dat2_nohet_nona)

# WHAT IF WE DO INDIVIDUALS? changes this line :)
# working_mat <- dat2_nohet_nona

# convert into a 1's and 0's matrix
working_mat[working_mat == 2] <- 1

building_mat <- matrix(data = 0L, nrow = nrow(working_mat), ncol = 0,
       dimnames = list(rownames(working_mat), NULL))

# intialize the distance matrix....
# columns are SNPs in building_mat, rows are SNPs that are to-be examined
dist_mat <- matrix(data = 0L, nrow = ncol(working_mat), ncol = 0,
                   dimnames = list(colnames(working_mat), NULL))

## INITIAILIZE
# calculate the MAFs
cSums <- colSums(working_mat)

# pick the current snp
this.SNP <- colnames(working_mat)[which.max(cSums)]

while (nrow(working_mat) > 1) {
  # update status
  print(this.SNP)
  print(dim(building_mat))
  
  # cbind the new SNP to the building_mat
  building_mat <- cbind(building_mat, working_mat[, this.SNP])
  colnames(building_mat)[ncol(building_mat)] <- this.SNP
  
  # drop this column from working_mat
  working_mat <- subset(working_mat, select = -c(get(this.SNP)))
  
  # drop this ROW from dist_mat
  dist_mat <- dist_mat[-c(which(rownames(dist_mat) %in% this.SNP)),]
  
  # calculate the distance_matrix for this SNP, and place into dist_mat
  
  # extract the new vector
  a_vec <- building_mat[, ncol(building_mat)]

  # now, attach the distance vectors to dist_mat
  dist_mat <- cbind(dist_mat, 
                    apply(working_mat, 2, 
                          function(x) dist(rbind(a_vec, x))))
  
  # and name it!
  colnames(dist_mat)[ncol(dist_mat)] <- this.SNP
  
  # remove individuals that can be uniquely identified;
  # adding columns can't change that!
  # https://stackoverflow.com/questions/38142890/find-unique-rows-in-a-data-frame-in-r
  
  # find unique by searching with duplicated up, then down....
  tf <- !(duplicated(building_mat) | duplicated(building_mat, 
                                                fromLast = TRUE))
  
  if (any(tf)) {
    # then drop that row!
    building_mat <- building_mat[!tf, ]
    working_mat <- working_mat[!tf, ]
    
    # should really recompute distance between SNPs here...
    # but that would really slow things down...
    # lets see what we get without doing that FIRST!
    
  }

  # now calculate rowSums, and find the MAX; this is the next THIS.SNP
  # this is the SNP that is most dissimilar to those we're already looking at
  this.SNP <- names(which.max(rowSums(dist_mat)))
}

# look at MAFs
dat2_nohet_nona[colnames(building_mat),] %>% rowMeans -> mafs

# now... lets check to see these are actually unique....
building_mat

t(dat2_nohet_nona[colnames(building_mat), ]) %>%
  apply(., 1, paste0, collapse = "") -> concated

# write the list of SNPs here
write.table(colnames(building_mat), file = "UniqueSNPs.txt", 
            sep = "\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE)

# test these data against cottongen

# import the cottongen dataset
dat.cg <- read.table("unimputed.noPD.traw", sep = "\t",
                  header = TRUE,
                  check.names = FALSE)[, -c(1, 3:6)] %>%
  column_to_rownames(., var = "SNP") %>%
  data.matrix(.)

dat.cg[is.na(dat.cg)] <- 0
dat.cg[dat.cg == 1] <- 0
dat.cg[dat.cg == 2] <- 1

# now test to see how many non-duplicated entries there are!
dat.cg.sub <- t(dat.cg[rownames(dat.cg) %in% colnames(building_mat), ])

dat.cg.sub
apply(dat.cg.sub, 1, paste0, collapse = "") -> cg.concat


tf <- !(duplicated(cg.concat) | duplicated(cg.concat, 
                                              fromLast = TRUE))

# okay, what about optimizing a different way?
# lets try to do it by using correlations....

# transpose the mat
working_mat <- t(dat2_nohet_nona)

# convert into a 1's and 0's matrix
working_mat[working_mat == 2] <- 1

# filter out low freq variants, or those with only HETs
working_mat <- working_mat[, which(colSums(working_mat) > 0.025*nrow(working_mat))]

building_mat <- matrix(data = 0L, nrow = nrow(working_mat), ncol = 0,
                       dimnames = list(rownames(working_mat), NULL))

# # intialize the distance matrix....
# # columns are SNPs in building_mat, rows are SNPs that are to-be examined
# dist_mat <- matrix(data = 0L, nrow = ncol(working_mat), ncol = 0,
#                    dimnames = list(colnames(working_mat), NULL))

library(DescTools)
combo_lst <- t(combn(1:ncol(working_mat), 2))
dim(combo_lst)

combo_sub <- combo_lst[,] %>% as.data.frame
colnames(combo_sub) <- c("x", "y")

combo_sub <- combo_sub[1:100000, ]
cortab <- map2_dbl(combo_sub[,1], combo_sub[,2], 
         function (x, y) cor(working_mat[, x],
                                               working_mat[, y]))

# instead, lets find the correlation between every combination of two columns
cor(working_mat[, 1], working_mat[, 10])

## INITIAILIZE
# calculate the MAFs
cSums <- colSums(working_mat)

# pick the current snp
this.SNP <- colnames(working_mat)[which.max(cSums)]

while (nrow(working_mat) > 1) {
  # update status
  print(this.SNP)
  print(dim(building_mat))
  
  # cbind the new SNP to the building_mat
  building_mat <- cbind(building_mat, working_mat[, this.SNP])
  colnames(building_mat)[ncol(building_mat)] <- this.SNP
  
  # drop this column from working_mat
  working_mat <- subset(working_mat, select = -c(get(this.SNP)))
  
  # drop this ROW from dist_mat
  dist_mat <- dist_mat[-c(which(rownames(dist_mat) %in% this.SNP)),]
  
  # calculate the distance_matrix for this SNP, and place into dist_mat
  
  # extract the new vector
  a_vec <- building_mat[, ncol(building_mat)]
  
  # now, attach the distance vectors to dist_mat
  dist_mat <- cbind(dist_mat, 
                    apply(working_mat, 2, 
                          function(x) dist(rbind(a_vec, x))))
  
  # and name it!
  colnames(dist_mat)[ncol(dist_mat)] <- this.SNP
  
  # remove individuals that can be uniquely identified;
  # adding columns can't change that!
  # https://stackoverflow.com/questions/38142890/find-unique-rows-in-a-data-frame-in-r
  
  # find unique by searching with duplicated up, then down....
  tf <- !(duplicated(building_mat) | duplicated(building_mat, 
                                                fromLast = TRUE))
  
  if (any(tf)) {
    # then drop that row!
    building_mat <- building_mat[!tf, ]
    working_mat <- working_mat[!tf, ]
    
    # should really recompute distance between SNPs here...
    # but that would really slow things down...
    # lets see what we get without doing that FIRST!
    
  }
  
  # now calculate rowSums, and find the MAX; this is the next THIS.SNP
  # this is the SNP that is most dissimilar to those we're already looking at
  this.SNP <- names(which.max(rowSums(dist_mat)))
}






