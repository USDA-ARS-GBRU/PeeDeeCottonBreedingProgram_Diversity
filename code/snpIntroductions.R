mat <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/full_datasets/unimputed/unimputed.PD.traw",
           header = TRUE,
           sep = "\t",
           check.names = FALSE)

rown <- mat$SNP

mat <- data.matrix(mat[, 7:ncol(mat)])
rownames(mat) <- rown

mat[is.na(mat)] <- 0

any(is.na(mat))

# check set all the 2't to 1's
mat[mat == 2] <- 1

mat %>%
  as.data.frame(.) %>%
  mutate_all(., as.logical) %>%
  as.matrix(.) -> mat

# true if contains the minor alleles
rownames(mat) <- rown

# import the group numbers

groups <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/Groups.txt",
                            header = FALSE,
                            sep = "\t",
                            check.names = FALSE)[, c(2, 1)]

names(groups) <- c("ID", "Group")

# now make the new table for summarizing
tab <- mat %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(., "ID") %>%
  right_join(groups, .)

tab$Group %>% as.factor(.) -> tab$Group
tab <- tab[, 2:ncol(tab)]

tab %>%
  group_by(., Group) %>%
  summarize_all(., any) -> tab_sum

tab_false <- tab
tab_false[, 2:ncol(tab_false)] <- ! tab_false[, 2:ncol(tab_false)]

tab_false %>%
  group_by(., Group) %>%
  summarize_all(., any) -> tab_sum2

# introduction of minor alleles
q <- sapply(2:ncol(tab_sum), 
       function(x) Position(isTRUE, tab_sum[, x, drop = TRUE]))

# introduction of major alleles
p <- sapply(2:ncol(tab_sum2), 
            function(x) Position(isTRUE, tab_sum2[, x, drop = TRUE]))


# where are the variants introduced in group eight?
rown[which(q == 8) - 1]

# i00706Gh
# i41792Gh
# i17702Gh
# i12411Gh
# i12409Gh
# i12408Gh
# i12406Gh
# i17700Gh
# i40932Gh
# i31241Gh
# i37866Gh
# i12405Gh
# i12404Gh
# i12401Gh
# i17697Gh
# i17696Gh
# i12396Gh
# i12395Gh
# i31240Gh
# i16780Gh
# i12695Gh

# View the table with the first introduction of each group
table(c(p, q))
# thinned
# 1    2    3    4    5    6    7    8 
# 8737  265   97   43   27   12    8    5 

# dataset1
# 1     2     3     4     5     6     7     8 
# 33341   831   335   182    60    86    25    21 

# write a table that says which group had minor/major introductions
out_df <- data.frame("SNP" = rown, 
                     "maj_intro" = p,
                     "min_intro" = q)

write.table(out_df,
            "~/1workingdir/grant_trimfilter/paper1_reanalysis/results/SNP_introductions.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# now check to see the number of SNPs FIXED in each group
# fixed means present in every individual
# tab %>%
#   group_by(., Group) %>%
#   summarize_all(., all) -> tab_all

# this is way faster
tab_all <- apply(tab[, 2:ncol(tab)], 
                 2, 
                 function(x) tapply(x, tab$Group, all)) %>%
  as.data.frame(.) %>%
  rownames_to_column(., "Group") %>%
  mutate(., Group = as.factor(Group))

# tab_false %>%
#   group_by(., Group) %>%
#   summarize_all(., all) -> tab_all2

tab_all2 <- apply(tab_false[, 2:ncol(tab_false)], 
                  2, 
                  function(x) tapply(x, tab_false$Group, all)) %>%
  as.data.frame(.) %>%
  rownames_to_column(., "Group") %>%
  mutate(., Group = as.factor(Group))

# maybe, figure out number of fixed alleles per group?
q2 <- apply(tab_all[, 2:ncol(tab_all)], 1, sum)
p2 <- apply(tab_all2[, 2:ncol(tab_all2)], 1, sum)

out_named <- rbind(p2 + q2)
colnames(out_named) <- 1:8

# thinned
# 1   2   3   4   5   6   7   8
# 457 764 854 702 816 561 569 273

# dataset 1
# 1    2    3    4    5    6    7   8
# 1539 2361 2590 2877 3398 2974 1572 792

q2
