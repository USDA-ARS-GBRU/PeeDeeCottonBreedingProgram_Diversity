# install.packages(c("dplyr", "tidyverse", "RColorBrewer", "ggplot2", "ggrepel", "ape"))
                 
library("dplyr")
library("tidyverse")
library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library("ape")
library("reticulate")


#### put your file path here
nwk_file <-  "~/1workingdir/grant_trimfilter/paper1_reanalysis/results/mega/Original_Tree1.nwk"

vert.tree <- read.tree(file = nwk_file)

# make the naming changes in python

# convert the r character vector into a list
py$names <- reticulate::r_to_py(vert.tree$tip.label)

# trim off the population label
py_run_string('names = [i.replace("|population:PD", "") for i in names]')

# fix the naming of the other objects
py_run_string('names = [i.replace("AHK2", " (A_HK2)") for i in names]')
py_run_string('names = [i.replace("AHK1", " (A_HK1)") for i in names]')
py_run_string('names = [i.replace("AHK", " (A_HK)") for i in names]')
py_run_string('names = [i.replace("A_HK", "AHK") for i in names]')

# use implicit conversion from list to character vector
vert.tree$tip.label <- py$names

# round, convert to percentages
vert.tree$node.label <- round(as.numeric(vert.tree$node.label),2) * 100

# find edges with bootstrap values under threshold, collapse to 0 length
thresh <- 50

## do a for loop to loop through all the edges
step_ = 0
for (i in 1:nrow(vert.tree$edge)) {
  node1 <- vert.tree$edge[i,1]
  node2 <- vert.tree$edge[i,2]
  
  # skip if the edge length is associated with a leaf
  if (node1 <= length(vert.tree$tip.label) | 
      node2 <= length(vert.tree$tip.label)) {
    next
  }
  
  step_ = step_ + 1
  
  # need to skip first node label, since its an unrooted tree
  if (vert.tree$node.label[step_ + 1] < thresh) {
    # extract the node length
    node_len <- vert.tree$edge.length[i]
    
    # set the node length to zero
    vert.tree$edge.length[i] <- 0
    
    # add the node length to the children, to maintain overall distance
    vert.tree$edge.length <- sapply(1:nrow(vert.tree$edge),
                                    function(x) if (node2 == vert.tree$edge[x,1]) vert.tree$edge.length[x] + node_len else vert.tree$edge.length[x])
  }
  print(vert.tree$edge[i,])
  print(step_)
}

# inspect edge lengths
vert.tree$edge.length
di2multi(vert.tree, tol = 50)
# lookup colors, based on numbers
vert.tree$tip.label
# import the group names
names_groups <- read.table("~/1workingdir/grant_trimfilter/paper1_reanalysis/data/Groups2.txt",
                          header = TRUE,
                          sep = "\t")

# check if they are all present
setdiff(vert.tree$tip.label, names_groups$Line)

# set the text color
colors <- sapply(vert.tree$tip.label,
                 function(x) brewer.pal(8, "Spectral")[names_groups[names_groups$Line == x, "Group"]])


dev.off()
tiff("~/1workingdir/grant_trimfilter/paper1_reanalysis/results/img/phylo.tiff",
     units="in", width=18, height=6, res=300)

plot.new()
# set the background for the figure
par(bg = "white", 
    oma = c(1,2,1,0.1),
    mar = c(4,3,1,0.1), 
    xpd = NA)  # prevent clipping on bottom edge

plot(vert.tree,
     type = "phylogram",
     no.margin = FALSE, 
     edge.width = 2,
     # show.node.label = TRUE,
     show.tip.label = FALSE,
     direction = "downwards",
     tip.color = colors)

# filter the list of bootstrap values
bootstrap_lst <- sapply(vert.tree$node.label, function(x) if (is.na(x) | x < thresh) NA else x)

# get the correct bootstrap values
bootstrap_lab <- bootstrap_lst[!is.na(bootstrap_lst)]

# get the positions of the boostrap values
bootstrap_pos <- c(1:length(bootstrap_lst))[! is.na(bootstrap_lst)]

# skip the nodes that are leafs
bootstrap_pos <- bootstrap_pos + length(vert.tree$tip.label)

nodelabels(bootstrap_lab,
           node = bootstrap_pos,
           cex = 0.75)

#par(mar = c(0,0,0,0))

tiplabels(vert.tree$tip.label,
          col = colors,
          cex = 0.75,
          bg = "black",
          srt = -90,
          offset = .035,
          xpd = NA,
          adj = c(0.5, 0.5))

axisPhylo(side = 2, cex = 2)

title(ylab = "Proportion Polymorphic Sites",
      main = "Phylogenetic Tree for 114 Pee Dee Genotypes",
      xlab = "Genotype ID")

dev.off()
