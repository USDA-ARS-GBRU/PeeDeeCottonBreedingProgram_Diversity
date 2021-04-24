library(tidyverse)

setwd("~/1workingdir/grant_trimfilter/paper1_reanalysis/")


# import the files with the groups
orig_grps <- read.table("./data/Groups.txt",
                        sep = "\t",
                        header = FALSE,
                        col.names = c("BreedingGroup", "ID"))[, c(2, 1)]

fs_grps <- read.table("./results/fsgrps.txt",
                      sep = "\t",
                      header = FALSE,
                      col.names = c("ID", "fsGroup", "x"))[, c(1,2)]

dapc_grps <- read.table("./results/dapcgrps.txt",
                        sep = "\t",
                        header = FALSE,
                        col.names = c("a", "b", "DAPCGroup", "ID"))[, c(4,3)]

# repair the names
dapc_grps$DAPCGroup <- dapc_grps$DAPCGroup %>%
  vapply(., function(x) paste("DAPC Group", x), FUN.VALUE = character(1))

orig_grps$BreedingGroup <- orig_grps$BreedingGroup %>%
  vapply(., function(x) paste("Breeding Group", x), FUN.VALUE = character(1))

# make sure all the names match up correctly
ID_repair <- read.table("./data/NameUpdates.txt",
                         sep = "\t",
                         header = TRUE)
dapc_grps$ID <- ID_repair[match(dapc_grps$ID, ID_repair$OLD), "NEW"]
orig_grps$ID <- ID_repair[match(orig_grps$ID, ID_repair$OLD), "NEW"]

grps <- left_join(fs_grps, orig_grps) %>%
  left_join(., dapc_grps)

write.table(grps,
            file = "./results/allgrps.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

test <- table(grps[,2:4])
View(test)

upper_left <- table(grps[,c("BreedingGroup", "DAPCGroup")])
lower_right <- table(grps[,c("fsGroup", "BreedingGroup")])

upper_left <- pivot_wider(as.data.frame(upper_left),
                          names_from = DAPCGroup,
                          values_from = Freq)

lower_right <- pivot_wider(as.data.frame(lower_right),
                           names_from = BreedingGroup,
                           values_from = Freq)

mat2 <- as.matrix(upper_left[2:ncol(upper_left)])
rownames(mat2) <- as.vector(upper_left$BreedingGroup)
upper_left <- mat2

mat2 <- as.matrix(lower_right[2:ncol(lower_right)])
rownames(mat2) <- as.vector(lower_right$fsGroup)
lower_right <- mat2

# diagonally bind the matrices

out_mat <- magic::adiag(upper_left, lower_right)

write.table(out_mat, file = "./results/sankey.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA) # use col.name NA to put na extra space for row names


## now need to actually make the diagram

library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)

library("RColorBrewer")
library(networkD3)


# load data
# data <- read.table("~/1workingdir/r/data/sansky.txt", header = TRUE)
data <- read.table("./results/sankey.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
colnames(data) <- gsub(" ", "_", colnames(data))
rownames(data) <- gsub(" ", "_", rownames(data))

# I need a long format
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>% 
  filter(value > 0)

colnames(data_long) <- c("source", "target", "value")

data_long <- data_long[order(data_long$source), ]

# data_long$target <- paste(data_long$target, " ", sep="")
data_long$target <- paste(data_long$target, "", sep="")

# From these flows we need to create a node data frame: it lists every entities involved in the flow

nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique() %>% sort())

nodes

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long$IDsource=match(data_long$source, nodes$name)-1 

data_long$IDtarget=match(data_long$target, nodes$name)-1

# prepare colour scale
ColourScal = 'd3.scaleOrdinal()
.domain(["Pop_1", "Pop_2", "Pop_3", "Pop_4", "Pop_5", "Pop_6", "Breeding_Group_1", "Breeding_Group_2", "Breeding_Group_3", "Breeding_Group_4", "Breeding_Group_5", "Breeding_Group_6", "Breeding_Group_7", "Breeding_Group_8", "DAPC_Group_1", "DAPC_Group_2", "DAPC_Group_3" ])
.range(["#D53E4F", "#fc8d59", "#fee08b", "#e6f598", "#99d594", "#3288bd", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#FC8D59", "#FFFFBF", "#99D594"])'

# Make the Network
tiff(filename = "sankey.tiff",
     width = 1200, height = 600, units = "px")

sankeyNetwork(Links = data_long, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE,
              nodeWidth=40, 
              colourScale = ColourScal,
              fontSize=13, nodePadding=20,
              iterations = 0)

# iterations = 0 does alphabetical order