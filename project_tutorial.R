# Tutorial: Creating a heat map that summarizes two different object's clustering. 
## We will create two data frames, each with list of cell ids and clusters from one object; then combine the two data frames into one; and use that as input to create a heat map that will summarize all the cluster/cell relationships between each object.  

# load these packages 
library(dplyr)
library(Seurat)

# First Goal: Make data frame for first object (ens_metadata) that just lists all cell IDs and cluster information. 

# load the first object (contains only the metadata - has same 10X input as the second object but used the Ens101 transcriptome for alignment)
ens_metadata <- readRDS('ens_metadata.rds')

# if working with an object that has all the data and want to copy just the metadata into a new object 
object_metadata <- your_object@meta.data

# add another column to ens_metadata for cell IDs - we'll just copy that first column
ens_metadata$cell_id <- row.names(ens_metadata)

# get a list of the seurat cluster for each cell: 
ens_clusters <- ens_metadata[['seurat_clusters']]

# and add back to ens_metadata as a new column 
ens_metadata['ens_clusters'] <- ens_clusters
head(ens_metadata)

# now ens_metadata contains the two important columns we need: cell_id and ens_clusters
# let's subset to just have those columns we need in the object 
ens_metadata <- subset(ens_metadata, select = c('cell_id', 'ens_clusters'))
head(ens_metadata)

# Second Goal: Make data frame for second object (NewV432) that just lists all cell IDs and cluster information (exact same procedure as before).

# load the second object (contains only the metadata - has same 10X input as the first object but used the LawsonV432 transcriptome for alignment)
lawson_metadata <- readRDS('lawson_metadata.rds')

# let's add another column to lawson_metadata for cell IDs - we'll just copy that first column
lawson_metadata$cell_id <- row.names(lawson_metadata)

# get a list of the seurat cluster for each cell in the metadata: 
lawson_clusters <- lawson_metadata[['seurat_clusters']]

# and add back to lawson_metadata as a new column 
lawson_metadata['lawson_clusters'] <- lawson_clusters
head(lawson_metadata)

# now lawson_metadata contains the two important columns we need: cell_id and lawson_clusters
# let's subset to just have those columns we need in the object 
lawson_metadata <- subset(lawson_metadata, select = c('cell_id', 'lawson_clusters'))
head(lawson_metadata)

# Third Goal: let's combine both of the data frames we just created into one matrix. 

# this command will combine them, the cell ids pf both objects are in one column and each object's cluster relationship with each cell is listed as two other columns (ens_clusters and lawson_clusters)
combined <- merge(ens_metadata, lawson_metadata, by = 'cell_id', all.x=TRUE, all.y=TRUE)
combined
# good - cell 5 exists in ens but not lawson, cell 9 doesn't exist in ens but does in lawson

# Final Goal: let's use the final matrix we just created as input for the heat map. 

# load package required to create the heat map
library(ggplot2)

# group by the object and the two data columns (lawson_clusters and ens_clusters), ggplot plots the data (n is the number of cells)
p1 <- group_by(combined,lawson_clusters,ens_clusters) %>% summarize(n=n()) %>% 
  ggplot(aes(lawson_clusters, ens_clusters,fill=n)) + geom_tile(aes(fill = n)) 

# geom_text overlays the raw data on top of the grid - the number of cells shared will be overlaid 
p2 <- p1 + geom_text(aes(label = n)) 

# scale_fill_distiller specifies the color pallette desired for heat map - let's use 'Spectral' as it gives greater range of colors. 
p3 <- p2 + scale_fill_distiller(palette = "Spectral")

# let's title the legend - labs is a function in ggplot2 used to modify axes and legends 
p4 <- p3 + labs(fill = "number of shared cells")
p4

# alternate view - entire heat map command together 
group_by(combined,lawson_clusters,ens_clusters) %>% summarize(n=n()) %>% 
  ggplot(aes(lawson_clusters, ens_clusters,fill=n)) + geom_tile(aes(fill = n)) +
  geom_text(aes(label = n)) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "number of shared cells")

## now you have a heat map that displays the relationships of the clusters between the two objects, based on the number of cells shared. 
# good - in ens_clusters, cluster 0 has 14 cells that belong only in the ens object and doesn't exist at all in the lawson object. everything is summarized. 
