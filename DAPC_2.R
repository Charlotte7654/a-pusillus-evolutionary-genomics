#Clear memory
rm(list = ls())

#Set working directory
setwd("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory.nosync")

#load adegenet for dapc  ####
if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

#make full genlight
full.vcf <-read.vcfR("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/full.vcf")
full.genlight <- vcfR2genlight(full.vcf)

#make neutral genlight
neutral.vcf <-read.vcfR("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/neutral.vcf")
neutral.genlight <- vcfR2genlight(neutral.vcf)

#make outlier genlight
outlier.db.vcf <-read.vcfR("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/outlier_db.vcf")
outlier.db.genlight <- vcfR2genlight(outlier.db.vcf)
outlier.d.vcf <-read.vcfR("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/outlier_d.vcf")
outlier.d.genlight <- vcfR2genlight(outlier.d.vcf)

#number of clusters
num_clust.full <- find.clusters(full.genlight) #99, 2
num_clust.neutral <- find.clusters(neutral.genlight) #99, 2
num_clust.outlier.db <- find.clusters(outlier.db.genlight, max.n.clust=14) #30, 8
num_clust.outlier.d <- find.clusters(outlier.d.genlight, max.n.clust=14) #99, 2

#dapc project 
full.dapc <- dapc(full.genlight, num_clust.full$grp) 
neutral.dapc <- dapc(neutral.genlight, num_clust.neutral$grp)
outlier.db.dapc <- dapc(outlier.db.genlight, num_clust.outlier.db$grp, n.pca = 30, n.da = 8) 
outlier.d.dapc <- dapc(outlier.d.genlight, num_clust.outlier.d$grp)

#define population map
popmap<-read.delim("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/info.nosync/popmap92", header = TRUE, sep = "\t")
View(popmap)
#file layout:
#individual	population
#P01_A01_K4	Kleinzee
#P01_A02_K5	Kleinzee
#...

print(outlier.db.dapc$eig)

#forced second axix
table(num_clust.outlier.db$grp)

#transfrom data
if (!require("tibble", quietly = TRUE))
  install.packages("tibble")
library(tibble)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

dapc_data_df <-
  as_tibble(populations_dapc$ind.coord, rownames = "individual") %>%
  mutate(population = popmap$V2,
         group = populations_dapc$grp)

dapc_data_df.neutral <-
  as_tibble(populations_dapc.neutral$ind.coord, rownames = "individual") %>%
  mutate(population = popmap$V2,
         group = populations_dapc.neutral$grp)

outlier.db.dapc.df <-
  as_tibble(outlier.db.dapc$ind.coord, rownames = "individual") %>%
  mutate(population = popmap$V2,
         group = outlier.db.dapc$grp)
View(outlier.db.dapc.df)

View(outlier.db.dapc)

outlier.db.dapc.df.LD1 <- data.frame(
  LD1 = outlier.db.dapc$ind.coord[, 1],  # Extract the first linear discriminant
  group = outlier.db.dapc$grp,
  population = outlier.db.dapc.df$population # Use group classifications
)

outlier.d.dapc.df <-
  as_tibble(outlier.d.dapc$ind.coord, rownames = "individual") %>%
  mutate(population = popmap$population,
         group = outlier.d.dapc$grp)

# Add population information to DAPC dataframe
outlier.db.dapc.df <- as_tibble(outlier.db.dapc$ind.coord, rownames = "individual") %>%
  mutate(population = popmap$population[match(rownames(outlier.db.dapc$ind.coord), popmap$individual)],
         group = outlier.db.dapc$grp)
View(outlier.db.dapc.df)

# a custom palette for plotting
palette <- c("darkslategray4","darkolivegreen","tomato4", "khaki3","steelblue4", "goldenrod3","darkseagreen4", "sienna3")

palette <- c("Pelican Point" = "darkslategray4", "Mossel Bay" = "darkolivegreen", "Seal Rocks" = "tomato4", "Sinclair Island" = "khaki3", "Kleinzee" = "steelblue4", "Lady Julia Percy" = "goldenrod3","Cape Cross" = "darkseagreen4", "The Skerries" = "sienna3")



palette <- c("steelblue","darkolivegreen","tomato4", "yellow3","darkslategrey", "goldenrod3","darkseagreen3", "sienna3")

palette <- c("Pelican Point" = "darkslategray4", "Mossel Bay" = "darkolivegreen", "Seal Rocks" = "tomato4", "Sinclair Island" = "khaki3", "Kleinzee" = "steelblue4", "Lady Julia Percy" = "goldenrod3","Cape Cross" = "darkseagreen4", "The Skerries" = "sienna3")



palette <- c("darkslategray4","darkolivegreen","tomato4", "khaki3","steelblue4", "goldenrod3","darkseagreen4", "sienna3")

palette <- c("Pelican Point" = "#414487FF", "Mossel Bay" = "#22A884FF", "Seal Rocks" = "#B41500", "Sinclair Island" = "#7AD151FF", "Kleinzee" = "#440154", "Lady Julia Percy" = "#F6AC00","Cape Cross" = "#2A788EFF", "The Skerries" = "#EE6900")

show_col(viridis(6))

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

#outlier
dapc.plot.outlier.db <-
  ggplot(outlier.db.dapc.df, aes(
    x = LD1,
    y = LD2,
    fill = population
  )) +
  ggtitle("Outlier DivBal")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette)+
  theme_bw(base_size = 16)
dapc.plot.outlier.db

ggplot(outlier.db.dapc.df.LD1, aes(x = LD1, fill = group)) +
  geom_density(alpha = 0.5) +  # For density plot
  labs(title = "DAPC 1D Visualization along LD1",
       x = "LD1",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = palette)  # Set your color palette

ggplot(outlier.db.dapc.df.LD1, aes(x = LD1, y = 0, fill = population)) +
  geom_jitter(shape = 21, size = 3, height = 0.05) +  # Spread points slightly to avoid overlap
  labs(title = "DAPC 1D Scatter Plot along LD1",
       x = "LD1",
       y = "") +
  theme_minimal() +
  scale_fill_manual(values = palette)

dapc.plot.outlier.db

unique(outlier.db.dapc.df.LD1$group)
names(palette)

# Check if outlier.db.dapc.df has more than one row for each group.
# If not, ensure your data frame has multiple observations per group.

# Make sure 'group' is a factor and properly assigned in your data.
outlier.db.dapc.df$group <- as.factor(outlier.db.dapc.df$group)

# Check the structure of the data
str(outlier.db.dapc.df)

# Redefine the plot, making sure to include 'fill = group' and possibly 'color = group' within geom_point()
dapc.plot.outlier.db <-
  ggplot(outlier.db.dapc.df, aes(
    x = LD1,
    y = LD2,
    fill = group
  )) +
  ggtitle("Outlier DivBal") +
  geom_point(shape = 21, size = 3, aes(color = group)) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  theme_bw(base_size = 16)

dapc.plot.outlier.db

#full
dapc_plot <-
  ggplot(dapc_data_df, aes(
    x = LD1,
    y = LD2,
    fill = population
  )) +
  ggtitle("Full")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette)+
  theme_bw(base_size = 16)

dapc_plot

dapc_plot_reordered <- dapc_data_df %>%
  mutate(populations = fct_relevel(population, 
                                   "NA", "BU","KA","KO","SW"))

dapc_plot_reordered_plot <- ggplot(dapc_plot_reordered, aes(
  x = LD1,
  y = LD2,
  fill = population
)) +
  ggtitle("Full")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette)+
  theme_bw(base_size = 16)

dapc_plot_reordered_plot


#neutral
dapc_plot.neutral <- ggplot(dapc_data_df.neutral, aes(
  x = LD1,
  y = LD2,
  fill = population
)) +
  ggtitle("Neutral")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette)+
  theme_bw(base_size = 16)

dapc_plot.neutral

#outlier
dapc_plot.outlier <-
  ggplot(dapc_data_df.outlier, aes(
    x = LD1,
    y = LD2,
    fill = population
  )) +
  ggtitle("Outlier")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette)+
  theme_bw(base_size = 16)

dapc_plot.outlier

#Plot axes 2 and 3
#Full 
dapc_plot2 <-
  ggplot(dapc_data_df, aes(
    x = LD2,
    y = LD3,
    fill = population
  )) +
  ggtitle("Full")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette) + 
  theme_bw(base_size = 16)

dapc_plot2

#neutral
dapc_plot.neutral <-
  ggplot(dapc_data_df.neutral, aes(
    x = LD2,
    y = LD3,
    fill = population
  )) +
  ggtitle("Neutral")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette)+
  theme_bw(base_size = 16)

dapc_plot.neutral

#outlier
dapc.plot.outlier.db <-
  ggplot(dapc_data_df.outlier, aes(
    x = LD2,
    y = LD3,
    fill = population
  )) +
  ggtitle("Outlier")+
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = palette)+
  theme_bw(base_size = 16)

#Compoplot 
postprobs <- as.data.frame(round(populations_dapc$posterior, 4))
View(postprobs)

zosteraclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>%
  mutate(population = pop_order$V2,
         group = populations_dapc$grp)

head(postprobs)

#Reorder data
populations.snps.vcf_reordered <- 
  read.vcfR("populations.snps_reordered.vcf") #reoredered by transposing in excel. 
#Error: std::bad_alloc means check the first 15 lines in text file
#include blank line at bottom

populations.snps.genlight_reordered <- vcfR2genlight(populations.snps.vcf_reordered)

samplelist_admix <- read_excel("samplelist_admix.xlsx")
samplelist_admix

num_clust_reordered <- find.clusters(populations.snps.genlight_reordered) #45, 6 

populations_dapc_reordered <- dapc(populations.snps.genlight_reordered, num_clust_reordered$grp) #50,3

#populations_dapc_reordered_matrix<-dapc.matrix(populations.snps.genlight_reordered)

dapc_data_df_reordered <-
  as_tibble(populations_dapc_reordered$ind.coord, rownames = "individual") %>%
  mutate(population = pop_order$pop,
         group = populations_dapc_reordered$grp)

long_palette <- c("#472d7b","#21918c","#fde725","#3b528b", "#440154","#2c728e", "#addc30","#28ae80","#5ec962")
compoplot(populations_dapc_reordered, 
          
          col = c("#2c728e","#3b528b","#addc30", "#440154","#28ae80","#fde725"), 
          xlab = "individuals",
          show.lab = TRUE,
          legend = TRUE,
          posi = "right",
          cex.lab = 1.5,
          cex.fac = 0.3,
          cex.axis = 1,
          cex.main=0.5, 
          cex.sub=0.3, 
          cex.names=0.5,
          txt.leg=NULL)

