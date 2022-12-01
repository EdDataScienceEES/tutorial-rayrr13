##%##########################################################%##
#                                                              #
####                     NMDS TUTORIAL                      ####
####                   By Ray Rubia Rankin                  ####
####                  Created on 29/11/2022                 ####
#                                                              #
##%##########################################################%##

# set working directory----
setwd("your_filepath")

# install packages (if you haven't already)-----
install.packages("vegan")
install.packages("permute")
install.packages("lattice")
install.packages("tidyverse")

# load libraries----
library(vegan) # for diversity analysis and ordination methods
library(permute) # necessary to run vegan package
library(lattice) # necessary to run vegan package
library(tidyverse) # for efficient data manipulation and visualization

# import data----
inverts <- read.csv("data/invert_data.csv")

# data manipulation----
# view dataset structure
head(inverts)
str(inverts)

# make community matrix by extracting abundance values from inverts dataframe
invert_community <- inverts[,6:20]

# turn data frame into matrix so it functions in the vegan package
invert_matrix <- as.matrix(invert_community)

# ordination----
# now we will proceed to create the ordination
inv.NMDS <- metaMDS(invert_matrix, distance = "bray", k = 3, autotransform = TRUE, trymax=100) 

# Bray-Curtis distance is chosen because it is not affected by zero values. 
# k represents the number of dimensions we want and is used to reduce stress.
# autotransform arguments ensures the data does a sqrt transformation
# trymax is the number of iterations the algorithm will run

# check stress value
inv.NMDS$stress # ideally stress value is < 0.2

# data visualization----
# using Base R
plot(inv.NMDS) # circles show different communities/sites, crosses show different species
ordiplot(inv.NMDS, type = "n") # create blank ordination plot
orditorp(inv.NMDS, display = "species", col="red", air = 0.1) # add species
orditorp(inv.NMDS, display = "sites", cex = 1.25, air = 0.1) # add sites

# by distance
# ellipse plot
ordiplot(inv.NMDS) # plot shows communities (circles) and species (crosses)
ordiellipse(inv.NMDS, inverts$Distance, label = FALSE, 
            col=c("darkorchid1", "darkslategray1", "bisque1", "brown1"), 
            draw = "polygon", alpha=120) # adding ellipses to the plot
legend("topright", title="Distance (m)",
       c("1","3","7","15"), fill=c("darkorchid1", "darkslategray1", "bisque1", 
                                   "brown1"), horiz=FALSE, cex=.9) # adding a legend

# save plot
png("figures/base_NMDSplot.png", width=6, height=5, units = "in", res = 300)
ordiplot(inv.NMDS)
ordiellipse(inv.NMDS, inverts$Distance, label = FALSE, 
            col=c("darkorchid1", "darkslategray1","bisque1", "brown1"), 
            draw = "polygon", alpha=120)
legend("topright", title="Distance (m)",
       c("1","3","7","15"), fill=c("darkorchid1", "darkslategray1",
                                   "bisque1", "brown1"), horiz=FALSE, cex=.9)
dev.off()

# polygon plot
ordiplot(inv.NMDS) #plot shows communities (circles) and species (crosses)
ordihull(inv.NMDS, groups = inverts$Distance, draw="polygon", col="grey90", label = TRUE) # adding polygons

# spider plot
ordiplot(inv.NMDS) #plot shows communities (circles) and species (crosses)
ordispider(inv.NMDS, groups = inverts$Distance, label = TRUE) # adding spider plot


# using ggplot to make an ellipse plot
# make new dataframe with by extracting NMDS scores
nmds.scores <- as.data.frame(scores(inv.NMDS)$sites) # for newest version of vegan package

# if you have a version of the vegan package <2.6-2 then the following code should also work:
# nmds.scores = as.data.frame(scores(nmds))

# add data from your original dataframe into new NMDS dataframe
# useful if you want to group based on different criteria
nmds.scores <- nmds.scores %>% 
  mutate(Site = as.factor(inverts$Site), Path.Type = as.factor(inverts$Path.Type), 
         Transect = as.factor(inverts$Transect), Distance = as.factor(inverts$Distance))

# check dataframe
head(nmds.scores) 
str(nmds.scores)

# define hidden vegan function that finds coordinates for drawing a covariance ellipse
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# create empty dataframe to combine NMDS data with ellipse data
ellipse_df <- data.frame()

# adding data for ellipse, in this case using distance as a grouping factor
for(g in levels(nmds.scores$Distance)){
  ellipse_df <- rbind(ellipse_df, cbind(as.data.frame(with(nmds.scores[nmds.scores$Distance==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                            wt=rep(1/length(NMDS1),length(NMDS1)))$cov,
                                                                     center=c(mean(NMDS1),mean(NMDS2)))))
                                  ,Distance=g))
}

# create ggplot
(inv_NMDS_plot <- ggplot(data = nmds.scores, aes(NMDS1, NMDS2)) + 
    geom_point(aes(color = Distance, shape = Distance)) + # adding different colours and shapes for points at different distances
    geom_path(data=ellipse_df, aes(x=NMDS1, y=NMDS2, colour=Distance), linewidth=1) + # adding covariance ellipses according to distance # use size argument if ggplot2 < v. 3.4.0
    guides(color = guide_legend(override.aes = list(linetype=c(NA,NA,NA,NA)))) + # removes lines from legend
    theme_bw() + # adding theme
    theme(panel.grid = element_blank()) + # remove background grid
    scale_color_manual(name = "Distance (m)", # legend title
                       labels = c("1","3","7","15"), # adjusting legend labels
                       values = c("gold", "chartreuse3", 
                                  "deepskyblue2", "salmon")) + # customising colours
    scale_shape_manual("Distance (m)", # legend title
                       labels = c("1","3","7","15"), # adjusting legend labels
                       values = c(17, 15, 3, 7))) # customising shapes
# save plot
ggsave(filename = "figures/inv_NMDS_plot.png", inv_NMDS_plot, device = "png")

# data analysis----
# using a PERMANOVA to test the differences in community composition
# This is a PERmutational Multivariate ANalysis Of VAriance and tests the differences between groups, like an ANOVA, but with lots of variables.
# it is essentially a multivariate analysis of variance used to compare groups of objects
inv_permanova <- adonis2(as.matrix(inverts [,6:20]) ~ Distance, inverts, 
                             permutations = 999, method = "bray")
# permutations is the number of possible arrangements the data could take
# using permutations = 999 is standard practice in science and across the literature
# can increase it if you want to get a ore thorough analysis
# use method = bray as it is what we previously used to calculate pairwise distances

# look at model output
inv_permanova

# simper analysis
# this is a SIMilarity PERcentage analysis and compares community differences 
# and reports what species (orders in this scenario)are driving those differences.

inv_simper <- with(inverts[,c(1,2,3,4,5,21)], simper(as.matrix(inverts[,5:20]), Distance, permutations = 100))
# with argument indicates what grouping variables to use to conduct the analysis 
# the number of permutations required to carry out the analysis

# see most influential species (orders in this case) contributing to community differences
inv_simper

# a more detailed summary showing the contributions of all species (or orders)
summary(inv_simper)

# shows species driving differences in community composition between different distances





