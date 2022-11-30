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

# variables of interest must be factors
inverts$Distance = as.factor(inverts$Distance)

str(inverts)

# ordination----
# now we will proceed to create the ordination
inv.NMDS <- metaMDS(inverts [,6:20], distance = "bray", k = 3, trymax=100) 

## Bray-Curtis distance is chosen because it is not affected by zero values. 
##k represents the number of dimensions we want and is used to reduce stress.

# check stress value
inv.NMDS$stress 

## stress values are: <0.05 = excellent, ,0.1 = great, <0.2 = okay, >0.2 = poor. 
## It is standard practice to report the stress value.\

# data visualization----
# using Base R
plot(inv.NMDS) # circles show different communities/sites, crosses show different species
ordiplot(inv.NMDS, type = "n") # create blank ordination plot
orditorp(inv.NMDS, display = "species", col="red", air = 0.1)
orditorp(inv.NMDS, display = "sites", cex = 1.25, air = 0.1)

# by distance
# ellipse plot
ordiplot(inv.NMDS) # plot shows communities (circles) and species (crosses)
ordiellipse(inv.NMDS, inverts$Distance, label = FALSE, col=c("darkorchid1", "darkslategray1",
                                                               "bisque1", "brown1"), draw = "polygon", alpha=120)
legend(1.7, 1.1, title="Distance (m)",
       c("1","3","7","15"), fill=c("darkorchid1", "darkslategray1",
                                        "bisque1", "brown1"), horiz=FALSE, cex=.9)

# save plot
png("figures/base_NMDSplot.png", width=6, height=5, units = "in", res = 300)
ordiplot(inv.NMDS) # plot shows communities (circles) and species (crosses)
ordiellipse(inv.NMDS, inverts$Distance, label = FALSE, col=c("darkorchid1", "darkslategray1",
                                                             "bisque1", "brown1"), draw = "polygon", alpha=120)
legend(1.45, 1.15, title="Distance (m)",
       c("1","3","7","15"), fill=c("darkorchid1", "darkslategray1",
                                   "bisque1", "brown1"), horiz=FALSE, cex=.9)
dev.off()

# polygon plot
ordiplot(inv.NMDS)
ordihull(inv.NMDS, groups = inverts$Distance, draw="polygon", col="grey90", label = TRUE)

# spider plot
ordiplot(inv.NMDS) #plot shows communities (circles) and species (crosses)
ordispider(inv.NMDS, groups = inverts$Distance, label = TRUE)

# using ggplot
# ellipse plot
# create dataframe with the NMDS scores
NMDS1 = data.frame(MDS1 = inv.NMDS$points[,1], MDS2 = inv.NMDS$points[,2], group = inverts$Distance)

str(NMDS1) # checking our group variable is a factor

# create function to create ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell1 <- data.frame() # create empty dataframe

# fill this dataframe with values for creating ellipses
for(g in levels(NMDS1$group)){
  df_ell1 <- rbind(df_ell1, cbind(as.data.frame(with(NMDS1[NMDS1$group==g,],
                                                     veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                  ,group=g))
}

# plot
(inv_NMDSplot<- ggplot(data = NMDS1, aes(MDS1, MDS2)) + 
  geom_point(aes(color = group, shape = group)) + # adding different colours and shapes for points at different distances
  geom_path(data=df_ell1, aes(x=MDS1, y=MDS2, colour=group), linewidth=1) + # use size argument if ggplot2 < v. 3.4.0
  guides(color = guide_legend(override.aes = list(linetype=c(NA,NA,NA,NA)))) + # removes lines from legend
  theme_bw() + 
  theme(panel.grid = element_blank()) + # remove background grid
  scale_color_manual(name = "Distance (m)", # legend title
                     labels = c("1","3","7","15"), # adjusting legend labels
                     values = c("gold", "chartreuse3", 
                                "deepskyblue2", "salmon")) + # customising colours
  scale_shape_manual("Distance (m)", # legend title
                     labels = c("1","3","7","15"), # adjusting legend labels
                     values = c(17, 15, 3, 7))) # customising shapes
# save plot
ggsave(filename = "figures/invNMDS_plot.png", inv_NMDSplot, device = "png")

# data analysis----
# using a PERMANOVA to test the differences in community composition
# This is a PERmutational Multivariate ANalysis Of VAriance and tests the differences between groups, like an ANOVA, but with lots of variables.
# it is essentially a multivariate analysis of variance used to compare groups of objects
inv_permanova <- adonis2(inverts [,6:20] ~ Distance, inverts, 
                             permutations = 999, method = "bray")
inv_permanova

# simper analysis
# this is a SIMilarity PERcentage analysis and compares community differences and reports what species are driving those differences.

(inv.dist.sim <- with(inverts[,c(1,2,3,4,5,21)], simper(inverts[,5:20], Distance, permutations = 100)))
summary(inv.dist.sim)

# shows species driving differences in community composition between different distances





