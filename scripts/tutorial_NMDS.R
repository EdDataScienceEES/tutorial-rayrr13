##%##########################################################%##
#                                                              #
####                     NMDS TUTORIAL                      ####
####                   By Ray Rubia Rankin                  ####
####                  Created on 29/11/2022                 ####
#                                                              #
##%##########################################################%##

# set working directory----
setwd("your_filepath")

# install packages-----
install.packages("vegan")
install.packages("permute")
install.packages("lattice")
install.packages("tidyverse")

# load libraries----
library(vegan) # for diversity analysis and ordination methods
library(permute) # may be necessary for vegan package
library(lattice) # may be necessary for vegan package
library(tidyverse) # for efficient data manipulation and visualization

# import data----
inverts <- read.csv("data/invert_data.csv")



