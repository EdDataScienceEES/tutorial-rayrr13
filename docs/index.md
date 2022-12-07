---
title: NMDS Plotting and Analysis Using the Vegan Package
subtitle: Analysing differences in ecological community compositions
date: 2022-12-06 20:00:00
author: Ray Rubia
tags: modelling, ecology, species composition
---
# Non-metric Multidimensional Scaling (NMDS)
*By Ray Rubia*

As an ecologist, have you ever wondered how to empirically visualise and test differences in community composition? Are you tired of constantly analysing complex diversity indices and their drivers in community ecology? Are you more interested in the assemblage of species/taxonomical units on which these indices build upon? If you've answered YES to any of the questions above then Christmas has come early for you! I present to you... Non-metric Multidimensional Scaling (NMDS)!!

This tutorial introduces you to the amazing world of NMDS ordinations, from basic visualization to statistical analysis of the results. This tutorial assumes good understanding of the basic principles of R and RStudio, so
if you have just started your R journey I suggest you work through <a href="https://ourcodingclub.github.io/tutorials/intro-to-r/" target="_blank">this tutorial</a> from the Coding Club. We will also be visualising data and building statistical models, so I recommend checking out the following tutorials if you aren't familiar with either of these things:
- Coding Club tutorial on data visualization using `ggplot2` package is available <a href="https://ourcodingclub.github.io/tutorials/datavis/" target="_blank"> here</a>
- Coding Club tutorial on introduction to model design is available <a href="https://ourcodingclub.github.io/tutorials/model-design/" target="_blank"> here</a>.

Advanced R users may also benefit and gain more from this tutorial by checking out <a href="https://ourcodingclub.github.io/tutorials/ordination/" target="_blank"> this tutorial</a> introducing ordinations and NMDS.

This tutorial will introduce the basics of NMDS ordinations. It will also teach you how to conduct and plot a simple NMDS. The second half of this tutorial is slightly more advanced and explores more advanced NMDS plotting as well as statistical modelling of NMDS results.

All files necessary to follow this tutorial can be found in <a href="https://github.com/EdDataScienceEES/tutorial-rayrr13" target="_blank"> this GitHub repository</a>. Download the repository as a ZIP file (clicking `Code/Download ZIP`) and unzip it in your local machine, or clone the repository to your own GitHub account.

### Tutorial Aims

#### <a href="#intro"> 1. Learn what Non-metric Multidimensional Scaling (NMDS) is</a>

#### <a href="#basicNMDS"> 2. Get familiar with conducting a basic NMDS using the `vegan` package</a>

#### <a href="#NMDSviz"> 3. Generate basic NMDS plots </a>

#### <a href="#advancedNMDSviz"> 4. Generate advanced NMDS plots</a>

#### <a href="#stats"> 5. Statistically analyse the results of an NMDS</a>

<a name="intro"></a>

## 1. Learn what Non-metric Multidimensional Scaling (NMDS) is
***

Often in community ecology, research is mainly concerned with studying community diversity and its drivers. However, we are not only interested in how single variables affect/describe communities, we are also intrigued by the community composition itself. The assemblage of species/taxonomic units making up a community is of monumental significance as it determines the functional diversity of a community. Functional diversity encompasses all of the organismal traits that rule ecosystem functioning, dynamics, productivity and stability. This in turn will determine the ecosystem services that us as humans can yield from different habitats. Thus, it is important to understand community assemblages and identify their influencing factors. Nevertheless, these differences in community assemblages are hard and tedious to analyse so they tend to be overlooked.

Non-metric Multidimensional Scaling (NMDS) is a great tool for ecologists to answer this sort of questions. It can condense a lot of information about community species composition into something visible and understandable. For instance, consider a data frame/matrix with relative abundances of different species/other taxonomical units in different communities, where each species/taxonomical unit abundance in the community is an axis and each axis is a dimension (see Figure 1). What NMDS does is summarise all of that information into a 2-dimensional representation, showing the differences in community composition. It tries to represent the original position of a community within a multidimensional space as accurately as possible while also minimising the number of dimensions to easily plot and visualise.

<center><img title = "Multidimensional data example" img src="report_figures/multidimensional_data.png" alt="Img"></center>
*Figure 1. Example of multidimensional data seen through the lens of community ecology, where each species abundance is an axis and each axis is a dimension. Source: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/*

The way NMDS does this is through an indirect gradient analysis ordination approach, meaning it ordinates variables based on scores it generates from a distance matrix. The distance matrix is generated from a community by species (or other taxonomical units) matrix and distances were calculated by comparing pairwise distances (measure of difference in community composition) between 2 different communities in a low dimensional space. This distance matrix represents the compositional dissimilarity between 2 different communities/sites using species counts. But remember, NMDS is a rank-based approach, and once the distances are calculated they are substituted by ranks. The closer two communities are ranked, the more similar composition they have.

If you are still not convinced enough to use NMDS to analyse differences in community composition, here is a list of further reasons why NMDS is such an amazing approach:
- It is non-metric
	- NMDS does not require a normal distribution nor a linear relationship
	- Doesn't require a very specific distribution
	- This is great in ecology as most distributions tend to be skewed (few abundant species and many with low abundance).
- It can tolerate missing pairwise distances
- It can use any dissimilarity measure
	- Doesn't rely on Euclidean distances (classic absolute distance between two points), meaning it can accomodate to more kinds of data
	- Not sensitive to data transformation as absolute distance doesn't matter
- Can group communities based on quantitative, qualitative or mixed variables
- It is a commonly used technique in community ecology and has a well-precedented use
- It can be used for any taxonomic group, not just species, as long as the data fits your research question

Hopefully, you understand some of this hard theory and have a better understanding of what an NMDS is and how it works. To help you grasp and comprehend these concepts even better, we are going to illustrate them by running an example NMDS.

Warnings about NMDS (maybe save for later section)

<a name="basicNMDS"></a>

## 2. Get familiar with conducting a basic NMDS using the `vegan` package
***

Let's get NMDSing!!

You can either run this tutorial by running the pre-written script for this tutorial (see `scripts/tutorial_NMDS.R` in the tutorial repository) or by creating an empty script in R by clicking on `File/New File/R Script` where you can add your own notes and annotations to the code.

Start by setting your working directory to the file path where you have all of the necessary files for this tutorial (probably the folder you just downloaded) and load the packages necessary for this tutorial. This tutorial will be mainly based around the `vegan` package in R which is a great tool for descriptive community ecology analysis e.g. ordinations, diversity analysis and community dissimilarity analysis. This package was developed in GitHub by Finnish community ecologist Jari Oksanen and his development team.
~~~r
# set working directory----
setwd("your_filepath")

# load libraries----
library(vegan) # for diversity analysis and ordination methods
library(permute) # necessary to run vegan package
library(lattice) # necessary to run vegan package
library(tidyverse) # for efficient data manipulation and visualization, contains ggplot2, dplyr and tidyr packages
~~~
*Note: If you haven't installed these packages on your R environment yet you can run the code `install.packages("package_name")` to install them. The packages `permute` and `lattice` are necessary to ensure smooth running of all of the functions on the `vegan` package*

Some brief background for our research question today: The IDH essentially states that diversity is maximised at intermediate levels of disturbance. At high levels of disturbance, good dispersers dominate the communities, not giving strong competitors a chance to settle. Conversely, at low rates of disturbance, stronger competitors dominate and competitively exclude weaker competitors (tend to be the good dispersers) and dominate the community. Thus, at different levels of disturbance, community assemblages should be different as at high levels we have rapid colonising species and at low levels we have competitively dominant species.

Now we can import our dataset. This data set was collected near Oban, Scotland during an undergraduate ecology field course by University of Edinburgh undergraduate students. It was collected to test the accuracy of the classic intermediate disturbance hypothesis (IDH) in predicting how the frequency and intensity of disturbances affects diversity. The data can be found in the repository under `data/invert_data.csv` and can be imported as follows:
~~~r
# import data----
inverts <- read.csv("data/invert_data.csv")
~~~

The data consists of invertebrate order abundance measured at varying distances from paths (1, 3, 7 and 15 m). Distance from paths was used as a proxy to quantify different levels of disturbance. Data was collected using 1 m<sup>2</sup> quadrats at each distance at 5 different grassland sites, along 3 different transects in each site, around Oban, Scotland, giving a total sample size of *n = 60*. The aim of the study is to see if communities at different distances from paths had different community composition and provided empirical evidence supporting the IDH. This leaves us with the following research question for our study: **Does the order composition of invertebrate communities vary at different distances from paths?**

*Note: In this case we can use order as our taxonomical unit of interest rather than species as it is the functional traits and their diversity (competitors vs dispersers) we are interested in, not the species themselves.*

***
### Conducting a basic NMDS

Before we carry out any sort of ordination, we first need to look at our data and analyse its structure as well as its variables.

~~~r
# view dataset structure
head(inverts)
str(inverts)
~~~

We can see that the dataset contains columns with information about environmental variables such as the site the data was collected at, the path type, the transect number and the distance from the path of each community. The "Orders" column indicates the total amount of orders found in each quadrat. The "Abundance" column indicates the total number of individuals found in each quadrat. However, we are more interested in all of the other columns which have information about invertebrate order abundance in the different communities. Each column is a different order and it contains information about that order's abundance in different communities. This will be the crucial data we use to construct our NMDS.

*Note: this data is in wide format. If your data is in long format and all variables regarding species abundance in a single column, you want to convert this data to wide format. This can be done using the `spread()` function from the `dplyr` package.*

Notice how the data is in the format of a dataset. However, as we said above, in order to conduct an NMDS we said that our dataset must take the shape of a community by species (or orders in this case) matrix. Thus, we have to extract the abundance values from our dataframe and turn them into a community by order matrix to ensure smooth running in the `vegan` package.
~~~r
# make community matrix by extracting abundance values from inverts dataframe
invert_community <- inverts[,6:20]

# turn data frame into matrix so it functions in the vegan package
invert_matrix <- as.matrix(invert_community)
~~~

We can now create our NMDS ordination by using the `metaMDS` function from the `vegan` package. We will save this NMDS as an object in our R environment (called "inv.NMDS").
~~~r
inv.NMDS <- metaMDS(invert_matrix, distance = "bray", k = 3, autotransform = TRUE, trymax=100)
~~~

Pay attention to the different arguments in the above line of code that illustrate the different decisions we have made. First, we call our newly created matrix object containing a community by order matrix.

Secondly, as classic Euclidean distances are sensitive to species/order abundances or absences, we have decided to use Bray-Curtis as our measure of distance to build our distance (or dissimilarity) matrix on which the NMDS runs on, indicated by the `distance = "bray"` argument. The distance measure used has a strong influence over the analysis and the final output, so selecting the right one is very important. Bray-Curtis is a distance measure commonly used in ecology in this sort of analysis for many reasons: (1) it is often used for raw count data (as is the case with our data), (2) it is invariant to changes in units, (3) it accounts for both species presence/absence as well as species abundance, (4) it is unaffected by additions or removals of species not ubiquitous across all communities, (5) it is unaffected by adding a new community, and (6) it recognizes differences in total species/order abundances when relative abundances are equal.

Next, we need to decide the reduced number of dimensions the NMDS algorithm will ordinate objects in. We communicate this to R through the `k = 3` argument. Increasing the number of dimensions will allow for a better representation of the original multidimensional distance matrix and lower NMDS stress values (goodness of fit). However, the more dimensions, the more challenging the interpretation becomes (especially if k > 3). Therefore, 2 or 3 dimensions are usually the consensus in ecology. We have decided to do 3 dimensions in order to provide a better fit to our data.

The `autotransform = TRUE` argument is simply to ensure that R performs a square-root transformation of the raw data before calculating Bray-Curtis distances (it is standard procedure in an NMDS). This reduces the relative influence of the most abundant species, so that they don't dominate the distance matrix.

Finally, it is worth noting that NMDS is an iterative algorithm, meaning it repeats the same process over and over again until it finds the best solution. The number of iterations you want it to run can be set by using the `trymax = ...` argument. In this case we have decided it to run 100 times. As it is an iterative algorithm, it is recommended to run an NMDS multiple times to ensure a stable solution has been reached.

Lets now assess how well our ordination has fitted our data by looking at its stress value:

~~~r
# check stress value
inv.NMDS$stress # ideally stress value is < 0.2
~~~

The stress of an NMDS indicates its goodness of fit. It is essentially the distance between the reduced dimensional space and the complete one. NMDS tries to optimise this as much as possible and it can be reduced by increasing the dimensionality of your ordination. Generally, stress values < 0.2 are considered a good fit, and stress values under 0.1 or 0.05 are ideal. In this case, we have a stres value of `~0.114`, indicating our ordination is a good fit for our data. The stress value of an NMDS should always be reported in the NMDS figure or NMDS figure caption.

And there you have it! You have now built your first NMDS, well done! We are now going to move on to the fun and interesting part which is plotting our NMDS results! But before, a few warnings to consider when working with an NMDS ordination:
#### **<span style="color:red">Warnings:</span>**
- Beware of NMDS interpretation if you have high stress values
- After running an NMDS, the generated distance matrix is no longer raw dat as pairwise distances have been calculated from that data, making the distance data interdependent
- Do not relate environmental variables to an NMDS ordination as it is only an approximate representation of the distance matrix. It is more relevant to relate variables to the distance matrix itself
- NMDS can be computationally demanding for large datasets
- As it is an iterative algorithm, each NMDS ordination or plot you generate from scratch can look slightly different.
- If you have an extremely low stress value your data might have too many zeros. R usually gives you a warning message about this anyways.

<a name="NMDSviz"></a>

## 3. Generate basic NMDS plots
***
Let's get straight into plotting our NMDS! We will first use the `vegan` package and `Base R` to generate some simple NMDS plots. More advanced NMDS plotting using `ggplot2` will be explored in <a href="#advancedNMDSviz"> section 4</a> of this tutorial.

Let's start by plotting the most basic of NMDS plots.
~~~r
# data visualization----
# using vegan package and Base R
plot(inv.NMDS) # circles show different communities/sites, crosses show different species
~~~

This should generate the following plot:

<center><img title = "Basic NMDS plot" img src="report_figures/NMDS_most_basic_plot.png" alt="Img"></center>
*Figure 2. Basic representation of an NMDS, where circles show different communities and crosses show different invertebrate orders.*

As you can see, this plot doesn't give away much information. Each dot represents a different community and each cross represents a different invertebrate order. The closer 2 points are, the more similar those communities are in terms of composition. However, here it is hard to identify the position of different communities or species because there are no labels or no other indicators. Let's add these! We will use the `ordiplot` and `orditorp` plotting functions embedded in the `vegan` package.

~~~r
ordiplot(inv.NMDS, type = "n") # create blank ordination plot
orditorp(inv.NMDS, display = "species", col="red", air = 0.1) # add order names in red
orditorp(inv.NMDS, display = "sites", cex = 1.25, air = 0.1) # add site numbers in black
~~~
<center><img title = "Messy basic NMDS plot" img src="report_figures/base_NMDSplot_messy.png" alt="Img"></center>
*Figure 3. NMDS representation with numbers representing the different sites and labels indicating the different species.*

We can now see what dot corresponds to each community and what cross
corresponds to each invertebrate order. The `air = ...` argument is simply just to add empty space around text labels. Nevertheless, the plot still looks a bit confusing and messy. It isn't telling us anything about the relationship between different environmental variables and community composition. It is simply telling us what communities are more similar to one another without giving any background environmental information for these communities.

Luckily enough, `vegan` has a solution to this! – It allows you to plot different shapes or lines around communities based on different grouping factors. Our research question is mainly concerned with how distance from disturbance (paths) affects invertebrate community composition. Thus, using the functionality of the `vegan` package we can add different shapes or lines grouping together communities (dots on the NMDS plot) at different distances. These plots can be ellipse plots, polygon plots or spider plots, among many others. This will allow us to understand our NMDS plots much better and will give a better visual representation of the effect distance from paths (disturbances) on community composition.

As it is my personal favourite and I think it is the best at depicting differences, I will thoroughly demonstrate how to create an ellipse plot using the `ordiellipse()` function from the `vegan` package and add all of the corresponding aesthetics. I will provide the key code to generate spider and polygon plots. If you prefer these and want to further develop them into a nicer plot, the aesthetics work exactly the same as for the ellipse plot.

#### Ellipse plot

We first create an NMDS ordination plot showing the dots (communities) and crosses (species, or orders in this case) and then overwrite this plot and add the ellipses grouping communities according to their distances from paths using the `ordiellipse()` functions. This will allow us to see if there is any distinction in community composition at different distances from paths. The less overlap there is among the ellipses, the more distinct communities at different distances from disturbances are. The rest of the arguments are for aesthetics: adding different colours to the ellipses, adjusting their transparency, removing title and adding a legend with the corresponding title and labels. You can adjust these aesthetics depending on personal preference.

~~~r
# ellipse plot
ordiplot(inv.NMDS) # plot shows communities (circles) and species (crosses)
ordiellipse(inv.NMDS, inverts$Distance, label = FALSE,
            col=c("darkorchid1", "darkslategray1", "bisque1", "brown1"),
            draw = "polygon", alpha=120) # adding ellipses to the plot
legend("topright", title="Distance (m)",
       c("1","3","7","15"), fill=c("darkorchid1", "darkslategray1", "bisque1",
                                   "brown1"), horiz=FALSE, cex=.9) # adding a legend

# save plot
png("your_filepath/image.png", width=6, height=5, units = "in", res = 300)
ordiplot(inv.NMDS)
ordiellipse(inv.NMDS, inverts$Distance, label = FALSE,
            col=c("darkorchid1", "darkslategray1","bisque1", "brown1"),
            draw = "polygon", alpha=120)
legend("topright", title="Distance (m)",
       c("1","3","7","15"), fill=c("darkorchid1", "darkslategray1",
                                   "bisque1", "brown1"), horiz=FALSE, cex=.9)
dev.off()
~~~
The plot should end up looking something like this:
<center><img title = "NMDS ellipse plot" img src="report_figures/base_NMDSplot.png" alt="Img"></center>
*Figure 4. NMDS ellipse plot showing invertebrate community order composition at different distances from disturbances (paths)*

#### Polygon plot
Polygon plots are created just the same as ellipse plots, the only difference is that now we are using the `ordihull()` function rather than `ordiellipse()` to add polygons surrounding communities at different distances from  paths. If you want to add aesthetics or save this plot the code and arguments work exactly the same as for ellipse plots.

~~~r
# polygon plot
ordiplot(inv.NMDS) #plot shows communities (circles) and species (crosses)
ordihull(inv.NMDS, groups = inverts$Distance, draw="polygon", col="grey90", label = TRUE) # adding polygons
~~~
The plot should end up looking something like this:
<center><img title = "NMDS polygon plot" img src="report_figures/base_NMDS_polygon_plot.png" alt="Img"></center>
*Figure 5. NMDS polygon plot showing invertebrate community order composition at different distances from disturbances (paths)*

#### Spider plot
Spider plots are created just the same as ellipse and polygon plots, the only difference is that now we are using the `ordispider()` function to add lines connecting all of the communities that are found at the same distance from disturbances. This way, we can see what communities are found at different distances and see if there are any distinctions in community assemblages. If you want to add aesthetics or save this plot the code and arguments work exactly the same as for ellipse plots.

~~~r
# spider plot
ordiplot(inv.NMDS) #plot shows communities (circles) and species (crosses)
ordispider(inv.NMDS, groups = inverts$Distance, label = TRUE) # adding spider plot
~~~
The plot should look something like this:
<center><img title = "NMDS spider plot" img src="report_figures/base_NMDS_spider_plot.png" alt="Img"></center>
*Figure 6. NMDS spider plot showing invertebrate community order composition at different distances from disturbances (paths)*

Phew! That was a lot of information to take in! We are now going to move on to more advanced NMDS plotting and some statistical analysis but let's take a breath first! Help yourself to a cup of tea or coffee or any other beverage of your choice. Get some fresh air and come prepared to tackle the second half of this tutorial!

*P.S. Don't be afraid to go over this first half again or revisit certain topics/definitions if you are struggling with grasping some of the concepts*

<a name="advancedNMDSviz"></a>

## 4. Generate advanced NMDS plots
***
More text, code and images.


<a name="stats"></a>

## 5. Statistically analyse the results of an NMDS</a>
***

This is the end of the tutorial. Summarise what the student has learned, possibly even with a list of learning outcomes. In this tutorial we learned:

##### - how to generate fake bivariate data
##### - how to create a scatterplot in ggplot2
##### - some of the different plot methods in ggplot2

We can also provide some useful links, include a contact form and a way to send feedback.

For more on `ggplot2`, read the official <a href="https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf" target="_blank">ggplot2 cheatsheet</a>.

Everything below this is footer material - text and links that appears at the end of all of your tutorials.
<hr>
<hr>

#### Check out our <a href="https://ourcodingclub.github.io/links/" target="_blank">Useful links</a> page where you can find loads of guides and cheatsheets.

#### If you have any questions about completing this tutorial, please contact us on ourcodingclub@gmail.com

#### <a href="INSERT_SURVEY_LINK" target="_blank">We would love to hear your feedback on the tutorial, whether you did it in the classroom or online!</a>

<ul class="social-icons">
	<li>
		<h3>
			<a href="https://twitter.com/our_codingclub" target="_blank">&nbsp;Follow our coding adventures on Twitter! <i class="fa fa-twitter"></i></a>
		</h3>
	</li>
</ul>

### &nbsp;&nbsp;Subscribe to our mailing list:
<div class="container">
	<div class="block">
        <!-- subscribe form start -->
		<div class="form-group">
			<form action="https://getsimpleform.com/messages?form_api_token=de1ba2f2f947822946fb6e835437ec78" method="post">
			<div class="form-group">
				<input type='text' class="form-control" name='Email' placeholder="Email" required/>
			</div>
			<div>
                        	<button class="btn btn-default" type='submit'>Subscribe</button>
                    	</div>
                	</form>
		</div>
	</div>
</div>
