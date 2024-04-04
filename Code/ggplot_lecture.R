rm(list=ls(all=T)) # clear global env.
gc()

# install packages:

# install.packages("tidyverse")
# install.packages("ggpmisc")
# install.packages("gridExtra")
# install.packages("ggpubr")

# load tidyverse (ggplot2 is in the tidyverse):

library(tidyverse)

# load data:

# dataset 1) df.model:

load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/df.model.Rdata")

# datset 2) df.CQ:

load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/df.CQ.Rdata")

#### ggplot introduction ####

# simple ggplot using dataset 1:

ggplot(data = df.model) +
  geom_point(mapping = aes(x = CSA.percent,
                           y = medC,
                           color = Landuse,
                           size = DA))

# take a look at first dataset df.model:

df.model

# practice: play around with different data columns in different aesthetics
# put your code here:










# take a look at dataset 2 - df.CQ

df.CQ 

# set up df.plot:

df.plot <- df.CQ %>% filter(site == "01349950")

# how do you make these plots
# put your code here:











# note on colors:

ggplot(df.plot) +
  geom_point(aes(x = Q, y = C, color = 'blue'))

ggplot(df.plot) +
  geom_point(aes(x = Q, y = C), color = 'blue')

# same aes, different geom:

ggplot(df.plot) +
  geom_point(aes(x = Q, y = C))

ggplot(df.plot) +
  geom_smooth(aes(x = Q, y = C))

# how would you make this plot
# put your code here:











# replace scatter plot with boxplot
# put your code here:










# make the density plot
# put your code here:











# mapping inside ggplot or geom:

# local:

ggplot(data = df.plot) +
  geom_point(mapping = aes(x = C,
                           y = Q,
                           color = Q_type))

# global:

ggplot(data = df.plot, 
       mapping = aes(x = C, y = Q, 
                     color = Q_type)) +
  geom_point()

# Any aesthetics in geom_* layer only apply to that layer:

ggplot(df.plot, mapping = aes(x = Q, y = C)) +
  geom_point(mapping = aes(color = Q_type)) +
  geom_smooth(method = 'lm')

# saving stock ggplot to variable:

p <- ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')

p

# scales - change axis ticks:

p + scale_x_continuous(breaks = seq(0,2500, 250))

# practice using scales - log x and y axis of CQ plot
# put your code here:










# practice using facet - make flow type x season facets
# put your code here:









# coords:

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  coord_polar()

# practice labels - add title, change x and y labels, change legend title for color aes
# put your code here:










# theme:

# dark:

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  theme_dark()

# minimal:

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal()

#### ~ advanced ggplots ~ ####

#### plot 1 - univariate plot ####

# here we are working with the df.CQ dataset 
# lets look at the scatterplot of the concentration-discharge relationship for one site 
# we can use the plotting df from the introduction since it was already filtered down to a single site

# we use this df to make the plot and save the plot as a variable called p1:

p1 <- ggplot(df.plot, aes(x = Q, y = C))+
  geom_point()

p1

# lets and a linear trend line, change the axis labels, and add a title:

p1 <- p1 + 
  geom_smooth(method = 'lm') +  
  labs(x = 'Discharge (cfs)',
       y = 'TP Concentration (mg/L)',
       title = 'Concentration-Discharge relationship for 01422747')

p1

# usually the C-Q relationship is displayed in log space:

p1 <- p1 +
  scale_x_log10() +
  scale_y_log10()

p1

# lets add the trendline equation and adj.R2 to the plot

# we will use a function call stat_poly_eq() from the ggpmisc package:

library(ggpmisc)

p1 <- p1 + stat_poly_eq(use_label(c("eq", "R2")))

p1

# note the equation is calculated for the log-log relationship in log base 10 

#### plot 2 - facet plots ####

# here we are working with the dataset df.model
# lets look at facet plots of the relationships between the response variables and the watershed percent CSA

# the response variables are in columns 2:5 and our predictor variable CSA.percent is in column 9

# set up a plotting df that contains just these columns:

df.plot <- df.model %>% select(c(2:5,9)) # note we are writing over df.plot on purpose

# ggplot likes data in long format because...

# "the entire tidyverse ecosystem [is built] on the concept of tidy data, which is essentially data in long format. 
# The basic reason for working with long-format data is that the same data can be represented by many wide formats, 
# but the long format is typically unique"

# we can use pivot_longer to make df.plot long:

df.plot <- df.plot %>% pivot_longer(cols = -CSA.percent)

# what we have now are three columns:

# CSA.percent, which will be our x variable
# value, which will be our y variable
# name, which will separate the facets

# now we are ready to make our plot:

p2 <- ggplot(df.plot, aes(x = CSA.percent, y = value)) +
  geom_point() + 
  geom_smooth(method = 'lm') + 
  facet_wrap(vars(name), scales = 'free')

p2

# we can add the trendline equation and adj.R2 to each facet:

p2 <- p2 + stat_poly_eq(use_label(c("eq", "R2")))

p2

# these facets are arranged by alphabetical order
# lets arrange them by adjR2 so we can easily look at which response variables have the best relationship with the predictor variable

# to do this we need to change the structure of the column we facet on
# class 'factor' lets us add another layer of information to a character data column
# with factor we can order a character string based on whatever we want

# *the reason I am having you do all of these steps is to show you what it may entail to come up with a solution to your ggplot question*

# here are the steps to set the levels of the column to facet on:
# (this is just one of probably many other simpler solutions, but this is how I though of coding it)

# 1) split df.plot into list of dfs for each response variable:

l.split <- split(df.plot, f = df.plot$name)

# 2) calculate a linear model for each response variable, which are the list elements:

l.split <- lapply(l.split, \(i) lm(value ~ CSA.percent, data = i))

# 3) extract the adjR2 from these models:

l.split <- lapply(l.split, \(i) summary(i)$adj.r.squared)

# 4) merge this list into a single dataframe of adjR2 values for each response variable:

df.adjR2 <- bind_rows(l.split) %>%        # combine list elements into df 
  pivot_longer(cols = everything()) %>%   # since class(l.split[[1]]) = "numeric", bind_rows returns a wide df with each list element being a column. if you input a list of dfs into bind_rows you get a long df
  rename(adjR2 = value)                   # rename 'value' column to 'adjR2' such that it is different from value in df.plot

# 5) now we are ready to convert our faceting variable, df.plot$name, to a factor and set its levels 
# we set the levels for df.plot$name based on the numeric order of df.adjR2$adj.R2 (https://stackoverflow.com/questions/40771675/set-levels-of-a-factor-based-on-numeric-value-of-another-column):
# what this is doing is arranging df.adjR2$name (which are the same names as in df.plot$name) by the descending numerical order of the adj.R2 column

df.plot$name <- factor(df.plot$name, levels=df.adjR2$name[order(df.adjR2$adjR2, decreasing = TRUE)], ordered=TRUE)

# check factor and its levels:

df.plot$name # looks good

# 6) now we can remake the plot

# since we just want to swap the dataset used in the ggplot, we can 'pipe' in our updated df.plot dataframe using a ggplot pipe operator '%+%' 
# this works because none of our column names changed, just the data that was in the columns:

p2 <- p2 %+% df.plot

p2

# now the facets arranged by adjR2 value

# I hope this wasnt too confusing. I use ordered factors all the time to 'trick' ggplot

#### plot 3 - timeseries plot ####

# here we are using df.CQ 
# lets look at the timeseries of flow and concentrations for a single site
# we can use the site we used in the introduction

# setting up a plotting df:

df.plot <- df.CQ %>% 
  filter(site == "01349950") %>% 
  pivot_longer(cols = c(Q,C))

# always check class of Date column before plotting:

class(df.plot$Date) 

# it is already "Date" which is good, but lets pretend it wasn't already

# what would happen if we read in the data and Date happened to be a character?:

df.plot$Date <- as.character(df.plot$Date)

ggplot(df.plot, aes(x=Date, y=value))+
  geom_point()

# the x-axis is messy

# this is why dates need to be class 'Date':

df.plot$Date <- as.Date(df.plot$Date)

# make plot:

p3 <- ggplot(df.plot,aes(x=Date, y=value, color = name)) +
  geom_point()

p3

# the dates look good, but, the y axis is dominated by the range of the discharges

# we will use a second y-axis for the values of discharge

# to do this, we need to scale the values of discharge down to the range of values of the concentrations
# and we will need to scale the axis of concentration up to the values of discharges

# I googled 'scale values to new maximum r' and found: 

# https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range

# this is not an r solution but it gives the algorithm to do this

# we can code it by hand:

# 1) set the values needed for the scaling approach: take the minimum and maximum of the concentrations and discharges:

t.min <- min(df.plot$value[which(df.plot$name=='C')])
t.max <- max(df.plot$value[which(df.plot$name=='C')])
r.min <- min(df.plot$value[which(df.plot$name=='Q')])
r.max <- max(df.plot$value[which(df.plot$name=='Q')])

# 2) create a function to scale down Q:

fun.scale.down <- function(Q){
  Q.scaled <- ((Q - r.min) / (r.max - r.min)) * (t.max - t.min) + t.min
  return(Q.scaled)
}

# 3) create a function to scale up C (for secondary y-axis):

fun.scale.up <- function(C){
  C.scaled <- ((C - t.min) / (t.max - t.min)) * (r.max - r.min) + r.min
  return(C.scaled)
}

# 4) re-set up plotting df with Q scaled down (done inside mutate):

df.plot <- df.CQ %>% 
  filter(site == "01349950") %>% 
  mutate(across(Q, fun.scale.down)) %>% 
  pivot_longer(cols = c(Q,C))

# 5) now we are ready to remake the plot - 
# the sec.axis argument in scale_y_continous is used with the fun.scale.up function to display labels
# that put Q back on its original scale:

p3 <- ggplot(df.plot, aes(x=Date, y=value, color = name))+
  geom_point()+
  scale_y_continuous(
    "TP Concentration (mg/L)",                         # this names the left hand y-axis
    sec.axis = sec_axis(~ fun.scale.up(.), 
                        name = "Discharge (cfs)")      # this names the right hand y-axis
  )

p3

# now we have two y-axis and we can see the full range of values for C and Q

#### plot 4 - combining ggplots ####

# we can combine the plots from this lecture using the grid.arrange function from the gridExtra package

# we will use the layout matrix argument to put the timeseries plot along the entire width of the combined plot space:

library(gridExtra)

p4 <- grid.arrange(p1, p2, p3, 
                   layout_matrix = matrix(c(1, 3, 2, 3), nrow = 2),
                   top = textGrob("Mashup of plots from this lecture",gp=gpar(fontsize=14, font=3, face = 'bold'), color = 'red'))

p4

# we an add a title to this plot using ggpubr::annotate_figure

library(ggpubr)

p4 <- annotate_figure(p4, top = text_grob("", color = "red", face = "bold", size = 14)) 

p4



