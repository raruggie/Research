rm(list=ls(all=T)) # clear global env.
gc()

# install packages:

# install.packages("tidyverse")
# install.packages("ggpmisc")
# install.packages("gridExtra")

# load tidyverse (ggplot2 is in the tidyverse):

library(tidyverse)

#### ggplot introduction ####

# load data:

# dataset 1) df.model:

load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/df.model.Rdata")

# datset 2) df.CQ:

load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/df.CQ.Rdata")

# dataset 1) 

# simple ggplot:

ggplot(data = df.model) +
  geom_point(mapping = aes(x = CSA.percent,
                           y = medC,
                           color = Landuse,
                           size = DA))

# dataset 2) 

# set up df.plot:

df.plot <- df.CQ %>% filter(site == "01349950")

# how do you make these plots:

ggplot(df.plot) +
  geom_point(aes(x = Q, 
                 y = C),
             color = 'blue')

ggplot(df.plot) +
  geom_point(aes(x = Q, 
                 y = C,
                 color = Q_type))

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

# how would you make this plot:

ggplot(df.plot) +
  geom_point(aes(x = Q, y = C)) +
  geom_smooth(aes(x = Q, y = C))

# replace scatter plot with boxplot:

ggplot(df.plot, aes(x = Season, y = C)) +
  geom_point()

ggplot(df.plot, aes(x = Season, y = C)) +
  geom_boxplot()

# density plot:

ggplot(df.plot) +
  geom_density(aes(x = Q, 
                   fill = Q_type), 
               alpha = 0.3)

# predict what this will do:

ggplot(data = df.plot) +
  geom_point(mapping = aes(x = Q,
                           y = C,
                           color = Q_type)) +
  geom_smooth(mapping = aes(x = Q,
                            y = C,
                            color = Q_type))

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

# scales:

ggplot(df.plot, aes(x = Q, y = C)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10()

# facet plot:

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(vars(Q_type, Season), nrow = 2)

# coords:

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  coord_polar()

# labels:

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'Seasonal log-log CQ plot for 01349950',
       y = 'Discharge (cfs)',
       x = 'TP Concentration (mg/L)',
       color = 'Season!')

# theme:

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  theme_dark()

ggplot(df.plot, aes(x = Q, y = C, color = Season)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal()

#### ~ advanced ggplots ~ ####

#### plot 1 - univariate plot ####

# here we are working with the dataset df.CQ
# lets look at the univariate plot of the concentration-discharge 
# relationship for one site 
# we can use the plotting df from the introduction since it was 
# already filtered down to a single site

# we use this df to make the plot and save the plot as a variable called p1:

p1 <- ggplot(df.plot, aes(x = Q, y = C))+
  geom_point()

p1

# the reason we save the plot as a variable is so that we can 
# add to the plot without having to call this setup code again

# lets change the axis labels, add a title, and a trend line:

p1 <- p1 + 
  geom_smooth(method = 'lm') +                                  # add trend line
  xlab('Discharge (cfs)') +                                     # change x-axis label
  ylab('TP Concentration (mg/L)') +                             # change y-axis label
  ggtitle('Concentration-Discharge relationship for 01422747')  # add title

p1

# note this is another way to add/modify these labels. we could have also
# used labs() as shown in the introduction

# here, geom_smooth with the argument method = 'lm' lets us add an OLS trendline

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
# lets look at facet plots of the relationships between the response 
# variables and the watershed percent CSA

# the response variables are in columns 2:4 and 
# our predictor variable CSA.percent is in column 9

# set up a plotting df that contains just these columns:

df.plot <- df.model %>% select(c(2:5,9))

# note we are writing ove df.plot. we will do this for all plots
# in this section

# ggplot likes data in long format anytime you want to make 
# facets, colors, fills, shapes,etc. by a certain identifier

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
  facet_wrap(~name, scales = 'free')

p2

# lets add the trendline equation and adj.R2 to each facet:

p2 <- p2 + stat_poly_eq(use_label(c("eq", "R2")))

p2

# these facets are arranged by alphabetical order
# lets arrange them by adjR2 so we can look at which response variables
# have the best relationships with the predictor variable

# the reason I am having you do all of these steps is to show you what it may
# entail to come up with a solution to your ggplot question

# 1) split df.plot into list of dfs for each response variable:

l.split <- split(df.plot, f = df.plot$name)

# 2) calculate lm for each response variable:

l.split <- lapply(l.split, \(i) lm(value ~ CSA.percent, data = i))

# 3) extract adjR2 from these models:

l.split <- lapply(l.split, \(i) summary(i)$adj.r.squared)

# 4) create a df from this list:

df.adjR2 <- bind_rows(l.split) %>% # combine list elements into df 
  pivot_longer(cols = everything()) %>% # since class(l.split[[1]]) = "numeric", bind_rows returns a wide df with each list element being a column. if you input a list of dfs into bind_rows you get a long df
  rename(adjR2 = value) # rename 'value' column to 'adjR2' such that it is different from value in df.plot

# 5) convert df.plot$name to factor and set its levels based on df.adjR2 (https://stackoverflow.com/questions/40771675/set-levels-of-a-factor-based-on-numeric-value-of-another-column):
# note that we still use df.adjR2 for the levels because in df.plot they are duplicated and that wont work

df.plot$name <- factor(df.plot$name, levels=df.adjR2$name[order(df.adjR2$adjR2, decreasing = TRUE)], ordered=TRUE)

# check factor and its levels:

df.plot$name # looks good

# now we can remake the plot:

p2 <- ggplot(df.plot, aes(x = CSA.percent, y = value)) +
  geom_point() + 
  geom_smooth(method = 'lm') + 
  facet_wrap(~name, scales = 'free') +
  stat_poly_eq(use_label(c("eq", "R2")))

p2

# now the facets arranged by adjR2 value

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

# it is already "Date" which is good, but lets pretend it wasnt already:

df.plot$Date <- as.character(df.plot$Date)

ggplot(df.plot, aes(x=Date, y=value))+
  geom_point()

# the x-axis is messy

# this is why Dates need to be class Date:

df.plot$Date <- as.Date(df.plot$Date)

# make plot:

p3 <- ggplot(df.plot,aes(x=Date, y=value, color = name)) +
  geom_point()

p3

# the dates look good, but, the y axis is dominated by the range of the discharges

# we can use a second y-axis for the values of discharge

# we will need to scale the values of discharge down to the range of values of the concentrations
# and we will need to scale the axis of concentration up to the values of discharges

# I googled 'scale values to new maximum r' and found: https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range

# this is not an r solution but it gives the algorithm to do this

# 1) setting the values needed for the scaling approach:

t.min <- min((df.plot <- df.CQ%>%filter(site == "01349950"))$C)
t.max <- max((df.plot <- df.CQ%>%filter(site == "01349950"))$C)
r.min <- min((df.plot <- df.CQ%>%filter(site == "01349950"))$Q)
r.max <- max((df.plot <- df.CQ%>%filter(site == "01349950"))$Q)

# 2) create function to scale down Q:

fun.scale.down <- function(Q){
  Q.scaled <- ((Q - r.min) / (r.max - r.min)) * (t.max - t.min) + t.min
  return(Q.scaled)
}

# 3) create function to scale up C (for secondary y-axis):

fun.scale.up <- function(C){
  C.scaled <- ((C - t.min) / (t.max - t.min)) * (r.max - r.min) + r.min
  return(C.scaled)
}

# 4) scale Q down to between t.min and t.max and re-set up plotting df:

df.plot <- df.CQ %>% 
  filter(site == "01349950") %>% 
  mutate(across(Q, fun.scale.down)) %>% 
  pivot_longer(cols = c(Q,C))

# now we are ready to remake the plot
# the sec.axis argument in scale_y_continous is used with the 
# fun.scale.up function to display labels that put Q 
# back on its original scale:

p3 <- ggplot(df.plot, aes(x=Date, y=value, color = name))+
  geom_point()+
  scale_y_continuous(
    sec.axis = sec_axis(~ fun.scale.up(.), name = "Discharge (cfs)"),
    "TP Concentration (mg/L)"
  )

p3

# now we have two y-axis

#### plot 4 - combining ggplots ####

# we can combine the plots from this lecture using the grid.arrange
# function from the gridExtra package

# we will use the layout matrix argument to put the timeseries plot
# along the entire width of the combined plot space:

library(gridExtra)

p4 <- grid.arrange(p1, p2, p3, layout_matrix = matrix(c(1, 3, 2, 3), nrow = 2))

p4

# lets add a title to this plot:

p4 <- annotate_figure(p4, top = text_grob("Mashup of plots from this lecture", color = "red", face = "bold", size = 14)) 

p4

