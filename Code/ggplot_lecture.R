rm(list=ls(all=T)) # clear global env.
gc()

# install packages:

# install.packages("tidyverse")
# install.packages("ggpmisc")
# install.packages("ggpubr")
# install.packages("scales")

# load data:

load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/df.Seg.Rdata")
load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/df.setup.Rdata")

#### plot 1 - univariate plot ####

# lets look at the univariate plot of the concentration-discharge 
# relationship for one USGS site in df.Seg

# first we set up a plotting df:

library(tidyverse)

df.plot <- df.Seg %>% filter(site == "01422747")

# then we use this df to make the plot:

p1 <- ggplot(df.plot, aes(x = Q_real, y = C))+
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

# here, geom_smooth with the argument method = 'lm' lets us add an OLS trendline to the relationship.

# usually the C-Q relationship is displayed in log space:

p1 <- p1 +
  scale_x_log10() +
  scale_y_log10()

p1

# lets add the trendline equation and adj.R2 to the plot:

library(ggpmisc)

p1 <- p1 + stat_poly_eq(use_label(c("eq", "R2")))

p1

# note the equation is calculated for the log-log relationship in log base 10 

# one last item to add to this plot

# df.Seg contains a column called Seg_C which contains the fitted values of 
# a breakpoint analysis, which determined the best two slope model for the
# CQ relationship

# the fitted values are in log space so we need to do e^(Seg_C)

# lets add the two-slope trend line to p1:

p1 <- p1 + geom_line(aes(x = Q_real, y = exp(Seg_C)),color = 'tomato')

p1

#### plot 2 - facet plots ####

# lets look at facet plots of the relationships between the response 
# variables and the watershed percent agriculture (R_PLANTNLCD06)

# taking a look at the column names of our df:

names(df.setup)

# the response variables are in columns 2:18 and R_PLANTNLCD06 is in column 23

# set up plotting df that contains just these columns:

df.plot <- df.setup %>% select(c(2:18,23))

# ggplot likes data in long format anytime you want to make 
# facets, colors, fills, shapes,etc. by a certain identifier

# we can use pivot_longer to make df.plot long:

df.plot <- df.plot %>% pivot_longer(cols = -R_PLANTNLCD06)

# what we have now are three columns:

# R_PLANTNLCD06, which will be our x variable
# value, which will be our y variable
# name, which will seperate the facets

# now we are ready to make our plot:

p2 <- ggplot(df.plot, aes(x = R_PLANTNLCD06, y = value)) +
  geom_point() + 
  geom_smooth(method = 'lm') + 
  facet_wrap(~name, scales = 'free')

p2

# note use of argument scales = 'free' in facet_wrap

# what happens if you dont use that argument? 

ggplot(df.plot, aes(x = R_PLANTNLCD06, y = value)) +
  geom_point() + 
  geom_smooth(method = 'lm') + 
  facet_wrap(~name)

# the variable with the highest values dictates the axis limits for all facets

# lets add the trendline equation and adj.R2 to each facet:

p2 <- p2 + stat_poly_eq(use_label(c("eq", "R2")))

p2

# these facets are arranged by alphabetical order
# lets arrange them by adjR2

# 1) split df.plot into list of dfs for each response variable:

l.split <- split(df.plot, f = df.plot$name)

# 2) calcualte lm for each response variable:

l.split <- lapply(l.split, \(i) lm(value ~ R_PLANTNLCD06, data = i))

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

p2 <- ggplot(df.plot, aes(x = R_PLANTNLCD06, y = value)) +
  geom_point() + 
  geom_smooth(method = 'lm') + 
  facet_wrap(~name, scales = 'free') +
  stat_poly_eq(use_label(c("eq", "R2")))

p2

# now the facets arranged by adjR2 value

#### plot 3 - distributions of response variables ####

# lets compare the distributions of the different AANY estimation methods

# setting up a plotting df with just these columns:

df.plot <- df.setup %>% select(contains('AANY')) %>%
  pivot_longer(cols = everything())

# violin plots visualize the distribution (i.e boxplot + histogram):

p3 <- ggplot(df.plot, aes(x=name, y=value))+
  geom_violin()

p3

# AANY.simple has a much larger range than the others...
# it looks like it maybe due to outliers:

boxplot(df.setup$AANY.simple)$out

# I'm sure if we remove these observations, the violin plots will look better...
# but justifying removing those is not easy, so lets work with them

# we can use a second y-axis for the values of AANY.simple

# we will need to scale AANY.simple down to the range of values of the other
# three variables, and scale the range of the other three variables up to
# the range of AANY.simple for the second y-axis

# I googled 'scale values to new maximum r' and found: https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range

# this is not an r solution but it gives the algorithm to do this

# 1) setting the values needed for the scaling approach:

t.min <- min((df.plot %>% filter(name !="AANY.simple"))$value)
t.max <- max((df.plot %>% filter(name !="AANY.simple"))$value)
r.min <- min((df.plot %>% filter(name =="AANY.simple"))$value)
r.max <- max((df.plot %>% filter(name =="AANY.simple"))$value)

# 2) create function to scale down AANY.simple:

fun.scale.down <- function(AANY.simple){
  AANY.scaled <- ((AANY.simple - r.min) / (r.max - r.min)) * (t.max - t.min) + t.min
  return(AANY.scaled)
}

# 3) create function to scale up all others (for secondary y-axis):

fun.scale.up <- function(AANY.method2){
  AANY.scaled <- ((AANY.method2 - t.min) / (t.max - t.min)) * (r.max - r.min) + r.min
  return(AANY.scaled)
}

# 4) scale AANY.simple down to between t.min and t.max and re-set up plotting df:

df.plot <- df.setup %>% 
  mutate(across(contains("simple"), fun.scale.down)) %>% 
  select(contains('AANY')) %>%
  pivot_longer(cols = everything())

# now we are ready to remake the plot - the sec.axis argument in scale_y_continous 
# is used with the fun.scale.up function to display labels that put AANY.simple back on its original scale:

p3 <- ggplot(df.plot, aes(x=name, y=value))+
  geom_violin()+
  scale_y_continuous(
    sec.axis = sec_axis(~ fun.scale.up(.), name = "value (for AANY.simple)"),
    "value (for method 2)"
  )

p3

# now we have two y-axis

# lets make sure the x-axis labels don't overlap:

p3 <- p3 + scale_x_discrete(guide = guide_axis(n.dodge=3))

p3

#### plot 4 - biplot with color gradient ####

# lets look at the magnitude of CQ slope as function of watershed 
# percent developed and percent agriculture

# setting up plotting df:

df.plot <- df.setup %>% select(OLS.Slope.1s, R_DEVNLCD06, R_PLANTNLCD06) %>% 
  pivot_longer(cols = -c(R_DEVNLCD06, R_PLANTNLCD06))

# make plot:

p4 <- ggplot(df.plot, aes(x = R_DEVNLCD06, y = R_PLANTNLCD06, color = value)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "yellow", high = "red") +
  xlab('Percent Developed') +
  ylab('Percent Agriculture') +
  ggtitle('CQ slope as a function of percent Ag and Developed Land')

p4

# what if the size argument was inside aes?

ggplot(df.plot, aes(x = R_DEVNLCD06, y = R_PLANTNLCD06, color = value)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "yellow", high = "red")

# it gets added to the legend

# lets change the legend title:

p4$labels$colour <- 'CQ Slope' # dont ask me why its colour and not color...

p4

# we can also do:

p4 <- p4 + labs(color='CQ Slope!')

p4

#### plot 5 - boxplots and barplots ####

# lets look at boxplots of the different HSG watershed characteristics:

# set up plotting df:

df.plot <- df.setup %>% select(contains('HSG')) %>% 
  pivot_longer(cols = everything())

# make plot:

p5 <- ggplot(df.plot, aes(x = name, y = value))+
  geom_boxplot() +
  xlab('HSG class') +
  ylab('Watershed % in HSG')

p5

# lets rotate the axis labels:

p5 <- p5 + theme(axis.text.x=element_text(angle=40,hjust=1))

p5

# lets look at a barplot of just the means of each HSG:

# we can use the df.plot from the boxplots to set up the plotting df for the barplot:

df.plot <- df.plot %>% group_by(name) %>% summarize(mean = mean(value))

# make plot:

p5 <- ggplot(df.plot, aes(x = name, y = mean))+
  geom_bar(stat = "identity") +
  xlab('HSG class') +
  ylab('Avg. Watershed % in HSG')

p5

#### plot 6 - combining ggplots ####

# make list of plots:

plist <- list(p1,p3,p4,p5) # not including p2 because facets take awhile to plot

# make plot:

library(ggpubr)

p6<-ggpubr::ggarrange(plotlist = plist, ncol=2, nrow=2) #, common.legend = TRUE, legend="bottom") # arrange plots:

p6

# lets add a title to this plot:

p6<-annotate_figure(p6, top = text_grob("Mashup of plots from this lecture", color = "red", face = "bold", size = 14)) 

p6

