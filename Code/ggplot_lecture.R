# load packages:

library(tidyverse)

# load data:

load('Processed_Data/df.setup.Rdata')

# plot 1: univariate plot of one of the response variables and one of the predictors:

names(df.setup)

p <- df.setup %>% select(2,23) %>% 
  ggplot(., aes(x = R_PLANTNLCD06, y = OLS.Int.1s))+
  geom_point()+
  geom_smooth(method = 'lm')

p

# add trendline equation to plot: https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph

# option 1: with function:

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df.setup %>% select(2,23)), parse = TRUE)

library(ggpmisc)

+
  stat_poly_line() +
  stat_poly_eq()
