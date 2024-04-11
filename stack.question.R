library(caret)
library(tidyverse)

# load data

load('Processed_Data/df1.Rdata')

# set up different CV methods:

trC.cv <- trainControl(method = "cv", number = 10)
trC.repeatedcv <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
trC.LGOCV <- trainControl(method = "LGOCV", p = 0.8)
trC.boot <- trainControl(method = "boot")
trC.LOOCV <- trainControl(method = "LOOCV")

# combine CV methods into list:

l.trCon <- list(trC.boot, trC.cv, trC.repeatedcv, trC.LGOCV, trC.LOOCV)
names(l.trCon) <- c("trC.boot", "trC.cv", "trC.repeatedcv", "trC.LGOCV", "trC.LOOCV")

# function to compare ncomp across different seeds/cv methods:

fun.compare.plsr.models <- function(df, seeds){
  df.i <- data.frame(seed = NA, CV.method = NA, ncomp = NA) # initialize df to cbind to for each i loop:
  for (i in seq_along(seeds)){ # loop through the seeds:
    v.CV.method <- NA # initialize vectors for Cv method and ncomp for each j loop:
    v.ncomp <- NA
    for(j in seq_along(l.trCon)){ # loop through CV methods:
      set.seed(i) # set seed:
      model.j <- train(term ~., data = df, # build model using CV method corresponding to j loop:
                       method = 'pls', 
                       scale = TRUE, 
                       trControl = l.trCon[[j]], 
                       tuneGrid = data.frame(ncomp = c(1:30)))
      v.CV.method[j] <- names(l.trCon)[j] # set vector element for j loop:
      v.ncomp[j] <- model.j$bestTune[1,1]
    }
    df.j <- data.frame(seed = i, CV.method = v.CV.method, ncomp = v.ncomp) # create df with resulting vectors from j loop:
    df.i <- rbind(df.i, df.j) # cbind j loop df to overall df:
  }
  df.i <- df.i[-1,] # remove first row since it was NA:
  return(df.i)
}

# run function (takes a minute to run):

system.time(df.compare10 <- fun.compare.plsr.models(df = df1, seeds = 1:10)) # 95 seconds

# look at results:

ggplot(df.compare10, aes(x = seed, y = ncomp, color = CV.method)) +
  geom_line()




