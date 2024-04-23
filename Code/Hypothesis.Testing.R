# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

library(FactoMineR)
library(factoextra)
library(leaps)
library(glmnet)
library(randomForest)
library(ggcorrplot)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(climateR)
library(ggpubr)
library(ggnewscale)
library(ggforce)
library(readxl)
library(MASS)
library(car)
library(CropScapeR)
library(FedData)
library(streamstats)
library(dataRetrieval)
library(sf)
library(raster)
library(terra)
library(kableExtra)
library(sfheaders)
library(mapview)
library(broom)
library(ggsignif)
library(ggpmisc)
library(segmented)
library(readxl)
library(tmap)
library(readxl)
library(zyp)
library(corrr)
library(gridExtra)
library(caret)
set.seed(123)
library(ggbiplot)
library(ggcorrplot)
library(tidyverse)

sf_use_s2(TRUE) # for sf

setTimeout(1000) # for streamstats api

meters_to_miles = 1/1609.334

####################### Functions #######################

source("Code/Ryan_functions.R")

####################### Goal of code #######################

####################### Workflow #######################

#### Read in Data ####

# read in l.resp.allpred for TP, TN, and SRP:

load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/TP.l.resp.allpred.Rdata")
l.TP <- l.resp.allpred
load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/TN.l.resp.allpred.Rdata")
l.TN <- l.resp.allpred
load("C:/Users/ryrug/OneDrive - SUNY ESF/Research/Processed_Data/SRP.l.resp.allpred.Rdata")
l.SRP <- l.resp.allpred




























#### Map of all Sites by consituent group ####

# get sites for each consituent:

TP.sites <- rownames(l.TP[[1]])
TN.sites <- rownames(l.TN[[1]])
SRP.sites <- rownames(l.SRP[[1]])

# find sites that are the intersection between all three consituents:

TP.TN <- intersect(TP.sites, TN.sites)
TP.TN.SRP <- intersect(TP.TN, SRP.sites)

TP.TN.SRP == TN.sites # all TN sites are in TP and SRP

# find sites that are intersection of TP and SRP:

TP.SRP <- intersect(TP.sites, SRP.sites)

TP.SRP == SRP.sites # all SRP sites are in TP

# set up df for mapping:

df.map <- data.frame(Name = TP.sites)

# add group column that identifes site with its consituent intersection:

df.map <- df.map %>% mutate(Group = ifelse(Name %in% TP.TN.SRP, 
                                           'TP.TN.SRP', 
                                           ifelse(Name %in% TP.SRP,
                                                  'TP.SRP',
                                                  'TP')))

# make map of sites for each consituent group:

# fun.map.DA(df.map$Name, df.map)

#



































#### CAFO count in Watershed ####

# read in CAFO locations and make sf dataframe

CAFOs<-read.csv('Raw_data/Copy of CAFO-list-from-NYSDEC March 2021 with lat long values.csv')%>%
  drop_na(longitude)%>%
  drop_na(latitude)%>%
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326)

# add CAFOS to watershed map:

# mapview(df.sf.NWIS.keep)+
#   mapview(CAFOs, zcol = 'SIZE')

# this file maps differently than the map on the DE website, but just going
# to proceed with this file for now

# load df of site name and watershed shapfile:

df.sf <- df.sf.NWIS.62 %>% filter(Name %in% TP.sites) %>% select(Name)

# determine how many CAFOs intersect each watershed:

df.sf$CAFO_count <- lengths(st_intersects(df.sf, CAFOs))

# make a map colored by CAFO count:

# mapview(df.sf, zcol = 'CAFO_count')

# create df of just site name and cafo count for merging later:

df.CAFO <- df.sf %>% st_set_geometry(NULL)

# 




































#### Set up df for H2 and H3 ####

# combine load response variables for each consituent:

l.load <- c(l.TP[3:6],l.TN[3:6],l.SRP[3:6])

names(l.load) <- c('TP.MAFWC', 'TP.FWAC', 'TP.AANY', 'TP.medC','TN.MAFWC', 'TN.FWAC', 'TN.AANY', 'TN.medC','SRP.MAFWC', 'SRP.FWAC', 'SRP.AANY', 'SRP.medC')

df.load <- bind_rows(l.load, .id = 'Term')

# set the rownames as first column:

df.load <- df.load %>% mutate(Name = gsub('\\..*', '', rownames(df.load)), .before = 1)

# add CAFO count:

df.load <- left_join(df.load, df.CAFO, by = 'Name')

# add constituent column to df.load:

df.load <- df.load %>% mutate(Consit = gsub("\\..*","",Term), .after = 1)

#





























#### H1 ####

# H1: Unless nearly 100% forested, most catchments display nutrient mobilization as Q increases

# combine slope response variables for each consituent:

l.slope <- c(l.TP[1:2],l.TN[1:2],l.SRP[1:2])

names(l.slope) <- c('TP.overal', 'TP.storm', 'TN.overal', 'TN.storm', 'SRP.overal', 'SRP.storm')

df.slope <- bind_rows(l.slope, .id = 'Term')

# make ggplot:

# ggplot(df.slope, aes(x = R_FORESTNLCD06, y = term))+
#   geom_point()+
#   geom_smooth(method = 'lm')+
#   facet_wrap(~Term, ncol = 2)+
#   stat_poly_eq(use_label(c("eq", "R2", 'p.value')))

#































#### H2 ####

# H2: Aggregate measures of the fraction developed land (ie. total ag, total urban) provide moderate explanation of watershed variations in nutrient load

# the predictor variables for this section are:

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06", 'FRAGUN')

# make correlation matrix plot for each consit:

l.plot <- df.load %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)

p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
           type = "lower",
           outline.col = "white",
           method = 'circle',
           lab = TRUE,
           p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
           insig = "blank",
           title = names(l.plot)[i]))

# ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')

# TN: remove WWTP.fraction and FRAGUN
# SRP/TP: remove HSG D and FRAGUN:

l.pred.vars <- list(SRP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN')],
                    TN = pred.vars[! pred.vars  %in% c('CAFO_count', 'FRAGUN')],
                    TP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN')])

# split df.load by consit:

l.load <- split(df.load, f=df.load$Consit)

# select the two term and just these predictors for each consit:

l.load <- lapply(names(l.load), \(i) l.load[[i]] %>% select(Term, term, l.pred.vars[[i]]))

# split this list into list ot list for each Term:

l.load <- lapply(l.load, \(i) split(i, f=i$Term))

# flatten list of list into just list of dfs:

l.load <- purrr::flatten(l.load)

# drop Term column from each df:

l.load <- lapply(l.load, \(i) i %>% select(-Term))

# build lm to predict:

l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))

# visualize model comparisons:

tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H2.html")

# lets look at VIF for these models:

l.VIF <- lapply(l.lm, vif)

df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% 
  pivot_longer(cols = -Term)

ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+
  geom_point()+
  geom_line()

# save l.load and l.lm as l.load.H2 and l.lm.H2:

l.load.H2 <- l.load

l.lm.H2 <- l.lm

#
























#### H3 ####

# H3: Incorporating measures of the fraction of developed land in close proximity to streams (i.e. land use within stream buffers) improves prediction of watershed variations in nutrient load

# the predictor variables for this section are:

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count',  'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06", "R_RIP100_PLANT", "R_RIP100_DEV")

# make correlation matrix plot for each consit:

l.plot <- df.load %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)

p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
                                      type = "lower",
                                      outline.col = "white",
                                      method = 'circle',
                                      lab = TRUE,
                                      p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
                                      insig = "blank",
                                      title = names(l.plot)[i]))

# ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')

# remove HSG_D, FRAGUN, R_PLANT, R_DEV:

l.pred.vars <- list(SRP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TN = pred.vars[! pred.vars  %in% c('FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")])

# split df.load by consit:

l.load <- split(df.load, f=df.load$Consit)

# select the two term and just these predictors for each consit:

l.load <- lapply(names(l.load), \(i) l.load[[i]] %>% select(Term, term, l.pred.vars[[i]]))

# split this list into list ot list for each Term:

l.load <- lapply(l.load, \(i) split(i, f=i$Term))

# flatten list of list into just list of dfs:

l.load <- purrr::flatten(l.load)

# drop Term column from each df:

l.load <- lapply(l.load, \(i) i %>% select(-Term))

# build lm to predict:

l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))

# visualize model comparisons:

tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H3.html")

# lets look at VIF for these models:

l.VIF <- lapply(l.lm, vif)

df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% 
  pivot_longer(cols = -Term)

ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+
  geom_point()+
  geom_line()

# save l.load and l.lm as l.load.H3 and l.lm.H3:

l.load.H3 <- l.load

l.lm.H3 <- l.lm

#

# RP improved Model adjR2:

l.H2.R2 <- lapply(l.lm.H2, \(i) summary(i)$adj.r.squared)
l.H3.R2 <- lapply(l.lm.H3, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H2.R2[[i]]-l.H3.R2[[i]])

l.compare.R2

#





























#### H3a: Include the WP ag/dev into H3 ####

l.pred.vars <- list(SRP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN')],
                    TN = pred.vars[! pred.vars  %in% c('FRAGUN')],
                    TP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN')])
l.load <- split(df.load, f=df.load$Consit)
l.load <- lapply(names(l.load), \(i) l.load[[i]] %>% select(Term, term, l.pred.vars[[i]]))
l.load <- lapply(l.load, \(i) split(i, f=i$Term))
l.load <- purrr::flatten(l.load)
l.load <- lapply(l.load, \(i) i %>% select(-Term))
l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))
tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H3a.html")
l.VIF <- lapply(l.lm, vif)
df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% 
  pivot_longer(cols = -Term)
ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+
  geom_point()+
  geom_line()

# RP improved model adjR2:

l.H3a.R2 <- lapply(l.lm, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H2.R2[[i]]-l.H3a.R2[[i]])

l.compare.R2

#

























#### Set up dfs for restrictions on sites ####

# load df.2 and df.3 site lists for TP, TN, and SRP:

load('Processed_Data/TP.df.2.sites.Rdata'); TP.restrict.2 <- unique(df.2$site_no)
load('Processed_Data/TP.df.3.sites.Rdata'); TP.restrict.3 <- unique(df.3$site_no)
load('Processed_Data/TN.df.2.sites.Rdata'); TN.restrict.2 <- unique(df.2$site_no)
load('Processed_Data/TN.df.3.sites.Rdata'); TN.restrict.3 <- unique(df.3$site_no)
load('Processed_Data/SRP.df.2.sites.Rdata'); SRP.restrict.2 <- unique(df.2$site_no)
load('Processed_Data/SRP.df.3.sites.Rdata'); SRP.restrict.3 <- unique(df.3$site_no)

# set up df for restriction 2 and 3:

l.restrict.2 <- mget(ls(pattern = 'restrict.2'))
l.restrict.3 <- mget(ls(pattern = 'restrict.3'))

# set names of list:

names(l.restrict.2) <- c('SRP', 'TN', 'TP')
names(l.restrict.3) <- c('SRP', 'TN', 'TP')

# split df.load by consit:

l.load <- split(df.load, f = df.load$Consit)

# filter each element of l.load by corresponding elements of restrict 2 and restirct 3 for 

l.load.restrict.2 <- lapply(names(l.load), \(i) l.load[[i]] %>% filter(Name %in% l.restrict.2[[i]]))
l.load.restrict.3 <- lapply(names(l.load), \(i) l.load[[i]] %>% filter(Name %in% l.restrict.3[[i]]))

# combine list back into single df:

df.load.restrict.2 <- bind_rows(l.load.restrict.2)
df.load.restrict.3 <- bind_rows(l.load.restrict.3)

# 






























#### H2 with restriction 2 ####

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06", 'FRAGUN')
l.plot <- df.load.restrict.2 %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)
p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
                                      type = "lower",
                                      outline.col = "white",
                                      method = 'circle',
                                      lab = TRUE,
                                      p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
                                      insig = "blank",
                                      title = names(l.plot)[i]))
ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')
pred.vars <- c('CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06")
l.load <- split(df.load.restrict.2, f=df.load.restrict.2$Term)
l.load <- lapply(l.load, \(i) i %>% select(Term, term, pred.vars)) # sticking with result 1 for now...
l.load <- lapply(l.load, \(i) i %>% select(-Term))
l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))
tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H2.restriction.2.html")
l.VIF <- lapply(l.lm, vif)
df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% 
  pivot_longer(cols = -Term)
ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+
  geom_point()+
  geom_line()+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))

# H2 w/ r2 improved H2 w/ r1:

l.H2.restrict.2.R2 <- lapply(l.lm, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H2.R2[[i]]-l.H2.restrict.2.R2[[i]])

l.compare.R2

#






























#### H2 with restriction 3 ####

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06", 'FRAGUN')
l.plot <- df.load.restrict.3 %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)
p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
                                      type = "lower",
                                      outline.col = "white",
                                      method = 'circle',
                                      lab = TRUE,
                                      p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
                                      insig = "blank",
                                      title = names(l.plot)[i]))
ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')
pred.vars <- c('WWTP.fraction', 'R_PLANTNLCD06', "R_DEVNLCD06")
l.load <- split(df.load.restrict.2, f=df.load.restrict.2$Term)
l.load <- lapply(l.load, \(i) i %>% select(Term, term, pred.vars)) # sticking with result 1 for now...
l.load <- lapply(l.load, \(i) i %>% select(-Term))
l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))
tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H2.restriction.3.html")
l.VIF <- lapply(l.lm, vif)
df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% 
  pivot_longer(cols = -Term)
ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+
  geom_point()+
  geom_line()+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))

# H2 w/ r3 improved H2 w/ r1

l.H2.restrict.3.R2 <- lapply(l.lm, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H2.R2[[i]]-l.H2.restrict.3.R2[[i]])

l.compare.R2

# H2 w/ r3 improved H2 w/ r2

l.compare.R2 <- sapply(seq_along(l.H2.restrict.2.R2), \(i) l.H2.restrict.2.R2[[i]]-l.H2.restrict.3.R2[[i]])

l.compare.R2

#

























#### H3 with restriction 2 ####

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count',  'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06", "R_RIP100_PLANT", "R_RIP100_DEV")
l.plot <- df.load.restrict.2 %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)
p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
                                      type = "lower",
                                      outline.col = "white",
                                      method = 'circle',
                                      lab = TRUE,
                                      p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
                                      insig = "blank",
                                      title = names(l.plot)[i]))
# ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')
l.pred.vars <- list(SRP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TN = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")])
l.load <- split(df.load.restrict.2, f=df.load.restrict.2$Consit)
l.load <- lapply(names(l.load), \(i) l.load[[i]] %>% select(Term, term, l.pred.vars[[i]]))
l.load <- lapply(l.load, \(i) split(i, f=i$Term))
l.load <- purrr::flatten(l.load)
l.load <- lapply(l.load, \(i) i %>% select(-Term))
l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))
tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H3.restriction.2.html")
l.VIF <- lapply(l.lm, vif)
df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% 
  pivot_longer(cols = -Term)
ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+
  geom_point()+
  geom_line()

# RP improved model adjR2:

l.H3.restrict.2.R2 <- lapply(l.lm, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H2.restrict.2.R2[[i]]-l.H3.restrict.2.R2[[i]])

l.compare.R2

# H3 w/ r2 improved H3 w/ r1 

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H3.R2[[i]]-l.H3.restrict.2.R2[[i]])

l.compare.R2

#
































#### H3 with restriction 3 ####

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count',  'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06", "R_RIP100_PLANT", "R_RIP100_DEV")
l.plot <- df.load.restrict.3 %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)
p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
                                      type = "lower",
                                      outline.col = "white",
                                      method = 'circle',
                                      lab = TRUE,
                                      p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
                                      insig = "blank",
                                      title = names(l.plot)[i]))
ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')
l.pred.vars <- list(SRP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TN = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TP = pred.vars[! pred.vars  %in% c('HSG_D', 'FRAGUN', 'CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06")])
l.load <- split(df.load.restrict.3, f=df.load.restrict.3$Consit)
l.load <- lapply(names(l.load), \(i) l.load[[i]] %>% select(Term, term, l.pred.vars[[i]]))
l.load <- lapply(l.load, \(i) split(i, f=i$Term))
l.load <- purrr::flatten(l.load)
l.load <- lapply(l.load, \(i) i %>% select(-Term))
l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))
tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H3.restriction.3.html")
l.VIF <- lapply(l.lm, vif)
df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% 
  pivot_longer(cols = -Term)
ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+
  geom_point()+
  geom_line()

# RP improved Model adjR2:

l.H3.restrict.3.R2 <- lapply(l.lm, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.restrict.3.R2), \(i) l.H2.restrict.3.R2[[i]]-l.H3.restrict.3.R2[[i]])

l.compare.R2

# H3 w/ r3 improved H3 w/ r1 

l.compare.R2 <- sapply(seq_along(l.H3.restrict.3.R2), \(i) l.H3.R2[[i]]-l.H3.restrict.3.R2[[i]])

l.compare.R2

# H3 w/ r3 improved H3 w/ r2: 

l.compare.R2 <- sapply(seq_along(l.H3.restrict.3.R2), \(i) l.H3.restrict.2.R2[[i]]-l.H3.restrict.3.R2[[i]])

l.compare.R2

#


























#### Set up df for overlapping sites across TP, TN, SRP ####

# for R1, filter SRP and TP sites down to just TN sites:

sites.to.keep <- unique((df.load %>% filter(Consit == 'TN'))$Name)

df.load.overlap <- df.load %>% filter(Name %in% sites.to.keep)

#






























#### H2 with overlapping sites + R1 ####

# H2 workflow:

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count', 'R_PLANTNLCD06', "R_DEVNLCD06", 'FRAGUN')
l.plot <- df.load.overlap %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)
p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
                                      type = "lower",
                                      outline.col = "white",
                                      method = 'circle',
                                      lab = TRUE,
                                      p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
                                      insig = "blank",
                                      title = names(l.plot)[i]))
ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')
l.pred.vars <- list(SRP = pred.vars[! pred.vars  %in% c( 'FRAGUN')],
                    TN = pred.vars[! pred.vars  %in% c('FRAGUN')],
                    TP = pred.vars[! pred.vars  %in% c('FRAGUN')])
l.load <- split(df.load.overlap, f=df.load.overlap$Consit)
l.load <- lapply(names(l.load), \(i) l.load[[i]] %>% select(Term, term, l.pred.vars[[i]]))
l.load <- lapply(l.load, \(i) split(i, f=i$Term))
l.load <- purrr::flatten(l.load)
l.load <- lapply(l.load, \(i) i %>% select(-Term))
l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))
tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H2.overlapping.sites.html")
l.VIF <- lapply(l.lm, vif)
df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% pivot_longer(cols = -Term)
ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+geom_point()+geom_line()

# compare adjR2 with H2 w/ R1:

l.H2.overlap.R2 <- lapply(l.lm, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H2.R2[[i]]-l.H2.overlap.R2[[i]])

l.compare.R2

#




























#### H3 with overlapping sites + R1 ####

# H3 workflow:

pred.vars <- c('HSG_D', 'WWTP.fraction', 'CAFO_count',  'FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06", "R_RIP100_PLANT", "R_RIP100_DEV")
l.plot <- df.load.overlap %>% group_by(Consit) %>% distinct(Name, .keep_all = T) %>% ungroup() %>%  split(., f = .$Consit)
p.list <- lapply(1:3, \(i) ggcorrplot(cor(l.plot[[i]] %>% select(pred.vars)), 
                                      type = "lower",
                                      outline.col = "white",
                                      method = 'circle',
                                      lab = TRUE,
                                      p.mat = cor_pmat(l.plot[[i]] %>% select(pred.vars)),
                                      insig = "blank",
                                      title = names(l.plot)[i]))
ggarrange(plotlist = p.list, ncol = 3, common.legend = T, legend = 'right')
l.pred.vars <- list(SRP = pred.vars[! pred.vars  %in% c('FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TN = pred.vars[! pred.vars  %in% c('FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")],
                    TP = pred.vars[! pred.vars  %in% c('FRAGUN', 'R_PLANTNLCD06', "R_DEVNLCD06")])
l.load <- split(df.load.overlap, f=df.load.overlap$Consit)
l.load <- lapply(names(l.load), \(i) l.load[[i]] %>% select(Term, term, l.pred.vars[[i]]))
l.load <- lapply(l.load, \(i) split(i, f=i$Term))
l.load <- purrr::flatten(l.load)
l.load <- lapply(l.load, \(i) i %>% select(-Term))
l.lm <- lapply(l.load, \(i) lm(term ~ ., data = i))
tab_model(l.lm, dv.labels = names(l.lm), title = paste('Comparison of MLR models'), file="H3.overlapping.sites.html")
l.VIF <- lapply(l.lm, vif)
df.VIF <- bind_rows(l.VIF, .id = 'Term') %>% pivot_longer(cols = -Term)
ggplot(df.VIF, aes(x = Term, y = value, color = name, group = name))+geom_point()+geom_line()

# compare with H2 :

l.H3.overlapping.R2 <- lapply(l.lm, \(i) summary(i)$adj.r.squared)

l.compare.R2 <- sapply(seq_along(l.H2.R2), \(i) l.H2.R2[[i]]-l.H3.overlapping.R2[[i]])

l.compare.R2



























