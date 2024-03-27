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
library(tidyverse)

sf_use_s2(TRUE) # for sf

setTimeout(1000) # for streamstats api

meters_to_miles = 1/1609.334

####################### Functions #######################

source("Code/Ryan_functions.R")

####################### Goal of code #######################

# 1) Process the NWIS database for the CQ analysis
# 2) Run through CQ metrics/watershed attribute correlation

####################### Workflow #######################

# TP pcode:

pcode<-'00665'

ncode<-'TP'























#### NWIS Query ####

# this code adapted from NYS_site_ranking.R

# Download the metadata for sites with daily flow data. To do this:

# Use dataRetrieval::whatNWISdata to query USGS gauge sites with daily flow data in NY

# df.NWIS.Q_sites<- whatNWISdata(stateCd = 'NY', parameterCd = "00060")

# write.csv(df.NWIS.Q_sites, "C:/PhD/CQ/Raw_Data/df.NWIS.Q_sites.csv", row.names=FALSE)

df.NWIS.Q_sites <- read.csv("C:/PhD/CQ/Raw_Data/df.NWIS.Q_sites.csv", colClasses = c(site_no = "character"))

# clean up this dataframe

# there are duplicates of site numbers with different numbers of daily flow observations
# so group by the site no numbers and just keep the observation with the highest number of flow observations
# I realize now that this is proably not needed since I am grouping by the site number so the resulting number used in the
# raw flow data download will have what it has, but I am just goin to keep this (aka i could just use distinct?)
# also rename some columns and convert them to numeric:
# also filter to just the sites with >= 20 nflow days since this is the number restriction on the number of samples:

df.NWIS.Q_sites<-df.NWIS.Q_sites%>%
  group_by(site_no)%>%
  slice(which.max(count_nu))%>%
  rename(latitude = dec_lat_va, longitude = dec_long_va, nflowdays = count_nu, begin_date_flow = begin_date, end_date_flow = end_date)%>%
  mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))%>%
  filter(nflowdays>=20)

# NOTE: whatNWISsites(stateCd = "NY", parameterCd = "00060") returns a dataframe with rnow = 1465
# while whatNWISdata(stateCd = "NY", parameterCd = "00060") returns a dataframe with nrow = 2249 ?!?!
# this is because there are duplicates. Using the group_by, slice functions gives a dataframe with
# nrow = 1452 for the whatNWISdata thing.

# From the sites with daily flow data, determine which sites also have over 99 discrete samples from one or more of the constituents of interests. To do this:

# use function I created (sourced) to get a dataframe of sites for just TP with #samples = 20 threshold:

# df.NWIS.TP_sites<-fun.df.Pair_consit_flow(pcode, df.NWIS.Q_sites, n_samples = 20, state = 'NY')

# write.csv(df.NWIS.TP_sites, "Raw_Data/df.NWIS.TP_sites.csv", row.names=FALSE)

df.NWIS.TP_sites<-read.csv("Raw_Data/df.NWIS.TP_sites.csv", colClasses = c(site_no = "character"))

# download the raw daily flow data for these sites. To do this:

# readNWISdv can be used to download a single dataframe with all the raw flow data for a vector of gauge numbers (this takes a looong time):

# df.NWIS.Q<-readNWISdv(siteNumbers = df.NWIS.TP_sites$site_no, parameterCd = '00060', startDate = "", endDate = "", statCd = "00003")

# write.csv(df.NWIS.Q, "Raw_Data/df.NWIS.Q.csv", row.names=FALSE)

df.NWIS.Q<-read.csv("C:/PhD/CQ/Raw_Data/df.NWIS.Q.csv", colClasses = c(site_no = "character"))

# download the raw discrete TP sample data:

# df.NWIS.TP<-readNWISqw(siteNumbers = df.NWIS.TP_sites$site_no, parameterCd = pcode)

# write.csv(df.NWIS.TP, "Raw_Data/df.NWIS.TP.csv", row.names=FALSE)

df.NWIS.TP<-read.csv("Raw_Data/df.NWIS.TP.csv", colClasses = c(site_no = "character"))

#























#### Building CQ df ####

# join the TP with the flow df:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP, df.NWIS.Q, by=c("site_no"="site_no", "sample_dt"="Date"))

# take average of multiple samples on the same day at the same site:

df.NWIS.TP_CQ<-df.NWIS.TP_CQ%>%
  group_by(site_no, sample_dt)%>%
  summarise_at(vars(result_va, X_00060_00003), funs(mean(., na.rm=TRUE)))%>%
  ungroup()

# remove sites where there are not >= 20 CQ pairs:

df.NWIS.TP_CQ<-df.NWIS.TP_CQ%>%
  drop_na(result_va, X_00060_00003)%>%
  group_by(site_no)%>%
  mutate(n=n())%>%
  ungroup()%>%
  filter(n>=20)%>%
  select(-n) # remove n column because going to get added back in next steps

# add column for number of paired observations by rank: to do this:

# create column for n_sample_rank:

temp<-df.NWIS.TP%>%group_by(site_no)%>%
  summarise(n=n())%>%
  arrange(desc(n))%>%
  mutate(n_sample_rank=rank(-n, ties.method='first'))

# merge this df with df.NWIS.TP_CQ and arrange by the new column:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP_CQ,temp, by='site_no')%>%
  arrange(n_sample_rank)

# next step is to determine which sites have draiange area listed
# if no draiange area is listed, there is no way to determine if delineation
# is correct. to do this:

# download the drainage areas from site metadata using readNWISsite

# df.NWIS.TP_site_metadata<-readNWISsite(siteNumbers = unique(df.NWIS.TP_CQ$site_no))

# write.csv(df.NWIS.TP_site_metadata, "Raw_Data/df.NWIS.TP_site_metadata.csv", row.names=FALSE)

df.NWIS.TP_site_metadata<-read.csv("Raw_Data/df.NWIS.TP_site_metadata.csv", colClasses = c(site_no = "character"))

# then select just the site number and DA column:

df.DA<-df.NWIS.TP_site_metadata%>%
  select(site_no, drain_area_va)

# finally merge with df.NWIS.TP_CQ and create a new Q column with area normalized flows (not worrying about units right now): 
# Note: Some sites returned NA on draiange areas in readNWISsite:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP_CQ, df.DA, by = 'site_no')%>%
  mutate(Q_yield = X_00060_00003/drain_area_va)%>%
  drop_na(drain_area_va)

# length(unique(df.NWIS.TP_CQ$site_no))  121 sites

















































#### Tradeoff matrix ####

# # build a matrix of the number of TP samples as a function of watershed size. To do this:
# 
# # I already have a df with paired CQ observations and DA, just need to get the distinct sites:
# 
# TP_sites<-df.NWIS.TP_CQ%>%distinct(site_no, .keep_all = T)
# 
# # set up df to populate:
# 
# m<-data.frame(Min_num_samples  = c(20,50,75,100,200), '25' = NA, '50' = NA, '100' = NA, '150'=NA, '250'=NA, '500'=NA, '1000'=NA, 'Unlimited'=NA)
# 
# # set up variables for thresholds for number of samples and DA size:
# 
# n_sam<- c(20,50,75,100,200)-1
# 
# min_DA<- c(25,50,100,150,250,500,1000)
# 
# # loop through number of samples (rows):
# 
# # i<-1
# 
# for (i in seq_along(n_sam)){
#   
#   temp.i<-TP_sites%>%filter(n>n_sam[i])
#   
#   m$Unlimited[i]<-dim(temp.i)[1]
# 
#   # j<-2
#   
#   # loop through the size of the DA (columns)
#   for(j in  seq_along(min_DA)){
#     
#     temp.j<-temp.i%>%filter(drain_area_va<=min_DA[j])
#     
#     m[i,j+1]<-dim(temp.j)[1]
#     
#   }
#   
# }
# 
# m
# 
# #


























#### Apply filters to sites ####

# number of sites with >=20 paired CQ observations:

length(unique(df.NWIS.TP_CQ$site_no))

# 1) sites with samples after 2001:

# note that you need to rerun the 
# number of sample filter again because some sites 
# have data pre and post 2001

temp1<-df.NWIS.TP_CQ%>%filter(sample_dt >= as.Date('2001-01-10'))%>%
  group_by(site_no)%>%
  mutate(n=n())%>%
  filter(n>=20)

length(unique(temp1$site_no))

# 64 sites

# 2) sites not on long island

# use latitude to filter:

temp2<-temp1%>%filter(site_no %in% filter(df.NWIS.TP_site_metadata, dec_lat_va >40.9364)$site_no)

length(unique(temp2$site_no))

# 63 sites

# 3) sites in gauges 2

# read in G2: to do this

# l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm.xlsx") # read in all sheets using function:
# l.G2[[27]]<-NULL # remove the last element (does notcomtain useful info)
# df.G2<-reduce(l.G2, full_join, by = "STAID") # convert to a single df
# save(df.G2, l.G2, file = 'Processed_Data/df.G2.Rdata') # save df.G2:
load('Processed_Data/df.G2.Rdata') # load df.G2:

# filter on G2:

temp3<-temp2%>%filter(site_no %in% df.G2$STAID)

length(unique(temp3$site_no))

# 42 sites

# one last test: redoing the filters but with no removing date restriction:

# number of sites in gauges 2:

x1<-filter(df.NWIS.TP_CQ, site_no %in%df.G2$STAID)

length(unique(x1$site_no))

# 83 sites

# number of sites in gauges 2 and not on LI:

x2<-x1%>%filter(site_no %in% filter(df.NWIS.TP_site_metadata, dec_lat_va >40.9364)$site_no)

length(unique(x2$site_no))

# 73 sites



























#### ! Finalizing df.TP_CQ to move on in code ####

# first pass 

df.TP_CQ<-temp3

# save:

# save(df.TP_CQ,file = 'Processed_Data/TP.Rdata')

# second pass

df.TP_CQ.SP<-temp2

# save:

# save(df.TP_CQ.SP, file = 'Processed_Data/TP.SP.Rdata')

# going to write over df.TP_CQ for this code to work smoothly with the these sites

df.TP_CQ<-df.TP_CQ.SP

#



















#### Fitting Breakpoints to CQ curves #### 

# use funciton in lapply to get list of df of segmented results:

l_Seg<-lapply(seq_along(unique(df.TP_CQ$site_no)), \(i) fun.CQ.BP(i, df.TP_CQ))

# now combine the list of dfs of the breakpoint analysis results (with fitted values, intercepts and slopes) into a single df:

df_Seg<-bind_rows(l_Seg)

# use replace and the BP yes column with a conditional statement to set the breakpoint Q and C column rows to NA, as to not plot the segmeneted line if davies test was false:

df_Seg<-df_Seg%>%
  mutate(across(c(1,2), ~replace(., BP_yes == 'FALSE', NA)))

# add a column for number of samples:

temp<-df_Seg%>%
  group_by(site)%>%
  summarise(n=n())%>%
  arrange(desc(n))%>%
  mutate(n_sample_rank=rank(-n, ties.method='first'))

df_Seg<-left_join(df_Seg, temp%>%select(site, n_sample_rank), by = 'site')%>%
  arrange(n_sample_rank)

df.Seg <- df_Seg

save(df.Seg, file = 'Processed_Data/df.Seg.Rdata')

# ready to plot:

# ggplot(df_Seg, aes(x = log(Q_real), y = log(C)))+
#   geom_point()+
#   geom_smooth(method = 'lm')+
#   geom_line(aes(x = Q, y = Seg_C), color = 'tomato')+
#   geom_vline(xintercept =log(median(temp$Q_real)), color='red')
#   facet_wrap(dplyr::vars(site), scales = 'free') #+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )


# one last thing: lets look at a map of these:

# map.NWIS.TP_sites<-df.NWIS.TP_site_metadata%>%
#   rename(longitude=8,latitude=7)%>%
#   drop_na(latitude,longitude)%>%
#   st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)%>%
#   mutate(site_id= paste(site_no, station_nm), n = NA, .before = 1)%>%
#   select(c(1:3))

# mapview(map.NWIS.TP_sites)

# save(map.NWIS.TP_sites, file = 'C:/PhD/CQ/Processed_Data/map.NWIS.TP_sites.Rdata')






























#### Setting up Response Variables ####

#### ~ Set up df for this section: ####

df.RP <- df.TP_CQ%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003, n_sample_rank)%>%
  mutate(sample_dt = as.Date(sample_dt), log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))

l.RP<-split(df.RP, f=df.RP$Name)

df.Q<-df.NWIS.Q%>% # this gives df of complete hydrographs for each site and includes columns to group by calender day later on
  filter(site_no %in% df.RP$Name)%>%
  mutate(Date = as.Date(Date))%>%
  mutate(year = year(Date), Date_AAH = as.Date(format(Date, "%m-%d"), "%m-%d"))

l.Q<-df.Q%>% # this gives list of dfs of average annual hydrograph for each site
  group_by(site_no, Date_AAH)%>%
  summarize(mean_Q = mean(X_00060_00003, na.rm= T))%>%
  split(., f = .$site_no)

l.Q.2<-df.Q%>% # this gives list of dfs of complete hydrographs for each site
  split(., f = .$site_no)

l.Q.v<-lapply(l.Q, \(i) i$mean_Q) # list of vectors of just the mean annual hydrographs

l.DA<-df.NWIS.TP_site_metadata %>% filter(site_no %in% names(l.Q))%>%select(site_no, drain_area_va)%>%split(., .$site_no)%>%lapply(., \(i) i$drain_area_va)

unit_conversion=28.3168*86400*(1/1000)*(1/1000)*(1/258.999) # mg/L * ft^3/sec * 28.3168 L/ft3 * 86400 sec /day * 1g/1000mg * 1kg/1000g * 1/mi^2 * 1mi^2/258.999 ha = kg/ha/day

# 1) CQ metrics and associated Average Annual Nutrient Yields:

# note: steps for calculating AANY:

# the first method employed used an average annual hydrograph. This worked for the
# 1 slope and 2 slope median Q segmentations, but a different method is need for the
# hydrosep and median C segmentataions.

# the second method will estimate C for everday of the sites POR, then an 
# average annual yieldograph will be computed, and then summed up to get AANY

#### ~ 1.1) OLS intercept, slope, AANY for single slope model: ####

# slope and intercept:

df.OLS<-df.RP%>%
  group_by(Name)%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.Int.1s = OLS.co[1], 
            OLS.Slope.1s = OLS.co[2])
  }) %>%
  ungroup

# AANY:

l.lm <- df.RP%>%split(., f=.$Name)%>%lapply(., \(i) lm(log_C ~ log_Q, i))
l.BCF<-lapply(l.lm, \(i) exp((summary(i)$sigma^2)/2))

# first approach:

l.pred.C<-lapply(names(l.lm), \(i) l.BCF[[i]]*exp(as.numeric(predict(l.lm[[i]], data.frame(log_Q = log(l.Q[[i]]$mean_Q))))))%>%purrr::set_names(names(l.Q))
l.Yield<-lapply(names(l.lm), \(i) l.pred.C[[i]]*l.Q.v[[i]]*(1/l.DA[[i]])*unit_conversion)%>%purrr::set_names(names(l.Q))
df.Yield <- lapply(l.Yield, sum)%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'OLS.AANY.1s') # units of kg/ha/year

df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')

# second approach:

l.pred.C<-lapply(names(l.lm), \(i) data.frame(Date_AAH = l.Q.2[[i]]$Date_AAH, Q = l.Q.2[[i]]$X_00060_00003, pred.C = l.BCF[[i]]*exp(as.numeric(predict(l.lm[[i]], data.frame(log_Q = log(l.Q.2[[i]]$X_00060_00003)))))))%>%purrr::set_names(names(l.Q.2))
l.Yield<-lapply(names(l.lm), \(i) l.pred.C[[i]]%>%mutate(pred.yield = pred.C*Q*(1/l.DA[[i]])*unit_conversion))%>%purrr::set_names(names(l.Q.2))
l.Yield.AAH <- lapply(l.Yield, \(i) i%>%group_by(Date_AAH)%>%summarise(mean_daily_yield=mean(pred.yield)))
df.Yield <- lapply(l.Yield.AAH, \(i) sum(i$mean_daily_yield, na.rm = T))%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'OLS.AANY.1s.method2') # units of kg/ha/year

df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')

#### ~ 1.2) OLS intercept, slope, AANY for two slope model - threshold by median flow rate ####

# slope and intercept: create two df.OLS, one for pre and post median flow rate:

temp.pre<-df.RP%>%
  group_by(Name)%>%
  filter(X_00060_00003<median(X_00060_00003))%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.Int.2s_medQ_pre = OLS.co[1], 
            OLS.Slope.2s_medQ_pre = OLS.co[2])
  }) %>%
  ungroup

temp.post<-df.RP%>%
  group_by(Name)%>%
  filter(X_00060_00003>=median(X_00060_00003))%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.Int.2s_medQ_post = OLS.co[1], 
            OLS.Slope.2s_medQ_post = OLS.co[2])
  }) %>%
  ungroup

# combine pre and post:

temp <- left_join(temp.pre,temp.post,by='Name')

# merge with df.OLS:

df.OLS <- left_join(df.OLS,temp,by='Name')%>%as.data.frame()

# AANY:

df.median.Q<-df.TP_CQ%>%group_by(site_no)%>%summarise(median_Q=median(X_00060_00003))

l.lm.pre <- df.RP%>%
  filter(X_00060_00003<median(X_00060_00003))%>%
  split(., f=.$Name)%>%
  lapply(., \(i) lm(log_C ~ log_Q, i))

l.lm.post <- df.RP%>%
  filter(X_00060_00003>=median(X_00060_00003))%>%
  split(., f=.$Name)%>%
  lapply(., \(i) lm(log_C ~ log_Q, i))

# make sure list orders are the same:

names(l.lm.pre)==names(l.lm.post)
names(l.lm.pre)==names(l.RP)
names(l.lm.pre)==df.median.Q$site_no

# back to workflow:

l.BCF.pre<-lapply(l.lm.pre, \(i) exp((summary(i)$sigma^2)/2))
l.BCF.post<-lapply(l.lm.post, \(i) exp((summary(i)$sigma^2)/2))

# first approach:

l.pred.C<-lapply(seq_along(l.lm.pre), \(i) l.Q[[i]]%>%mutate(pred.C = ifelse(mean_Q < df.median.Q$median_Q[i],
                     l.BCF.pre[[i]]*exp(as.numeric(predict(l.lm.pre[[i]], data.frame(log_Q = log(l.Q[[i]]$mean_Q))))),
                     l.BCF.post[[i]]*exp(as.numeric(predict(l.lm.post[[i]], data.frame(log_Q = log(l.Q[[i]]$mean_Q)))))
                     )))%>%
  purrr::set_names(names(l.Q))

l.Yield<-lapply(names(l.lm), \(i) l.pred.C[[i]]$pred.C*l.Q.v[[i]]*(1/l.DA[[i]])*unit_conversion)%>%purrr::set_names(names(l.Q))

df.Yield <- lapply(l.Yield, sum)%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'OLS.AANY.2s.medQ') # units of kg/ha/year

df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')

# second approach:

l.pred.C<-lapply(seq_along(l.lm.pre), \(i) l.Q.2[[i]]%>%mutate(pred.C = ifelse(X_00060_00003 < df.median.Q$median_Q[i],
                                                                l.BCF.pre[[i]]*exp(as.numeric(predict(l.lm.pre[[i]], data.frame(log_Q = log(l.Q.2[[i]]$X_00060_00003))))),
                                                                l.BCF.post[[i]]*exp(as.numeric(predict(l.lm.post[[i]], data.frame(log_Q = log(l.Q.2[[i]]$X_00060_00003)))))
)))%>%
  purrr::set_names(names(l.Q.2))

l.Yield<-lapply(names(l.lm), \(i) data.frame(Date_AAH =l.pred.C[[i]]$Date_AAH, pred.yield = l.pred.C[[i]]$pred.C*l.Q.2[[i]]$X_00060_00003*(1/l.DA[[i]])*unit_conversion))%>%purrr::set_names(names(l.Q.2))
l.Yield.AAH <- lapply(l.Yield, \(i) i%>%group_by(Date_AAH)%>%summarise(mean_daily_yield=mean(pred.yield)))
df.Yield <- lapply(l.Yield.AAH, \(i) sum(i$mean_daily_yield, na.rm = T))%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'OLS.AANY.2s.medQ.method2') # units of kg/ha/year

df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')

# check: I want to plot CQ curves with pre and post median OLS lines. to do this:

# use models in ifelse to predict concentration column:

l.test<-lapply(seq_along(l.lm.pre), \(i) l.RP[[i]]%>%mutate(Seg_C.medQ = ifelse(Q < df.median.Q$median_Q[i], 
          predict(l.lm.pre[[i]], data.frame(log_Q = l.RP[[i]]$log_Q)),
          predict(l.lm.post[[i]], data.frame(log_Q = l.RP[[i]]$log_Q))
))%>%
  mutate(Seg_type = ifelse(Q < df.median.Q$median_Q[i],
                           'Pre_Q_med',
                           'Post_Q_med'))
)%>%purrr::set_names(names(l.RP))

# check with plot:

i <- 2

ggplot(l.test[[i]], aes(x = log_Q, y = log_C))+
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = log_Q, y = Seg_C.medQ, color = Seg_type))

df.OLS[i,6:9]

# looks good!

#### ~ 1.3) OLS intercept, slope, AANY for two slope model - threshold by median concentration ####

# not going to run since I cant figure out how to calculate AANY using median C...

# # create two df.OLS, one for pre and post median conc:
# 
# temp.pre<-df.RP%>%
#   group_by(Name)%>%
#   filter(C<median(C))%>%
#   do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
#   summarize(., OLS.Int.2s_medC_pre = OLS.co[1], 
#             OLS.Slope.2s_medC_pre = OLS.co[2])
#   }) %>%
#   ungroup
# 
# temp.post<-df.RP%>%
#   group_by(Name)%>%
#   filter(C>=median(C))%>%
#   do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
#   summarize(., OLS.Int.2s_medC_post = OLS.co[1], 
#             OLS.Slope.2s_medC_post = OLS.co[2])
#   }) %>%
#   ungroup
# 
# # combine pre and post:
# 
# temp <- left_join(temp.pre,temp.post,by='Name')
# 
# # merge with df.OLS:
# 
# df.OLS <- left_join(df.OLS,temp,by='Name')%>%as.data.frame()
# 
# # AANY:
# 
# df.median.C<-df.TP_CQ%>%group_by(site_no)%>%summarise(median_C=median(result_va))
# 
# l.lm.pre <- df.RP%>%
#   filter(C<median(C))%>%
#   split(., f=.$Name)%>%
#   lapply(., \(i) lm(log_C ~ log_Q, i))
# 
# l.lm.post <- df.RP%>%
#   filter(C>=median(C))%>%
#   split(., f=.$Name)%>%
#   lapply(., \(i) lm(log_C ~ log_Q, i))
# 
# # make sure list orders are the same:
# 
# names(l.lm.pre)==names(l.lm.post)
# names(l.lm.pre)==names(l.RP)
# names(l.lm.pre)==df.median.Q$site_no
# 
# # back to workflow:
# 
# l.BCF.pre<-lapply(l.lm.pre, \(i) exp((summary(i)$sigma^2)/2))
# l.BCF.post<-lapply(l.lm.post, \(i) exp((summary(i)$sigma^2)/2))
# 
# l.pred.C<-lapply(seq_along(l.lm.pre), \(i) l.Q[[i]]%>%mutate(pred.C = ifelse(mean_Q < df.median.Q$median_Q[i],
#                                                                              l.BCF.pre[[i]]*exp(as.numeric(predict(l.lm.pre[[i]], data.frame(log_Q = log(l.Q[[i]]$mean_Q))))),
#                                                                              l.BCF.post[[i]]*exp(as.numeric(predict(l.lm.post[[i]], data.frame(log_Q = log(l.Q[[i]]$mean_Q)))))
# )))%>%
#   purrr::set_names(names(l.Q))
# 
# l.Yield<-lapply(names(l.lm), \(i) l.pred.C[[i]]$pred.C*l.Q.v[[i]]*(1/l.DA[[i]])*unit_conversion)%>%purrr::set_names(names(l.Q))
# 
# df.Yield <- lapply(l.Yield, sum)%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'OLS.AANY.2s.medC') # units of kg/ha/year
# 
# df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')
# 
# # check: 
# 
# # make sure list orders are the same:
# 
# names(l.lm.pre)==names(l.lm.post)
# names(l.lm.pre)==names(l.RP)
# names(l.lm.pre)==df.median.Q$site_no
# 
# # use models in ifelse to predict concentration column:
# 
# l.test<-lapply(seq_along(l.lm.pre), \(i) l.RP[[i]]%>%mutate(Seg_C.medC = ifelse(C < df.median.C$median_C[i], 
#                                                                               predict(l.lm.pre[[i]], data.frame(log_Q = l.RP[[i]]$log_Q)),
#                                                                               predict(l.lm.post[[i]], data.frame(log_Q = l.RP[[i]]$log_Q))
# ))%>%
#   mutate(Seg_type = ifelse(C < df.median.C$median_C[i],
#                            'Pre_C_med',
#                            'Post_C_med'))
# )%>%purrr::set_names(names(l.RP))
# 
# # check with plot:
# 
# i <- 2
# 
# ggplot(l.test[[i]], aes(x = log_Q, y = log_C))+
#   geom_point()+
#   geom_smooth(method = 'lm')+
#   geom_line(aes(x = log_Q, y = Seg_C.medC, color = Seg_type))
# 
# df.OLS[i,10:14]
# 
# # I want to look at facet plot of all these:
# 
# # x <- bind_rows(l.test, .id = "Name")
# 
# # ggplot(x, aes(x = log_Q, y = log_C))+
# #   geom_point()+
# #   geom_smooth(method = 'lm')+
# #   geom_line(aes(x = log_Q, y = Seg_C.medC, color = Seg_type))+
# #   facet_wrap('Name', scales = 'free') +
# #   theme(
# #     strip.background = element_blank(),
# #     strip.text.x = element_blank()
# #   )
# 
# # I just realized that separation by C doesnt scale well in this workflow

#### ~ 1.4) baseflow separation ####

library(grwat)

# example site in workflow:

# Calculate baseflow:

# i<-2

# hdata = df.NWIS.Q %>%
#   filter(site_no == df.OLS$Name[i])%>%
#   mutate(Qbase = gr_baseflow(X_00060_00003, method = 'lynehollick', a = .925, passes = 3),
#          Date=as.Date(Date))%>%
#   mutate(Q_type=ifelse(X_00060_00003<=Qbase, 'Baseflow', 'Stormflow'))

# Visualize:

# ggplot(hdata) +
#   geom_area(aes(Date, X_00060_00003), fill = 'steelblue', color = 'black') +
#   geom_area(aes(Date, Qbase), fill = 'orangered', color = 'black')+
#   geom_point(aes(Date, X_00060_00003, color=Q_type))+
#   scale_x_date(limits = c(ymd(20070601), ymd(20070731)))

# pretty sick

# now scale up:

df.Q.w.base <- df.Q%>%
  split(., f=.$site_no)%>%
  lapply(., \(i) i%>%mutate(Q_base = gr_baseflow(X_00060_00003, method = 'lynehollick', a = 0.925, passes = 3)))%>%
  lapply(., \(i) i%>%mutate(Q_type=ifelse(X_00060_00003<=Q_base, 'Baseflow', 'Stormflow')))%>%
  bind_rows(.)

# join this to df.RP:

df.RP <- left_join(df.RP, df.Q.w.base%>%select(site_no,Date,Date_AAH,Q_type), by = c('Name'='site_no', 'sample_dt'='Date'))

# now add OLS columns for Int and slope of base and Storm flow:

# create two df.OLS, one for base and storm:

temp.base<-df.RP%>%
  group_by(Name)%>%
  filter(Q_type == 'Baseflow')%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.Int.2s_hydrosep_pre = OLS.co[1], 
            OLS.Slope.2s_hydrosep_pre = OLS.co[2])
  }) %>%
  ungroup

temp.storm<-df.RP%>%
  group_by(Name)%>%
  filter(Q_type == 'Stormflow')%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.Int.2s_hydrosep_post = OLS.co[1], 
            OLS.Slope.2s_hydrosep_post = OLS.co[2])
  }) %>%
  ungroup

# combine base and storm:

temp <- left_join(temp.base, temp.storm,by='Name')

# merge with df.OLS:

df.OLS <- left_join(df.OLS,temp,by='Name')%>%as.data.frame()

# AANY:

l.lm.pre <- df.RP%>%
  filter(Q_type == 'Baseflow')%>%
  split(., f=.$Name)%>%
  lapply(., \(i) lm(log_C ~ log_Q, i))

l.lm.post <- df.RP%>%
  filter(Q_type == 'Stormflow')%>%
  split(., f=.$Name)%>%
  lapply(., \(i) lm(log_C ~ log_Q, i))

# make sure list orders are the same:

names(l.lm.pre)==names(l.lm.post)
names(l.lm.pre)==names(l.RP)
names(l.lm.pre)==df.median.Q$site_no

# back to workflow:

l.BCF.pre<-lapply(l.lm.pre, \(i) exp((summary(i)$sigma^2)/2))
l.BCF.post<-lapply(l.lm.post, \(i) exp((summary(i)$sigma^2)/2))

# first approach:

l.pred.C<-lapply(seq_along(l.lm.pre), \(i) l.Q[[i]]%>%mutate(pred.C = ifelse(mean_Q < df.median.Q$median_Q[i],
                                                                             l.BCF.pre[[i]]*exp(as.numeric(predict(l.lm.pre[[i]], data.frame(log_Q = log(l.Q[[i]]$mean_Q))))),
                                                                             l.BCF.post[[i]]*exp(as.numeric(predict(l.lm.post[[i]], data.frame(log_Q = log(l.Q[[i]]$mean_Q)))))
)))%>%
  purrr::set_names(names(l.Q))

l.Yield<-lapply(names(l.lm), \(i) l.pred.C[[i]]$pred.C*l.Q.v[[i]]*(1/l.DA[[i]])*unit_conversion)%>%purrr::set_names(names(l.Q))

df.Yield <- lapply(l.Yield, sum)%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'OLS.AANY.2s.hydrosep') # units of kg/ha/year

df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')

# second approach:

l.Q.2.w.base<-split(df.Q.w.base, f=df.Q.w.base$site_no)

l.pred.C<-lapply(seq_along(l.lm.pre), \(i) l.Q.2.w.base[[i]]%>%mutate(pred.C = ifelse(Q_type == 'Baseflow',
                                                                               l.BCF.pre[[i]]*exp(as.numeric(predict(l.lm.pre[[i]], data.frame(log_Q = log(l.Q.2.w.base[[i]]$X_00060_00003))))),
                                                                               l.BCF.post[[i]]*exp(as.numeric(predict(l.lm.post[[i]], data.frame(log_Q = log(l.Q.2.w.base[[i]]$X_00060_00003)))))
)))%>%
  purrr::set_names(names(l.Q.2.w.base))

l.Yield<-lapply(names(l.lm), \(i) data.frame(Date_AAH =l.pred.C[[i]]$Date_AAH, pred.yield = l.pred.C[[i]]$pred.C*l.Q.2.w.base[[i]]$X_00060_00003*(1/l.DA[[i]])*unit_conversion))%>%purrr::set_names(names(l.Q.2.w.base))
l.Yield.AAH <- lapply(l.Yield, \(i) i%>%group_by(Date_AAH)%>%summarise(mean_daily_yield=mean(pred.yield)))
df.Yield <- lapply(l.Yield.AAH, \(i) sum(i$mean_daily_yield, na.rm = T))%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'OLS.AANY.2s.hydrosep.method2') # units of kg/ha/year

df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')

# check with plot:

# plots:

i<-5

ggplot(df.RP%>%filter(Name == df.OLS$Name[i]), aes(x = log_Q, y = log_C,color = Q_type))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap('n_sample_rank', scales = 'free') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

df.OLS[i,12:17]

#### ~ 1.5) advanced baseflow separation ####

# add precip and temperature to df.NWIS.Q:

# use function to download climate data for each sites watershed boundary
# and summarize it into a single mean daily value for each flow observation...
# that is a lot of data!

# i am going to wait on this since I dont know if it is worth it

#### ~ 2) non-CQ metrics: CV_C/CV_Q, median and mean concentrations ####

df.CV <- df.RP%>%
  group_by(Name)%>%
  summarise(CV_C =  cv(C), CV_Q = cv(Q), CV_CQ =cv(C)/cv(Q))

df.medC <- df.RP%>%
  group_by(Name)%>%
  summarise(medC = median(C), meanC = mean(C))

#### ~ 3) WRTDS ####

# library(EGRET)

# write flow dfs to csv files and read back in:

# x<-df.NWIS.Q%>%
#   filter(site_no %in% df.RP$Name)%>%
#   split(., f=.$site_no)
# 
# lapply(seq_along(x),\(i) x[[i]]%>%
#            select(Date,X_00060_00003)%>%
#            mutate(Date = as.Date(Date))%>%
#            rename(Q=X_00060_00003)%>%
#            write.csv(., file = paste0('Raw_Data/Q_csv/', names(x)[i], '.csv'), row.names = F)
#          )
# 
# l.readUserDaily<-list()
# 
# for(i in seq_along(l.Q)){
#   
#   l.readUserDaily[[i]]<-readUserDaily('Raw_Data/Q_csv', paste0(names(l.Q)[i], '.csv'))
#   
# }
# 
# names(l.readUserDaily) <- names(l.Q)
# 
# save(l.readUserDaily, file = "Raw_Data/Q_csv/l.readUserDaily.Rdata")

# load("Raw_Data/Q_csv/l.readUserDaily.Rdata")

# write sample dfs to csv files and read back in:

# x<-df.RP%>%
#   ungroup()%>%
#   split(., f=.$Name)
# 
# lapply(seq_along(x),\(i) x[[i]]%>%
#          rename(Date=sample_dt)%>%
#            select(Date,C)%>%
#            mutate(Date = as.Date(Date))%>%
#            mutate(Remarks = NA, .after = 1)%>%
#            write.csv(., file = paste0('Raw_Data/C_csv/', names(x)[i], '.csv'), row.names = F))
# 
# l.readUserSample<-list()
# 
# for(i in seq_along(l.Q)){
# 
#   l.readUserSample[[i]]<-readUserSample('Raw_Data/C_csv', paste0(names(l.Q)[i], '.csv'))
# 
# }
# 
# names(l.readUserSample) <- names(l.Q)
# 
# save(l.readUserSample, file = "Raw_Data/Q_csv/l.readUserSample.Rdata")

# load("Raw_Data/Q_csv/l.readUserSample.Rdata")

# download info df and save:

# l.readNWISInfo <- lapply(names(l.Q), \(i) readNWISInfo(i,  pcode, interactive = FALSE))
# 
# save(l.readNWISInfo , file = "Raw_Data/Q_csv/l.readNWISInfo.Rdata")

# load("Raw_Data/Q_csv/l.readNWISInfo.Rdata")

# build l.eList:

# l.eList <- lapply(seq_along(l.Q), \(i) mergeReport(l.readNWISInfo[[i]],l.readUserDaily[[i]],l.readUserSample[[i]]))

# run WRTDS using modified function:

# l.eList <- lapply(l.eList, fun.R.modelEstimation)

# save l.elist with WRTDS results:

# save(l.eList , file = "Processed_Data/l.eList.Rdata")

# load in l.elist:

# load("Processed_Data/l.eList.Rdata")

# look at annual results:

# ApMayJuneResults <- setupYears(l.eList[[11]]$Daily, paLong = 12, paStart = 1)

# plotFluxHist(l.eList[[11]])

# x<-l.eList[[11]]$INFO

#### ~ 4) Simple AANY ####

# calculate yield for each CQ observation (i.e. day):

l.AANY.simple <- lapply(names(l.RP), \(i) l.RP[[i]]$C*l.RP[[i]]$Q*(1/l.DA[[i]])*unit_conversion)%>%purrr::set_names(names(l.RP))

# take the average over all observations to get mean daily yield:

l.AANY.simple <- lapply(l.AANY.simple, mean)

# then multiply by 365:

l.AANY.simple <- lapply(l.AANY.simple, \(i) i*365)

# convert to df:

df.Yield <- bind_rows(l.AANY.simple) %>% pivot_longer(cols = everything(), names_to = 'Name', values_to = 'AANY.simple')

# merge with df.OLS:

df.OLS<-left_join(df.OLS, df.Yield, by = 'Name')
 
#### ~ Putting response variables together into single df: ####

df.Response <- left_join(df.OLS, df.CV, by = 'Name')%>%
  left_join(.,df.medC, by = 'Name')

# look at distribution of response variables:

# set up plotting df:

df.plot <- df.Response%>%select(-c(Name, CV_C, CV_Q))%>%
  pivot_longer(cols = everything())%>%
  mutate(group = case_when(grepl("AANY", name) ~ "AANY",
                           grepl("Int", name) & grepl("pre", name) ~ "Int_pre",
                           grepl("Int", name) & grepl("1s", name) ~ "Int_post",
                           grepl("Int", name) & grepl("post", name) ~ "Int_post",
                           grepl("Slope", name) & grepl("pre", name) ~ "Slope_pre",
                           grepl("Slope", name) & grepl("1s", name) ~ "Slope_post",
                           grepl("Slope", name) & grepl("post", name) ~ "Slope_post"))%>%
  mutate(group = ifelse(is.na(group), 'Other', group))%>%
  mutate(name2 = case_when(grepl("1s", name) ~ "1s",
                          grepl("hydrosep", name) ~ "hydrosep",
                          grepl("medC", name) ~ "medC",
                          grepl("medQ", name) ~ "medQ",
                          grepl("CV_CQ", name) ~ "CV_CQ",
                          grepl("meanC", name) ~ "meanC",
                          grepl("simple", name) ~ "simple"
                          ))%>%
  mutate(name2 = ifelse(grepl('method2',name), paste(name2, 'method2'), name2))

# set group to factor with specific order:

df.plot$group <- factor(df.plot$group, levels = c("AANY", "Other", "Int_pre", "Int_post", "Slope_pre", "Slope_post"))

# make plot:

# density plot:

# ggplot(df.plot, aes(x=value, color = name2, fill = name2)) +
#   geom_density(alpha=0.5, position="identity")+
#   facet_wrap(~group, ncol=2, scales = 'free')

# violin plot:

ggplot(df.plot, aes(x=name2, y=value))+
  geom_violin()+
  theme(axis.title.x=element_blank())+
  geom_hline(yintercept=0,color='red')+
  facet_wrap(~group, ncol=2, scales='free')

# alternatively, can use gridextra package if you wanted individual legends
# for each facet: https://stackoverflow.com/questions/14840542/place-a-legend-for-each-facet-wrap-grid-in-ggplot2

# library(gridExtra)
# 
# # set up list (to be each facet using %+%):
# 
# l.df.plot <- split(df.plot,f = df.plot$group)
# 
# # density plot:
# 
# # p1 <- ggplot(l.df.plot[[1]],aes(x=value, color = name2, fill = name2)) + 
# #   geom_density(alpha=0.5, position="identity")+
# #   facet_wrap(~group, ncol=3)
# 
# # violin plot:
# 
# p1<-ggplot(l.df.plot[[1]], aes(x=name2, y=value))+
#   geom_violin()+
#   theme(axis.title.x=element_blank())+
#   geom_hline(yintercept=0,color='red')+
#   facet_wrap(~group, ncol=3)
#   
# # then create plots for each facet:
# 
# p2 <- p1 %+% l.df.plot[[2]]
# p3 <- p1 %+% l.df.plot[[3]]
# p4 <- p1 %+% l.df.plot[[4]]
# p5 <- p1 %+% l.df.plot[[5]]
# p6 <- p1 %+% l.df.plot[[6]]
# 
# # and combine:
# 
# grid.arrange(p1,p2,p3,p4,p5,p6)

#






























#### Set up Predictors ####

# first pass predictors are from G2:

# reduce the predictor set down to those we believe are defensible:

# temp1<-read.csv('Processed_Data/G2_pred_to_keep_2.csv')[,1] # I went in to the variable description excel sheet and isolated these variables as well:

# pred_to_keep<-c("STAID",
#                 "HYDRO_DISTURB_INDX",
#                 "pre1990_STOR",
#                 "FRAGUN_BASIN", # names(l.G2$Landscape_Pat),
#                 "DEVNLCD06","FORESTNLCD06","PLANTNLCD06", # names(l.G2$LC06_Basin),
#                 "CDL_CORN", # names(l.G2$LC_Crops),
#                 "PDEN_2000_BLOCK", "RIP800_DEV", "RIP800_FOREST", "RIP800_PLANT", # temp1,
#                 'HGA', 'HGB', 'HGC', 'HGD',
#                 "ELEV_MEDIAN_M_BASIN","ELEV_STD_M_BASIN","RRMEDIAN")

# df.G2.reduced<-df.G2%>%select(c('STAID', pred_to_keep))%>%filter(STAID %in% df.TP_CQ$site_no) # filtering df.G2 to this predictor list and sites:

# df.G2.reduced<-df.G2.reduced%>%
#   rowwise() %>%
#   mutate(MAIN_RIP_DEV_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("DEV")), na.rm = TRUE),
#          MAIN_RIP_FOREST_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("FOREST")), na.rm = TRUE),
#          MAIN_RIP_PLANT_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("PLANT")), na.rm = TRUE),
#          .keep = 'unused'
#          ) # take the average of the RIP and MAIN for each land use:

# df.G2.reduced<-df.G2.reduced%>%mutate(sum_of_misc_crops = sum(c_across(ends_with(c("COTTON", "RICE", "SORGHUM", "SUNFLOWERS", "PEANUTS", "BARLEY", "DURUM_WHEAT", "OATS", "DRY_BEANS", "POTATOES", "ORANGES", 'OTHER_CROPS', 'IDLE'))), na.rm = TRUE), .keep = 'unused') # adding up the other crops:

# names(df.G2.reduced) # check

# save(df.G2.reduced, file = 'Processed_Data/df.G2.reduced.Rdata') # export this predictor set for the 42 sites:

load('Processed_Data/df.G2.reduced.Rdata')

# second pass predictors are from datalayers:

load('Processed_Data/df.datalayers.62.Rdata') # read in datalayers predictors:

df.datalayers<-df.datalayers.62 # rename to just df.datalayers to work with correlations workflow:

# convert CSA columns to 0-1 %:

df.datalayers <- df.datalayers%>%
  mutate(CSA_perc=CSA_perc/100,
         Ag.CSA_perc=Ag.CSA_perc/100,
         RIP.CSA.100=RIP.CSA.100/100,
         RIP.CSA.800=RIP.CSA.800/100)

df.datalayers<-df.datalayers%>%filter(!is.na(CSA_perc)) # *note* that one site has NA for pedictors that require soils since it had no sSurgo coverage. I will remove that site:

# # *note* df.datalayers from second pass does not have column "Dbl_Crop_WinWht/Soybeans"
# # *note* df.datalayers has 61 sites now. 
# # It did have 62 sites while df.TP_CQ (for second pass) has 63 sites
# # this is not a big deal for the correlaitons because it will just omit the
# # row with NAs in all the predictor columns, but this row will need to be omitted
# # when running the MLR
# # as per a first run of this workflow, remove emergentwetlands, soybeans, winter wheat, grassland
# # and replace NA with zero (this is ok since I removed the site with NA for soils where a zero would be wrong):

df.datalayers<-df.datalayers%>%
  select(-c(R_EMERGWETNLCD06, R_SNOWICENLCD06, R_SHRUBNLCD06, R_BARRENNLCD06,R_GRASSNLCD06, Soybeans, Winter_Wheat, Spring_Wheat, `Dbl_Crop_WinWht/Soybeans`, R_WOODYWETNLCD06))%>%
  replace(is.na(.), 0)

#












#### Correlations and MLR ####

# merge response and predictors:

# first pass:

df.setup<-left_join(df.Response, df.G2.reduced, by = c('Name'='STAID')) # merge predictor variables with response variable df:

# second pass:

df.setup<-left_join(df.Response, df.datalayers, by = 'Name')%>%
  drop_na(R_CROPSNLCD06)# merge predictor variables with response variable df:

# remove CV_C and CV_Q,and n

df.setup <- select(df.setup, -c(CV_C, CV_Q, n)) 

# remove method 1 AANY:

df.setup <- select(df.setup, -c(matches("AANY")&!contains(c("method2", 'simple'))))

# Remove water, HSG.coverage, forest types:

df.setup <- select(df.setup, -c(R_WATERNLCD06, HSG.coverage, R_MIXEDFORNLCD06, R_EVERGRNLCD06, R_DECIDNLCD06))

# save(df.setup, file='Processed_Data/df.setup.Rdata')

# set up seq of column numbers of response and predictors in df.setup:

names(df.setup)

v.resp.cols <- 2:18

v.pred.cols <- 19:ncol(df.setup)

# look at correlation matrix of response and predictors:

# response :

corr.resp <- round(cor(df.setup[v.resp.cols], method = 's'), 1)

ggcorrplot(corr.resp, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 2, 
           method="square", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Correlogram of Response Variables", 
           ggtheme=theme_bw)+
  theme(axis.text.x=element_text(angle=40,hjust=1, size = 7), axis.title.x=element_blank())

# predictors:

corr.pred <- round(cor(df.setup[v.pred.cols], method = 's'), 1)

ggcorrplot(corr.pred, 
           hc.order = T, 
           hc.method = "mcquitty",
           type = "lower", 
           lab = TRUE, 
           lab_size = 2,
           method="square", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Correlogram of Predictor Variables", 
           ggtheme=theme_bw)+
  theme(axis.text.x=element_text(angle=40,hjust=1, size = 7), axis.title.x=element_blank())

# look at univariate plots among response variables: to do this:

# # create df of response vars with strong correlations:
# 
# df.corr.resp<-as.data.frame(corr.resp)%>% 
#   tibble::rownames_to_column(., "Resp")%>%
#   pivot_longer(., cols = -Resp, names_to = 'Resp.2', values_to = 'Cor')%>%
#   filter(., abs(Cor)>0.7)%>%filter(., Resp != Resp.2)%>%
#   select(1,2)
# 
# df.corr.pred<-as.data.frame(corr.pred)%>% 
#   tibble::rownames_to_column(., "Pred")%>%
#   pivot_longer(., cols = -Pred, names_to = 'Pred.2', values_to = 'Cor')%>%
#   filter(., abs(Cor)>0.7)%>%filter(., Pred != Pred.2)%>%
#   select(1,2)
# 
# # reduce this df down to just unique pairs, then group by and summarize column of Resp.2 for each Resp 1:
# 
# df.corr.resp <- unique(t(apply(as.matrix(df.corr.resp), 1, sort)))%>%
#   as.data.frame(.)%>%
#   rename(Resp=1,Resp.2=2)%>%
#   arrange(Resp)%>%
#   group_by(Resp)%>%
#   summarise(Resp.2=list(Resp.2))
# 
# df.corr.pred <- unique(t(apply(as.matrix(df.corr.pred), 1, sort)))%>%
#   as.data.frame(.)%>%
#   rename(Pred=1,Pred.2=2)%>%
#   arrange(Pred)%>%
#   group_by(Pred)%>%
#   summarise(Pred.2=list(Pred.2))
# 
# # loop through resp. vars. and create plotting df:
# 
# df.plot.resp<-data.frame(Resp.1.name=NA,Resp.1.value=NA,Resp.2.name=NA,Resp.2.value=NA)
# 
# for (i in seq_along(df.corr.resp$Resp)){
#   
#   df <- df.setup%>%
#     select(df.corr.resp$Resp[i], any_of(df.corr.resp$Resp.2[i][[1]]))%>%
#     pivot_longer(cols = -c(df.corr.resp$Resp[i]), names_to = 'Resp.2.name', values_to = 'Resp.2.value')%>%
#     pivot_longer(cols = c(df.corr.resp$Resp[i]), names_to = 'Resp.1.name', values_to = 'Resp.1.value')%>%
#     select(c(3,4,1,2))
#   
#   df.plot.resp <- bind_rows(df.plot.resp, df)
#   
# }
# 
# df.plot.resp <- df.plot.resp[-1,]
# 
# df.plot.pred<-data.frame(Pred.1.name=NA,Pred.1.value=NA,Pred.2.name=NA,Pred.2.value=NA)
# 
# for (i in seq_along(df.corr.pred$Pred)){
#   
#   df <- df.setup%>%
#     select(df.corr.pred$Pred[i], any_of(df.corr.pred$Pred.2[i][[1]]))%>%
#     pivot_longer(cols = -c(df.corr.pred$Pred[i]), names_to = 'Pred.2.name', values_to = 'Pred.2.value')%>%
#     pivot_longer(cols = c(df.corr.pred$Pred[i]), names_to = 'Pred.1.name', values_to = 'Pred.1.value')%>%
#     select(c(3,4,1,2))
#   
#   df.plot.pred <- bind_rows(df.plot.pred, df)
#   
# }
# 
# df.plot.pred <- df.plot.pred[-1,]
# 
# # make plot:
# 
# ggplot(df.plot.resp,aes(x=Resp.1.value, y=Resp.2.value, color = Resp.2.name))+
#   geom_point()+
#   geom_smooth(method = 'lm')+
#   facet_wrap('Resp.1.name', scales = 'free')
# 
# ggplot(df.plot.pred,aes(x=Pred.1.value, y=Pred.2.value, color = Pred.2.name))+
#   geom_point()+
#   geom_smooth(method = 'lm')+
#   facet_wrap('Pred.1.name', scales = 'free')

#

# now run correlations between response and predictor variables. to do this: (I orginally did this workflow using n_months (C:\PhD\Research\Mohawk\Code\Mohawk_Regression-analyizing_predictor_df.R)

# set up variable for number of sites:

n_sites<-dim(df.setup)[1]

# look at names of df.setup to see which columns are respose variables:

names(df.setup)

# use the corrr package to correlate() and focus() on your variable of choice.

df.cor <- df.setup %>%
  corrr::correlate(method = 'spearman') %>%
  corrr::focus(2:18)%>%
  pivot_longer(cols= -term, names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# I also want to try standardizing the predictors:

# df.cor.0to1<-df.OLS_Sens %>%
#   mutate(across(7:last_col(), ~ (.-min(.))/(max(.)-min(.))))%>% # could add term: .names = "{.col}_standarized")
#   correlate(method = 'spearman') %>%
#   focus(c(OLS.I, OLS.S, Sen.I, Sen.S, Yield))%>%
#   pivot_longer(cols= c(2:6), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
#   mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
#   mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
#   drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero
# 
# df.cor$Spearman_Correlation==df.cor.0to1$Spearman_Correlation # looking at the regualr and 0 to 1 standarized spearman correlaitons, they are the same
# so not going to proceed with 0 to 1

# then plot results: to do this:

# filter df.cor to just the signficant correlates:

df.cor<-filter(df.cor, sig_0.05=='sig')

# extract just the top 10 correlates for each response variable:

df.cor.top10 <- df.cor%>%group_by(CQ_Parameter)%>%arrange(desc(abs(Spearman_Correlation)))%>%slice_max(abs(Spearman_Correlation), n=3) %>% ungroup()
  
# see how many unique predictors are signficant:

unique(df.cor.top10$term)

# make plot:

ggplot(df.cor.top10, aes(x=term, y=Spearman_Correlation, color=term))+
  geom_bar(stat = "identity")+
  # theme(axis.title.y=element_blank())+
  facet_wrap(~CQ_Parameter, ncol = 4,scales="free")

#### Model Building ####

# set up list for each response variable:

l.setup <- lapply(v.resp.cols, \(i) df.setup[,c(i, v.pred.cols)] %>% rename(term = 1)) %>% purrr::set_names(names(df.setup[v.resp.cols]))

# set up list for just AANY, med/mean C, and CV C/Q response variables:

l.setup.AANY <- l.setup[grep("AANY|C", names(l.setup))] 

#### ~ MLR ####

#### ~| VIF proof of multicolineatiry ####

# run full model and calculate vif:

l.lm.full <- lapply(l.setup, \(i) lm(term ~ ., data = i))

l.vif.full <- lapply(l.lm.full, car::vif)

# doesnt run: 'Error: there are aliased coefficients in the model'

#### ~| real v.s. log space models ####

# wont work, zeros in data

#### ~| Remove all but one highly Correlated ####

# find correlation:

highCorr <- findCorrelation(corr.pred, .75, names = T) # find highly correlated columns:

# loop through these variables and save kable table:

leave <- 1

for (leave in seq_along(highCorr)){
  
  # look at correlaogram again:
  
  # x <- select(df.setup, -highCorr[-leave]) # remove all but one of these from predictor list
  # corr.pred <- round(cor(x[19:ncol(x)], method = 's'), 1) # rerun correlation matrix:
  # ggcorrplot(corr.pred,
  #            hc.order = T,
  #            hc.method = "mcquitty",
  #            type = "lower",
  #            lab = TRUE,
  #            lab_size = 2,
  #            method="square",
  #            colors = c("tomato2", "white", "springgreen3"),
  #            title="Correlogram of Predictor Variables",
  #            ggtheme=theme_bw)+
  #   theme(axis.text.x=element_text(angle=40,hjust=1, size = 7), axis.title.x=element_blank())
  
  # looks good
  
  # remove all but one of the highly correlated predictors from l.setup:
  
  l.setup.MLR <- lapply(l.setup.AANY, \(i) select(i, -highCorr[-leave]))
  
  # look at example correlation between resp. and predictors:
  
  # x <- l.setup.MLR[[1]] %>% pivot_longer(cols=-term)
  # 
  # ggplot(x, aes(x=value, y=term))+
  #   geom_point()+
  #   geom_smooth(method='lm')+
  #   facet_wrap(~name,scales='free')
  
  # build and compare models using full lm, regsubsets, and stepAIC:
  
  df.m.compare <- lapply(l.setup.MLR, fun.compare.linear.model.selections) %>% 
    purrr::set_names(names(l.setup.MLR)) %>% 
    bind_rows(., .id='Resp') 
  
  # look at just the top models for each response variable:
  
  df.table <- df.m.compare %>% 
    group_by(Resp) %>% 
    slice_max(order_by = adjR2, n = 1) %>% 
    slice(1) %>% 
    mutate(Resp = paste(Resp, '-', round(adjR2,2))) %>% 
    mutate(Resp = gsub("OLS.", "", Resp)) %>% 
    arrange(desc(adjR2)) %>% 
    select(-c(Type, adjR2)) %>% 
    select_if(~sum(!is.na(.)) > 0) %>% 
    rename('Resp - AdjR2'=1)
  
  df.table %>%
    kbl(caption = "Table ?") %>%
    kable_classic(full_width = F, html_font = "Times") %>% 
    row_spec(seq(1,nrow(df.table),2), background="#FF000020") %>% 
    column_spec(which(names(df.table)==highCorr[leave]), background = 'yellow')
  
}


# subsitute other highly correlated predictors for the one that was 
# kept (i.e. highCorr[1]):

highCorr[leave]

#






# Train the models:

l.step.model <- lapply(l.setup.MLR, \(i) train(term ~., data = i,
                                method = "leapForward", 
                                tuneGrid = data.frame(nvmax = 1:5),
                                trControl = trainControl(method = "cv", number = 10)
))

# build MLR lm using function to build lm from terms in step.model$finalmodel final model :

l.m.final <- lapply(names(l.step.model), \(i) fun.train.MLR.to.lm(l.step.model[[i]], l.setup.MLR[[i]])) %>% purrr::set_names(names(l.step.model))

# extract significant model terms:

l.terms.sig <- lapply(l.m.final, fun.extract.sig.model.terms)

# extract model R2:

l.R2 <- lapply(l.m.final, \(i) summary(i)$adj.r.squared)

# extract variable importance:

l.varImp <- lapply(l.step.model, \(i) varImp(i, scale = FALSE)$importance %>% arrange(desc(Overall)) %>% tibble::rownames_to_column(., "Predictor"))

# make table of results:

# for each response variable (row), make columns of model significant terms, AdjR2, and top three most important variables, and variable importance values:

# look at model results:

i<-1

step.model <- l.m[[i]]
df <- l.setup.MLR[[i]]

step.model$results
step.model$bestTune
summary(step.model$finalModel)
coef(step.model$finalModel, as.numeric(step.model$bestTune))
lmImp <- varImp(step.model, scale = FALSE)$importance %>% arrange(desc(Overall))



pred<-names(coef(step.model$finalModel, as.numeric(step.model$bestTune)))[-1]



pred<-gsub("`", "", pred)

# select down df:

df.2<-df%>%select(term, pred)

# make lm:

m<-lm(term~., data=df.2)

# m.list[[i]]<-m

summary(m)

# look at univarite plots:

p<-df.2%>%pivot_longer(cols = 2:last_col(), names_to = 'Predictor', values_to = 'Value')%>%
  rename(!!names(l.cor.MLR.full)[i]:=term)%>%
  ggplot(., aes(x=Value, y=!!sym(names(l.cor.MLR.full)[i])))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap('Predictor', scales = 'free')

# look at model residuals:

# plot(m)




#### ~ PLS ####

# split l.setup into two for training and testing:

training.samples <- l.setup[[1]]$term %>%
  createDataPartition(p = 0.8, list = FALSE)

l.train <- lapply(l.setup, \(i) i[training.samples, ])
l.test <- lapply(l.setup, \(i) i[-training.samples, ])


# Build the model on training set:

l.model <- lapply(l.train, \(i) train(term~., data = i, 
                                              method = "pls",
                                              scale = TRUE,
                                              trControl = trainControl("cv", number = 10),
                                              tuneLength = 10)
)

i<-1

# Plot model RMSE vs different values of components

plot(l.model[[i]])

# Print the best tuning parameter ncomp that minimize the cross-validation error, RMSE:

l.model[[i]]$bestTune

# Summarize the final model:

summary(l.model[[i]]$finalModel)

# Make predictions:

l.predictions <- lapply(seq_along(l.model), \(i) l.model[[i]] %>% predict(l.test[[i]]))

# Model performance metrics:

data.frame(
  RMSE = caret::RMSE(l.predictions[[i]], l.test[[i]]$term),
  Rsquare = caret::R2(l.predictions[[i]], l.test[[i]]$term)
)


#















# create a list of each CQ parameter and format it for ggplotting:

l.cor<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>%
           arrange(Spearman_Correlation)%>%
           slice_head(n = 10)%>%
           bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           filter(!between(Spearman_Correlation, -0.25,.25))) # Order by correlation strength

x <- l.cor[[1]]

# make plot list using lapply:

plist<-lapply(l.cor, \(i) i%>%ggplot(aes(x = term, y = Spearman_Correlation, color = sig_0.05)) +
                geom_bar(stat = "identity") +
                scale_color_manual(values = c("not" = "red", "sig" = "blue"),na.value = NA, drop = FALSE)+
                facet_wrap('CQ_Parameter')+
                ylab(paste('Spearman Correlation')) +
                xlab("Watershed Attribute")+
                theme(axis.text.x=element_text(angle=40,hjust=1, size = 7), axis.title.x=element_blank())+
                theme(legend.position="bottom"))

plist[[15]]

# arrange plots:

p1<-ggpubr::ggarrange(plotlist = plist, ncol=5, nrow=4, common.legend = TRUE, legend="bottom")

# add title:

p1<-annotate_figure(p1, top = text_grob(paste(ncode, "Correlated Against Gauges 2"),
                                        color = "red", face = "bold", size = 14))
# plot:

p1

# save(p1, file = 'Processed_Data/p1.Rdata')

# univariate plots of the top correlates/small land use types:

# Lets do plots of OLS intercept: to do this:

# determine the top correlates (the numberof the list determined which CQ parameter)
# here 1 = OLS intercept:

OLS<-l.cor[[13]]%>%arrange(desc(Spearman_Correlation))

# make univariate plots (facets) of Y (determined in [[]]) above and predictors:
# note the predictor column isturned into an ordered factor based on thespearman correlation values
# to makethe facet plots in order:

x <- df.setup%>%
  select(Name, OLS$CQ_Parameter[1], OLS$term)%>%
  pivot_longer(cols = c(3:last_col()), names_to = 'Type', values_to = 'Value')%>%
  drop_na(Value)%>%
  # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
  left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
  mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
  ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Type', scales = 'free')+
  ggtitle(paste('Univariate plots of', OLS$CQ_Parameter[1], 'against Top Gauges 2 Correlates for', OLS$CQ_Parameter[1], 'ordered from highest to lowest Spearman Rank Correlation' ))

# create a function to use lapply to plot all 4 univariate sets:

# number<-5

# make_plot<-function(number){
#
#   OLS<-l.cor[[number]]%>%arrange(desc(Spearman_Correlation))
#
#   x<-df.OLS_Sens%>%
#     select(Name, OLS$CQ_Parameter[1], OLS$term)%>%
#     pivot_longer(cols = c(3:last_col()), names_to = 'Type', values_to = 'Value')%>%
#     drop_na(Value)%>%
#     # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
#     left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
#     mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
#     ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
#     geom_smooth(method = 'lm')+
#     geom_point()+
#     facet_wrap('Type', scales = 'free')+
#     ggtitle(paste('Univariate plots of', OLS$CQ_Parameter[1], 'against Top Gauges 2 Correlates for', OLS$CQ_Parameter[1], 'ordered from highest to lowest Spearman Rank Correlation' ))
#
#   x
#
# }

# use function in lapply (clear plot list first):

# lapply(c(1:5), \(i) make_plot(i))

#

# end of comment out



























#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Correlations with 'top50' ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# note: only run for first pass

# start comment out:

# # determine the median flow rate for each sites CQ observations: to do this:
# 
# df.TP_CQ.top50<-split(df.TP_CQ, f = df.TP_CQ$site_no)%>% # split the df into list of dfs for each site, 
#   lapply(., \(i) i%>%filter(X_00060_00003>median(i$X_00060_00003)))%>% # filter to the median flow rate for each site,
#   bind_rows(.) # bind back into single dataframe:
# 
# # now run through the Correlations workflow for the top50 OLS model:
# 
# df.Correlation.top50<-df.TP_CQ.top50%>%
#   rename(Name = site_no)%>%
#   select(Name, sample_dt,result_va, X_00060_00003)%>%mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
#   filter(is.finite(log_C))%>%
#   filter(is.finite(log_Q))%>%
#   left_join(., df.G2.reduced, by = c('Name'='STAID'))
# df.OLS.top50<-df.Correlation.top50%>% # create a dataframe with OLS and sens slope intercept and slopes:
#   group_by(Name)%>%
#   do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
#   summarize(., OLS.I = OLS.co[1], 
#             OLS.S = OLS.co[2])
#   }) %>%
#   ungroup
# df.Sens.top50<-df.Correlation.top50%>%
#   group_by(Name)%>%
#   do({ Sens.co<-zyp.sen(log_C~log_Q,.)
#   summarize(., Sen.I = Sens.co$coefficients[[1]],
#             Sen.S= Sens.co$coefficients[[2]])
#   }) %>%
#   ungroup
# df.OLS_Sens.top50<-left_join(df.OLS.top50, df.Sens.top50, by = 'Name') # merge OLS and Sens:
# # here is where top 50 workflow differes: need to use the recalcualted df.Yield using two OLS models 
# df.OLS_Sens.top50<-left_join(df.OLS_Sens.top50, df.Yield.2model, by = 'Name') # merge constiuent yield:
# df.OLS_Sens.top50<-left_join(df.OLS_Sens.top50, df.G2.reduced, by = c('Name'='STAID')) # merge back the watershed characteristic data to this dataframe:
# n_sites<-dim(df.OLS_Sens.top50)[1] # set up variable for number of sites:
# df.cor.top50 <- df.OLS_Sens.top50 %>% # use the corrr package to correlate() and focus() on your variable of choice.
#   correlate(method = 'spearman') %>%
#   focus(c(OLS.I, OLS.S, Sen.I, Sen.S, Yield))%>%
#   pivot_longer(cols= c(2:6), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
#   mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
#   mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
#   drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero
# l.cor.top50<-df.cor.top50 %>% # create a list of each CQ parameter (4: OLS and Sens slope and intercept)and format it for ggplotting:
#   split(., df.cor.top50$CQ_Parameter)%>%
#   lapply(., \(i) i%>% 
#            arrange(Spearman_Correlation)%>%  
#            slice_head(n = 7)%>%
#            bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
#            mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
#            mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
#            filter(!between(Spearman_Correlation, -0.25,.25))) # Order by correlation strength
# # note*: need to add thisstep for top50 workflow:
# # when removing CDL_WWHT_SOY_DBL_CROP from pred to keep (it is clearly not good to include when looing at the univariate plots)
# # OLS.S and SEN.S have no correlates between -0.25 and 0.25
# # thus I need toremove those from l.cor:
# l.cor.top50<-purrr::keep(l.cor.top50, ~ nrow(.x) > 0)
# plist<-lapply(l.cor.top50, \(i) i%>%ggplot(aes(x = term, y = Spearman_Correlation, color = sig_0.05)) + # make plot list using lapply:
#                 geom_bar(stat = "identity") +
#                 scale_color_manual(values = c("not" = "red", "sig" = "blue"),na.value = NA, drop = FALSE)+
#                 facet_wrap('CQ_Parameter')+
#                 ylab(paste('Spearman Correlation')) +
#                 xlab("Watershed Attribute")+
#                 theme(axis.text.x=element_text(angle=40,hjust=1))+
#                 theme(legend.position="bottom"))
# p1<-ggpubr::ggarrange(plotlist = plist, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") # arrange plots:
# p1<-annotate_figure(p1, top = text_grob(paste(ncode, "Correlated Against Gauges 2"), color = "red", face = "bold", size = 14)) # add plot title:
# p1 # plot:
# 
# # look at univarite plot for AANY:
# 
# OLS<-l.cor.top50[[3]]%>%arrange(desc(Spearman_Correlation))
# 
# # make univariate plots (facets) of Y (determined in [[]]) above and predictors: 
# # note the predictor column isturned into an ordered factor based on thespearman correlation values
# # to makethe facet plots in order:
# 
# df.OLS_Sens.top50%>%
#   select(Name, OLS$CQ_Parameter[1], OLS$term)%>%
#   pivot_longer(cols = c(3:last_col()), names_to = 'Type', values_to = 'Value')%>%
#   drop_na(Value)%>%
#   # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
#   left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
#   mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
#   ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
#   geom_smooth(method = 'lm')+
#   geom_point()+
#   # scale_x_log10(
#   #   breaks = scales::trans_breaks("log10", function(x) 10^x),
#   #   labels = scales::trans_format("log10", scales::math_format(10^.x))
#   # )+
#   # scale_y_log10(
#   #   breaks = scales::trans_breaks("log10", function(x) 10^x),
#   #   labels = scales::trans_format("log10", scales::math_format(10^.x))
#   # )+
#   facet_wrap('Type', scales = 'free')+
#   ggtitle(paste('Univariate plots of', OLS$CQ_Parameter[1], 'against Top Gauges 2 Correlates for', OLS$CQ_Parameter[1], 'ordered from highest to lowest Spearman Rank Correlation' ))
# 
# #

# end comment out



























#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Correlations with datalayers ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# read in datalayers predictors:

load('Processed_Data/df.datalayers.62.Rdata')

# rename to just df.datalayers to work with code:

df.datalayers<-df.datalayers.62

# *note* that one site has NA for pedictors that require soils since it had no sSurgo coverage
# I will remove that site:

df.datalayers<-df.datalayers%>%filter(!is.na(CSA_perc))

# *note* df.datalayers from second pass does not have column "Dbl_Crop_WinWht/Soybeans"

# *note* df.datalayers has 61 sites now. 
# It did have 62 sites while df.TP_CQ (for second pass) has 63 sites
# this is not a big deal for the correlaitons because it will just omit the
# row with NAs in all the predictor columns, but this row will need to be omitted
# when running the MLR

# as per a first run of this workflow, remove emergentwetlands, soybeans, winter wheat, grassland
# and replace NA with zero (this is ok since I removed the site with NA for soils where a zero would be wrong):

df.datalayers<-df.datalayers%>%
  select(-c(R_EMERGWETNLCD06, R_SNOWICENLCD06, R_SHRUBNLCD06, R_BARRENNLCD06,R_GRASSNLCD06, Soybeans, Winter_Wheat, `Dbl_Crop_WinWht/Soybeans`, R_WOODYWETNLCD06))%>%
  replace(is.na(.), 0)

# look at correlaiton matrix:

library(corrplot)

corrplot(cor(df.datalayers[,-c(1, 10:30)]))

# now ready to run through correlaitons workflow:

# combined the reduced predictors df with the already filtered CQ data:

df.Correlation<-df.TP_CQ%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))%>%
  left_join(., df.datalayers, by = 'Name')

# create a dataframe with OLS and sens slope intercept and slopes:

df.OLS<-df.Correlation%>%
  group_by(Name)%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.I = OLS.co[1], 
            OLS.S = OLS.co[2])
  }) %>%
  ungroup

df.Sens<-df.Correlation%>%
  group_by(Name)%>%
  do({ Sens.co<-zyp.sen(log_C~log_Q,.)
  summarize(., Sen.I = Sens.co$coefficients[[1]],
            Sen.S= Sens.co$coefficients[[2]])
  }) %>%
  ungroup

# merge OLS and Sens:

df.OLS_Sens<-left_join(df.OLS, df.Sens, by = 'Name')

# merge constiuent yield:

df.OLS_Sens<-left_join(df.OLS_Sens, df.Yield, by = 'Name')

# merge back the watershed characteristic data to this dataframe:

df.OLS_Sens<-left_join(df.OLS_Sens, df.datalayers, by = 'Name')

# now run correlations between intercepts and slopes and watershed characteristics. to do this: (I orginally did this workflow using n_months (C:\PhD\Research\Mohawk\Code\Mohawk_Regression-analyizing_predictor_df.R)

# set up variable for number of sites:

n_sites<-dim(df.OLS_Sens)[1] 

# use the corrr package to correlate() and focus() on your variable of choice.

df.cor <- df.OLS_Sens %>% 
  correlate(method = 'spearman') %>%
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S, Yield))%>%
  pivot_longer(cols= c(2:6), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# then plot results: to do this:

# create a list of each CQ parameter (4: OLS and Sens slope and intercept)and format it for ggplotting:

l.cor<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(Spearman_Correlation)%>%  
           slice_head(n = 7)%>%
           bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           filter(!between(Spearman_Correlation, -0.25,.25))) # Order by correlation strength

# make plot list using lapply:

plist<-lapply(l.cor, \(i) i%>%ggplot(aes(x = term, y = Spearman_Correlation, color = sig_0.05)) +
                geom_bar(stat = "identity") +
                scale_color_manual(values = c("not" = "red", "sig" = "blue"),na.value = NA, drop = FALSE)+
                facet_wrap('CQ_Parameter')+
                ylab(paste('Spearman Correlation')) +
                xlab("Watershed Attribute")+
                theme(axis.text.x=element_text(angle=40,hjust=1))+
                theme(legend.position="bottom"))

# arrange plots:

p1<-ggpubr::ggarrange(plotlist = plist, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

# add title:

p1<-annotate_figure(p1, top = text_grob(paste(ncode, "Correlated Against Gauges 2"),
                                        color = "red", face = "bold", size = 14))
# plot:

p1

# save(p1, file = 'Processed_Data/p1.Rdata')

# univariate plots of the top correlates/small land use types:

# Lets do plots of OLS intercept: to do this:

# determine the top correlates (the numberof the list determined which CQ parameter)
# here 1 = OLS intercept:

OLS<-l.cor[[1]]%>%arrange(desc(Spearman_Correlation))

# make univariate plots (facets) of Y (determined in [[]]) above and predictors: 
# note the predictor column isturned into an ordered factor based on thespearman correlation values
# to makethe facet plots in order:

df.OLS_Sens%>%
  select(Name, OLS$CQ_Parameter[1], OLS$term)%>%
  pivot_longer(cols = c(3:last_col()), names_to = 'Type', values_to = 'Value')%>%
  drop_na(Value)%>%
  # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
  left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
  mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
  ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Type', scales = 'free')+
  ggtitle(paste('Univariate plots of', OLS$CQ_Parameter[1], 'against Top Gauges 2 Correlates for', OLS$CQ_Parameter[1], 'ordered from highest to lowest Spearman Rank Correlation' ))

# create a function to use lapply to plot all 4 univariate sets:

# number<-5

# make_plot<-function(number){
#   
#   OLS<-l.cor[[number]]%>%arrange(desc(Spearman_Correlation))
#   
#   x<-df.OLS_Sens%>%
#     select(Name, OLS$CQ_Parameter[1], OLS$term)%>%
#     pivot_longer(cols = c(3:last_col()), names_to = 'Type', values_to = 'Value')%>%
#     drop_na(Value)%>%
#     # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
#     left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
#     mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
#     ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
#     geom_smooth(method = 'lm')+
#     geom_point()+
#     facet_wrap('Type', scales = 'free')+
#     ggtitle(paste('Univariate plots of', OLS$CQ_Parameter[1], 'against Top Gauges 2 Correlates for', OLS$CQ_Parameter[1], 'ordered from highest to lowest Spearman Rank Correlation' ))
#   
#   x
#   
# }

# use function in lapply (clear plot list first):

# lapply(c(1:5), \(i) make_plot(i))

#


























#### Categorizing land use - need to recode to accomdate G2 and datalayers land use predictors...  ####

# # note* this section should come after regular old correlations such that df.OLS_Sens has the G2 predictors:
# 
# # USGS criteria:
# # Agricultural sites have >50% agricultural land and ≤5% urban land;
# # urban sites have >25% urban and ≤25% agricultural land;
# # undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land;
# # all other combinations of urban, agricultural, and undeveloped lands are classified as mixed
# 
# # initate a dataframe with the landuse values from gauges 2:
# 
# df.NLCD06<-distinct(df.OLS_Sens, Name, .keep_all = T)
# 
# # I confirmed that pasture + crops = plant:
# 
# (df.NLCD06$PASTURENLCD06+df.NLCD06$CROPSNLCD06)==df.NLCD06$PLANTNLCD06
# 
# # and that the ones that start with DEV sum up to the first DEV one.
# 
# # create the land use class column based on USGS critiera:
# # adjusted thresholds:
# 
# df.NLCD06<-df.NLCD06%>%
#   mutate(USGS.LU.Adjusted = case_when(.default = 'Mixed',
#                                       PLANTNLCD06 > 30 & DEVNLCD06 <= 10 ~ 'Agriculture',
#                                       DEVNLCD06 > 10 & PLANTNLCD06 <= 30 ~ 'Urban',
#                                       DEVNLCD06+PLANTNLCD06 <= 10 ~ 'Undeveloped'))
# 
# # frequency table:
# 
# table(df.NLCD06$USGS.LU.Adjusted)
# 
# # boxplots of
# 
# p<-df.NLCD06%>%
#   pivot_longer(cols = 2:5, names_to = 'CQ_parameter', values_to = 'Value')%>%
#   ggplot(., aes(x=USGS.LU.Adjusted, y=Value, color =USGS.LU.Adjusted ))+
#   geom_boxplot(varwidth = TRUE, alpha=0.2)+
#   # scale_x_discrete(labels=my_xlab)+
#   facet_wrap('CQ_parameter', scales = 'free')+
#   stat_compare_means(method = "anova")+      # Add global p-value # , label.y = max(df.datalayers$Value) ???
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "0.5") +
#   ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')
# 
# # p
# 
# #




































#### Grouping CQ curves (stationary, mobilization, dilutionary, complex) ####

# cant run next few section until categorizing landuse works (need dfNLCD...)

# # create a list of dataframes for each sites CQ observations:
# 
# df.TP_CQ<-df.TP_CQ%>%
#   rename(Name = site_no)%>%
#   mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
#   filter(is.finite(log_C))%>%
#   filter(is.finite(log_Q))
# 
# l.TP_CQ<-df.TP_CQ%>%
#   split(., .$Name)
# 
# # create lm models for each site:
# 
# l.TP_CQ.lm<-lapply(l.TP_CQ, \(i) lm(log_C~log_Q, data=i))
# 
# # save the model coef ad pvals:
# 
# TP_CQ.coef<-as.data.frame(bind_rows(lapply(l.TP_CQ.lm, \(i) summary(i)$coefficients[,1]), .id = 'site_no'))%>%
#   rename(Intercept = 2, Slope = 3)
# 
# # save the pvalues 
# 
# TP_CQ.pvals<-as.data.frame(bind_rows(lapply(l.TP_CQ.lm, \(i) summary(i)$coefficients[,4]), .id = 'site_no'))%>%
#   rename(Intercept.pval = 2, Slope.pval = 3)
# 
# # merge the two dfs:
# 
# df.lm<-left_join(TP_CQ.coef,TP_CQ.pvals,by='site_no')%>%left_join(., df.TP_CQ%>%select(Name, n_sample_rank)%>%distinct(., .keep_all = T), by = c('site_no' = 'Name'))
# 
# # add column for CQ type:
# 
# df.lm<-mutate(df.lm, Type = ifelse(Slope.pval>0.05, 'Stationary', ifelse(Slope>0, 'Mobilization', 'Dilution')))
# 
# # merge CQ type labels to the df for plotting
# 
# df.TP_CQ<-left_join(df.TP_CQ,df.lm%>%select(site_no, Type),by=c('Name'='site_no'))
# 
# # create a df for CQplot of all sites with BP analysis: make sure to filter out pre 2001 samples:
# 
# df_Seg.2<-filter(df_Seg, site %in% df.lm$site_no)%>%
#   left_join(.,df.lm%>%select(site_no, Type),by=c('site'='site_no'))
# 
# # # make plot:
# # 
# # p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
# #   geom_point(aes(color = Type))+
# #   geom_smooth(method = 'lm')+
# #   geom_line(aes(x = Q, y = Seg_C), color = 'purple', size = 2)+
# #   facet_wrap('n_sample_rank', scales = 'free')+
# #   theme(
# #     strip.background = element_blank(),
# #     strip.text.x = element_blank()
# #   )
# 
# # p
# 
# # looking at this plot I want to add a fourth CQ type for complex, if the slopes of the BP analysis look widely different. 
# # I will start with calcuating the angle between pre-post BP slope and see which sites get signled out:
# 
# # add a new column with the angle between the two lines:
# 
# df_Seg.2<-df_Seg.2%>%
#   mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))
# 
# # I want to add land use as color as well. Merge the plotting df with the land use df:
# 
# df_Seg.2<-left_join(df_Seg.2, df.NLCD06%>%select(Name, USGS.LU.Adjusted), by = c('site'='Name'))
# 
# # create color pallete for the slope angle:
# 
# hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)
# 
# # make plot:
# 
# p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
#   geom_point(aes(color = Type))+
#   scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
#   ylab('log(TP)')+
#   geom_smooth(method = 'lm')+
#   new_scale_color() +
#   geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
#   geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
#   scale_color_manual(name = "Slope Angle", values = hc)+
#   facet_wrap(dplyr::vars(n_sample_rank), scales = 'fixed')+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )+
#   geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .35)+
#   scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","yellow", "green"))
# 
# # p
# 
# # export df_Seg.2:
# 
# # save(df_Seg.2, file = 'Processed_Data/TP.df_Seg.2.Rdata')
# 
# # based on this plot, I think I need to zoom in on each one:
# 
# # create a plotting function for each site:
# 
# # p.fun<-function(df){
# #   ggplot(df, aes(x = log(Q_real), y = log(C)))+
# #     geom_point(aes(color = Type))+
# #     scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
# #     geom_smooth(method = 'lm')+
# #     new_scale_color() +
# #     geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
# #     geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
# #     scale_color_manual(name = "Slope Angle", values = hc)+
# #     ggtitle(df_Seg.2$site[df_Seg.2$n_sample_rank==i])
# # }
# 
# # use function in lappy to make lots of plots (clear plot list first)
# 
# # lapply(sort(unique(df_Seg.2$n_sample_rank)), \(i) df_Seg.2%>%filter(n_sample_rank==i)%>%p.fun(.))
# 
# # now looking at the expanded plots for each site, I feel that
# # non really could be justified as complex! it feels like the breakpoint
# # analysis lines aren't 'real'
# 
# # My first thought is to color the CQ points based on season, AMC, etc
# # this is so overwhleming. 
# 
# # I think I'm going to pause here in the analysis, and carry on with the
# # other constituents. Once I have those up to this point I can check in with Chuck and Steve






















#### Grouping CQ curves using top 50 ####

# # determine the median flow rate for each sites CQ observations: to do this:
# 
# df.TP_CQ.top50<-split(df.TP_CQ, f = df.TP_CQ$Name)%>% # split the df into list of dfs for each site, 
#   lapply(., \(i) i%>%filter(X_00060_00003>median(i$X_00060_00003)))%>% # filter to the median flow rate for each site,
#   bind_rows(.) # bind back into single dataframe:
# 
# # now run through grouping CQ workflow with this truncated df:
# 
# l.TP_CQ<-df.TP_CQ.top50%>%split(., .$Name) # resplit the df into list:
# l.TP_CQ.lm<-lapply(l.TP_CQ, \(i) lm(log_C~log_Q, data=i)) # create lm models for each site:
# TP_CQ.coef<-as.data.frame(bind_rows(lapply(l.TP_CQ.lm, \(i) summary(i)$coefficients[,1]), .id = 'site_no'))%>%rename(Intercept = 2, Slope = 3) # save the model coef ad pvals:
# TP_CQ.pvals<-as.data.frame(bind_rows(lapply(l.TP_CQ.lm, \(i) summary(i)$coefficients[,4]), .id = 'site_no'))%>%rename(Intercept.pval = 2, Slope.pval = 3) # save the pvalues 
# df.lm.top50<-left_join(TP_CQ.coef,TP_CQ.pvals,by='site_no')%>%left_join(., df.TP_CQ%>%select(Name, n_sample_rank)%>%distinct(., .keep_all = T), by = c('site_no' = 'Name')) # merge the two dfs:
# df.lm.top50<-mutate(df.lm.top50, Type_top50 = ifelse(Slope.pval>0.05, 'Stationary', ifelse(Slope>0, 'Mobilization', 'Dilution'))) # add column for CQ type:
# df.TP_CQ.top50<-left_join(df.TP_CQ.top50,df.lm.top50%>%select(site_no, Type_top50),by=c('Name'='site_no')) # merge CQ type labels to the df for plotting
# 
# # here is where it changes: add another column for the y coordinates of the OLS line for just the top 50 percentile:
# 
# temp<-lapply(seq_along(l.TP_CQ.lm), \(i) data.frame(Seg_C_log_top50 = fitted(l.TP_CQ.lm[[i]]), l.TP_CQ[[i]]$sample_dt))%>%
#   purrr::set_names(names(l.TP_CQ.lm))%>%
#   bind_rows(., .id = 'Name')%>%
#   rename(site = Name, Date = 3)%>%
#   mutate(Date = as.Date(Date))
# 
# # now add upper 50 CQ type and Seg_C_log_top50 columns to df_Seg.2:
# 
# df_Seg.2<-left_join(df_Seg.2, df.lm.top50%>%select(site_no, Type_top50), by = c('site'='site_no'))%>%
#   mutate(Date = as.Date(Date))%>%
#   left_join(., temp, by = c('site', 'Date'))
# 
# # make plot:
# 
# p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
#   geom_point()+
#   geom_smooth(aes(color = Type),method = 'lm')+
#   scale_color_manual(name = "Full CQ Type", values = c("red", "blue", "green"))+
#   ylab('log(TP)')+
#   new_scale_color() +
#   geom_line(aes(x = log(Q_real), y = Seg_C_log_top50), size = 2.5, color = 'black')+
#   geom_line(aes(x = log(Q_real), y = Seg_C_log_top50, color = Type_top50), size = 2)+
#   scale_color_manual(name = "Top 50 CQ Type", values = c("blue", "green"))+
#   facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )
# 
# # p
# 
# # export df_Seg.2 for use in FingerLakesPresentation.R
# 
# # save(df_Seg.2, file = 'Processed_Data/TP.df_Seg.2.w_top50.Rdata')
# 
# #
























#### Triangle plot ####

# # set up dataframe:
# 
# # add percent ag and urban to df.lm:
# 
# df.tri<-left_join(df.lm, df.NLCD06%>%select(Name, DEVNLCD06, PLANTNLCD06), by = c('site_no'='Name'))%>%select(site_no, Slope, Type, DEVNLCD06, PLANTNLCD06)%>%mutate(Slope = round(Slope, 2), DEVNLCD06=round(DEVNLCD06, 2), PLANTNLCD06=round(PLANTNLCD06, 2))
# 
# # create plot:
# 
# p<-ggplot(df.tri, aes(x=DEVNLCD06, y=PLANTNLCD06)) +
#   geom_point(aes(shape=as.factor(Type), fill=abs(Slope)), size = 4) +
#   scale_shape_manual(values = c("Mobilization" = 24, "Dilution" = 25, "Stationary" = 22),
#                      guide = guide_legend(override.aes = list(fill = "pink")))+
#   scale_fill_gradient(low = "yellow", high = "red")+
#   xlab('Percent Developed')+
#   ylab('Percent Agriculture')+
#   ggtitle(paste(ncode, 'CQ type and slope magnitude as a function of percent Ag and Developed Land'))
# 
# p$labels$fill <- "Slope Magnitude"
# p$labels$shape <- "CQ Type"
# 
# p
# 
# # want to recreate this plot in Code/FingerLakesPresentation.R so exporting df.tri:
# 
# # save(df.tri, file = 'Processed_Data/df.tri.Rdata')
# 
# #


























#### Triangle plot with top 50 ####

# df.tri<-left_join(df.lm.top50, df.NLCD06%>%select(Name, DEVNLCD06, PLANTNLCD06), by = c('site_no'='Name'))%>%select(site_no, Slope, Type_top50, DEVNLCD06, PLANTNLCD06)%>%mutate(Slope = round(Slope, 2), DEVNLCD06=round(DEVNLCD06, 2), PLANTNLCD06=round(PLANTNLCD06, 2))
# 
# # create plot:
# 
# p<-ggplot(df.tri, aes(x=DEVNLCD06, y=PLANTNLCD06)) +
#   geom_point(aes(shape=as.factor(Type_top50), fill=abs(Slope)), size = 4) +
#   scale_shape_manual(values = c("Mobilization" = 24, "Dilution" = 25, "Stationary" = 22),
#                      guide = guide_legend(override.aes = list(fill = "pink")))+
#   scale_fill_gradient(low = "yellow", high = "red")+
#   xlab('Percent Developed')+
#   ylab('Percent Agriculture')+
#   ggtitle(paste(ncode, 'CQ type and slope magnitude as a function of percent Ag and Developed Land'))
# 
# p$labels$fill <- "Slope Magnitude"
# p$labels$shape <- "CQ Type"
# 
# p
# 
# # want to recreate this plot in Code/FingerLakesPresentation.R so exporting df.tri:
# 
# # save(df.tri, file = 'Processed_Data/df.tri.Rdata')
# 
# #

























#### MLR ####

# note that this will run with the based on the last run correlations section
# thus, it will run withdf.datalayers:

# create a list of dfs, one for each cQ parameter, with the values of the spearman correlations for each predictor:

l.cor.MLR<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(desc(abs(Spearman_Correlation)))%>%  
           # slice_head(n = 7)%>%
           # bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           # distinct(term, .keep_all = T)%>%
           # filter(!between(Spearman_Correlation, -0.25,.25))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           # mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           # filter(sig_0.05 == 'sig')%>%
           as.data.frame())


x<-l.cor.MLR[[1]]

# create a list of dfs from df.OLS.Sens, which contains the values of the predictors, subseting using the names in each l.cor[[i]]$term:

l.cor.MLR.full<-lapply(1:5, \(i) df.OLS_Sens%>%
                         select(i+1, l.cor.MLR[[i]]$term)%>%
                         as.data.frame()%>%
                         rename(term = 1))

names(l.cor.MLR.full)<-names(df.OLS_Sens)[2:6]

y<-l.cor.MLR.full[[1]]

y.round<-round(y, 2)

# create another list but keep the site names for use at end of loop in site outliers:

l.cor.MLR.full.w_name<-lapply(1:5, \(i) df.OLS_Sens%>%
                                select(1, i+1, l.cor.MLR[[i]]$term)%>%
                                as.data.frame())

names(l.cor.MLR.full.w_name)<-names(df.OLS_Sens)[2:6]

z<-l.cor.MLR.full.w_name[[1]]

# Trying no feature selection:

# loop through the 5 CQ parameters:

i<-1

m.list<-list()

for (i in 1:5){
  
  df<-l.cor.MLR.full[[i]]%>%
    drop_na(CSA_perc) # in the second pass, this is needed to remove the 1 site that did not have soils data nd did not carry over into df.datalayers (but is in df.TP_CQ)
  
  library(caret)
  
  set.seed(123)
  
  # Set up repeated k-fold cross-validation
  
  train.control <- trainControl(method = "cv", number = 10)
  
  # Train the model
  
  step.model <- train(term ~., data = df,
                      method = "leapForward", 
                      tuneGrid = data.frame(nvmax = 1:5),
                      trControl = train.control
  )
  
  
  # look at model results:
  
  step.model$results
  step.model$bestTune
  summary(step.model$finalModel)
  coef(step.model$finalModel, as.numeric(step.model$bestTune))
  
  # get df of just the predictors in this model:
  
  pred<-names(coef(step.model$finalModel, as.numeric(step.model$bestTune)))[-1]
  
  # remove potential backticks:
  
  pred<-gsub("`", "", pred)
  
  # select down df:
  
  df.2<-df%>%select(term, pred)
  
  # make lm:
  
  m<-lm(term~., data=df.2)
  
  m.list[[i]]<-m
  
  # summary(m)
  
  # look at univarite plots:
  
  p<-df.2%>%pivot_longer(cols = 2:last_col(), names_to = 'Predictor', values_to = 'Value')%>%
    rename(!!names(l.cor.MLR.full)[i]:=term)%>%
    ggplot(., aes(x=Value, y=!!sym(names(l.cor.MLR.full)[i])))+
    geom_point()+
    geom_smooth(method = 'lm')+
    facet_wrap('Predictor', scales = 'free')
  
  # look at model residuals:
  
  # plot(m)
    
}

# set m.list names:

names(m.list)<-names(l.cor.MLR.full)

tab_model(m.list, dv.labels = names(m.list), title = paste('Comparison of MLR models for',ncode, 'using forward selection implemented in caret::train'), file="temp.html")

# export m.list for use in Code/FingerLakesPresentation.R:

# save(m.list, file = 'Processed_Data/m.list.TP.Rdata')

#



























#### MLR with observations greater than the median ####

# # determine the median flow rate for each sites CQ observations: to do this:
# 
# df.TP_CQ.top50<-split(df.TP_CQ, f = df.TP_CQ$Name)%>% # split the df into list of dfs for each site, 
#   lapply(., \(i) i%>%filter(X_00060_00003>median(i$X_00060_00003)))%>% # filter to the median flow rate for each site,
#   bind_rows(.) # bind back into single dataframe:
# 
# # now run through the MLR workflow:
# 
# l.cor.MLR<-df.cor.top50 %>% # create a list of dfs, one for each cQ parameter, with the values of the spearman correlations for each predictor:
#   split(., df.cor.top50$CQ_Parameter)%>%
#   lapply(., \(i) i%>% 
#            arrange(desc(abs(Spearman_Correlation)))%>%  
#            # slice_head(n = 7)%>%
#            # bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
#            # distinct(term, .keep_all = T)%>%
#            # filter(!between(Spearman_Correlation, -0.25,.25))%>%
#            mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
#            # mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
#            # filter(sig_0.05 == 'sig')%>%
#            as.data.frame())
# l.cor.MLR.full<-lapply(1:5, \(i) df.OLS_Sens.top50%>% # create a list of dfs from df.OLS.Sens, which contains the values of the predictors, subseting using the names in each l.cor[[i]]$term:
#                          select(i+1, l.cor.MLR[[i]]$term)%>%
#                          as.data.frame()%>%
#                          rename(term = 1))
# names(l.cor.MLR.full)<-names(df.OLS_Sens.top50)[2:6]
# l.cor.MLR.full.w_name<-lapply(1:5, \(i) df.OLS_Sens.top50%>% # create another list but keep the site names for use at end of loop in site outliers:
#                                 select(1, i+1, l.cor.MLR[[i]]$term)%>%
#                                 as.data.frame())
# names(l.cor.MLR.full.w_name)<-names(df.OLS_Sens.top50)[2:6]
# m.list<-list() # create list to store model objects
# i<-1
# for (i in 1:5){ # loop through the 5 CQ parameters:
#   df<-l.cor.MLR.full[[i]]
#   library(caret)
#   set.seed(123)
#   train.control <- trainControl(method = "cv", number = 10) # Set up repeated k-fold cross-validation
#   step.model <- train(term ~., data = df, method = "leapForward", tuneGrid = data.frame(nvmax = 1:5),trControl = train.control) # Train the model
#   step.model$results  # look at model results:
#   step.model$bestTune
#   summary(step.model$finalModel)
#   coef(step.model$finalModel, as.numeric(step.model$bestTune))
#   pred<-names(coef(step.model$finalModel, as.numeric(step.model$bestTune)))[-1] # get df of just the predictors in this model:
#   df.2<-df%>%select(term, pred)
#   m<-lm(term~., data=df.2) # make lm:
#   m.list[[i]]<-m # append model to list
#   summary(m)
#   p<-df.2%>%pivot_longer(cols = 2:last_col(), names_to = 'Predictor', values_to = 'Value')%>% # look at univarite plots:
#     rename(!!names(l.cor.MLR.full)[i]:=term)%>%
#     ggplot(., aes(x=Value, y=!!sym(names(l.cor.MLR.full)[i])))+
#     geom_point()+
#     geom_smooth(method = 'lm')+
#     facet_wrap('Predictor', scales = 'free')
#   p
#   # plot(m)  # look at model residuals:
#   
# }
# 
# # set m.list names:
# 
# names(m.list)<-names(l.cor.MLR.full)
# 
# tab_model(m.list, dv.labels = names(m.list), title = paste('Comparison of MLR models for',ncode, 'using forward selection implemented in caret::train'), file="temp.top50.html")
# 
# #

























#### ~PLSR~ ####

# install.packages("pls")
library(pls)
library(nortest)
library(flipFormat)
library(caret)

# run PLSR for 1 CQ metric: to do this:

# set up df:

df.OLS.I<-df.OLS_Sens[,c(6, 7:ncol(df.OLS_Sens))]%>%
  drop_na(CSA_perc)%>%
  as.data.frame()%>%
  select_if(colSums(.) != 0)

# standardize:

df.OLS.I <- as.data.frame(lapply(df.OLS.I, scale, center = TRUE, scale = TRUE))
  

#### PLS package ####

# test for normality:

df.norm.test<-lapply(df.OLS.I, lillie.test)%>% # use nortest::lillie.test in lapply,
  lapply(., \(i) round(i$p.value,3))%>% # extract just pvalue
  enframe(.) # convert list to df

# I cant boxcox transform variables with significant p-values (reject null hypothesis that data is normally distirubed):
# because of zeros and negatives...
# I was doing this because Musolff et al. 2015 said they normalized variable prior to using plsr
# however the pls user manual doesnt say this and the example data they use is not normally distributed.

# so I am just going to run the plsr:

# run:

test <- plsr(Yield ~ ., data = df.OLS.I, validation = "CV")

# summary(test)

# plot(RMSEP(test), legendpos = "topright")

# plot(test, plottype = "correlation", ncomp=1:2, legendpos = "bottomleft",labels = "numbers")

plot(test, labels = "numbers")

# if you rerun test, consistently 31, 16, 11, 40, and 60 are outliers from the pack
# somtimes 11 is in cloud...

# but removing these and trying again:

df.OLS.I<-df.OLS.I[-c(31, 16, 11, 40, 60),]

load('Processed_Data/NWIS_Watershed_Shapefiles.62.Rdata')
temp<-df.sf.NWIS.62%>%filter(Name %in% df.OLS_Sens$Name[c(31, 16, 11, 40, 60)])
mapview(temp)

# another method for ncomp determinitaiton:

ncomp.onesigma <- selectNcomp(test, method = "onesigma", plot = TRUE)
ncomp.permut <- selectNcomp(test, method = "randomization", plot = TRUE)

# extract the ideal number of comps (minimum RMSEP):

ncomp<-which.min(RMSEP(test)$val[estimate = "adjCV", , ]) - 1

# rerun:

# test <- plsr(OLS.I ~ ., ncomp = ncomp, data = df.OLS.I, validation = "CV")

# extract VIP:

coef<-as.data.frame.table(coef(test))%>%
  select(1,4)%>%
  mutate(VIP = round(Freq/sum(Freq), 3), Value = round(Freq, 3))%>%
  rename(Variable = 1)%>%
  select(Variable, Value, VIP)%>%
  arrange(desc(Value))

#### Caret package ####

# Build the model on training set
set.seed(123)

model <- train(
  Yield ~ ., data = df.OLS.I, method = "rf",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)

# Plot model RMSE vs different values of components

plot(model)

# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE

model$bestTune

# Summarize the final model

summary(model$finalModel)

# get variable importance

v <- varImp(model,scale = TRUE)[["importance"]]
v$Overall <- v$Overall / sum(v$Overall)
v<-v%>%arrange(Overall)

#























#### Checking Dam Storage and Corn ####

# It was found that dam storage was showing up in the MLR model (done later in the code)
# And that corn was shwoing up in the later editions of the MLR with the reduced datasets
# so this seciton it to see if thetwo correlate and also make maps of dam storage, corn, etc.






























#### 













# I would chose the following sites as complex:

complex_sites<-unique(df_Seg.2$site)[c(3,16,21,22,29,37)] 

# old complex site (when n=53): c(3,18,24,25,26,35,42,43)

df_Seg.2<-mutate(df_Seg.2, Type = ifelse(site %in% complex_sites, 'Complex', Type))

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p


# plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

# plot of trouble sites from MLR analysis:


p<-filter(df_Seg.2, site %in% look_at_sites)%>%
  ggplot(., aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(site), scales = 'free')+
  # theme(
  #   strip.background = element_blank(),
  #   strip.text.x = element_blank()
  # )+
  geom_rect(data = .%>%distinct(.$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

# save(p1, file = 'Processed_Data/p1.Rdata') 

# add CQ type to the mapping df with polygon trasperncy based on number of samples:

map.CQ_Type<-df.sf.NWIS.keep%>%
  filter(Name %in% df.datalayers$Name)%>%
  left_join(.,distinct(df_Seg.2, site, .keep_all = T)%>%select(.,c(site, Type, n_sample_rank)),  by = c('Name'='site'))%>%
  select(Name, Type, n_sample_rank)%>%
  arrange(n_sample_rank)%>%
  mutate(NEW = (1/row_number())*40)
# left_join(.,df.NWIS.TP%>%group_by(site_no)%>%summarise(n=n()), by = c('site'='site_no'))%>%
# mutate(NEW = case_when(n < 50 ~ .05,
#                        n >=50 & n < 100 ~ .25,
#                        n >=100 & n < 500 ~ .5,
#                        n >=500 ~ .9))
# map:

mapview(map.CQ_Type, zcol = 'Type', alpha.regions = 'NEW')

# map of trouble sites in MLR section:

x<-map.CQ_Type%>%filter(Name %in% look_at_sites)

mapview(x, zcol = 'Type')

######################################################
#### ~~~~ XY plot of Slopes and Intercepts ~~~~ ####
####################################################

# plot of slopes and intercepts for each CQ type:

m<-m%>%mutate(Type = ifelse(site_no %in% complex_sites, 'Complex', Type))%>%
  left_join(., df.datalayers%>%select(Name, USGS.LU.Adjusted), by = c('site_no'='Name'))

ggplot(m, aes(x=`(Intercept)`, y=log_Q, color = Type))+
  geom_point(size = 2)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  facet_wrap('USGS.LU.Adjusted', scales = 'fixed')

# I want to do this with more sites























#### Land Use changes across watersheds (not run) ####

# # Look at changes in NLCD and CDL calculated Ag %land between 2008-2020 for each 53 watershed
# # to do this:
# 
# # pivot longerto getthe land use catgories in a single column:
# 
# df.LU.2<-pivot_longer(df.LU, cols = 6:12, names_to = 'LandUse', values_to = 'values')
# 
# # pivot wider to get the land use source (data layers) in their own columns:
# 
# df.LU.2<-pivot_wider(df.LU.2, names_from = 'LU_source', values_from = 'values')
# 
# # pivot longer to put CDL and NLCD in own columns:
# 
# df.LU.CDL<-pivot_longer(df.LU.2, cols = starts_with('CDL'), names_to = 'CDL', values_to = 'CDL_values')
# df.LU.NLCD<-pivot_longer(df.LU.2, cols = starts_with('NLCD'), names_to = 'NLCD', values_to = 'NLCD_values')
# 
# # group by land use and calculate the difference between the years:
# 
# df.LU.2<-df.LU.2%>%mutate(CDL_diff = abs(round(CDL_2020-CDL_2008,2)),NLCD_diff = abs(round(NLCD_2019-NLCD_2001,2)))%>%
#   arrange(desc(CDL_diff))
# 
























#### Seasonal ANanlysis ####

# is the CQ relaitonship different between seasons at these sites?

# I need to filter down to sites with over 100 samples for this, 20 samples isnt going tocut it:

df_Seg.3<-filter(df_Seg.2, n>100)

# add season column:

df_Seg.3$Season<-getSeason(df_Seg.3$Date)

# plot:

p1<-ggplot(df_Seg.3, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Season), size = 4)+
  # scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.title=element_text(size=14)
  )+
  geom_rect(data = df_Seg.3%>%distinct(df_Seg.3$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = Type), alpha = .15)+
  scale_fill_manual(name = "CQ Type", values = c("red", "blue","purple", "green"))

p1

# save(p1, file = 'Processed_Data/p1.Rdata') 

# 


























############################################
#### ~~~~ Q Exceedence Proability ~~~~ ####
############################################




























####~~~~ Calcualte EP ~~~~####

# subset the daily flows to the 40 sites:

temp<-filter(df.NWIS.Q, site_no %in% df.datalayers$Name)

# split into lists by site:

temp<-split(temp, f = temp$site_no) 

# sort each dataframe in the list by the flow, add a column for m and n, convert back to a dataframe using bind_rows, and select only the needed columns for the left join (next step):

temp<-bind_rows(lapply(temp, \(i) i[order(i$X_00060_00003,decreasing = T),]%>%mutate(m = 1:n(), n = n())))%>%select(site_no, Date, m, n)

# append the value of M and n for each C-Q observation in the C-Q dataframe:
# need to rename the n column in the left dataframe as well:

df_Seg.2<-left_join(df_Seg.2%>%rename(n_samples = n), temp, by = c("site" = "site_no", "Date" = "Date"))

# calcualte the exceednce proability for Q for each C observation:

df_Seg.2$EP<-round(df_Seg.2$m/(df_Seg.2$n+1), 4)

# plot:

p1<-ggplot(df_Seg.2, aes(x = EP, y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "Log-Log\nCQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  # geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = .5, y =-2, label = CAFO_count), inherit.aes = FALSE, size=15)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  scale_x_reverse()+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p1

# save(p1, file = 'Processed_Data/p1.Rdata')


























####~~~~ Estimate Break point analysis using EP-Q ~~~~####

#  rerun the BP analysis using EP-Q:

l_Seg.EP<-list()
davies.test.matrix.EP<-NULL
EP<-df.NWIS.TP_CQ%>%
  filter(site_no %in% df.datalayers$Name)%>%
  left_join(., df_Seg.2%>%select(site, Date, EP), by = c('site_no'='site', 'sample_dt'='Date'))%>%
  drop_na(result_va, Q_yield)%>%
  rename(site=site_no, Date = sample_dt)
temp.loop<-sort(unique(EP$site))

i<-1

for (i in seq_along(temp.loop)){
  
  tryCatch({
    
    # print the site name for loop debugging:
    
    print(i)
    print(temp.loop[i])
    
    # create a dataframe that will work with segmented. To do this: 
    # filter for the site the loop is in
    # add log transformed C and Q columns, as well as duplicated columns for renamed C and Q
    # filter for real log C and Q values so breakpoint analysis works smoothly:
    
    df<-EP%>%
      filter(site == temp.loop[i])%>%
      mutate(C = result_va)%>%
      mutate(log_C = log(C))%>%
      filter(is.finite(log_C))%>%
      filter(is.finite(EP))
    
    # build a single slope lm for log C and Q. Tis model is also used inthe breakpoint analysis inthenext step:
    
    m<-lm(log_C~EP, df)
    
    # perform breakpoint regression:
    
    m_seg<-segmented(obj = m, npsi = 1)
    
    # perform davies test for constant linear predictor:
    # the results are saved as a string with the the site name and true/false:
    
    x<-paste(temp.loop[i], '-', davies.test(m)$p.val<0.05)
    
    # add the results of davies test to the matrix made prior to this for loop:
    
    davies.test.matrix.EP<-c(davies.test.matrix.EP,x)
    
    # save the breakpoints
    
    bp<-m_seg$psi[1]
    
    # save the slopes: To do this:
    # a conditional statement is needed since sometimes the segmented function wont fit a two slope model at all and will return a object that doesnt work with the slope function used here:
    
    if(length(class(m_seg))==2){
      s<-as.data.frame(slope(m_seg))
    } else{
      s<-NA
    }
    
    # get the intercepts (again conditional statement is needed):
    
    if(length(class(m_seg))==2){
      inter<-as.data.frame(intercept(m_seg))
    } else{
      inter<-NA
    }
    
    
    # get the model fitted data, put in a dataframe, and reformat this dataframe to export out of the loop
    
    if(length(class(m_seg))==2){
      fit <- data.frame(site = df$site, Date = df$Date, Seg_C_EP = fitted(m_seg), I1_EP = inter$Est.[1], I2_EP = inter$Est.[2], Slope1_EP = s[1,1], Slope2_EP = s[2,1], BP_EP = bp)
      result_df<-left_join(df, fit, by = c('site', 'Date'))
    } else{
      fit <- data.frame(site = df$site, Date = df$Date, Seg_C_EP = fitted(m_seg), I1_EP = inter$Est.[1], I2_EP = NA, Slope1_EP = NA, Slope2_EP = NA, BP_EP = NA)
      result_df<-left_join(df, fit, by = c('site', 'Date'))
    }
    
    l_Seg.EP[[i]]<-result_df
    
  }, 
  
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

davies.test.matrix.EP
df.davies.EP<-as.data.frame(davies.test.matrix.EP)%>%
  separate_wider_delim(1, "-", names = c("site", "BP_EP_yes"))%>%
  mutate(across(c(1,2), trimws))
df_Seg.EP<-bind_rows(l_Seg.EP)%>%
  mutate(Seg_EP=EP)
df_Seg.EP<-df_Seg.EP%>%
  left_join(., df.davies.EP, by = 'site')%>%
  mutate(across(c(Seg_C_EP, Seg_EP), ~replace(., BP_EP_yes == 'FALSE', NA)))

# merge the land use, CQ type, data to this df:

df_Seg.EP<-left_join(df_Seg.EP, m%>%select(site_no, Type, USGS.LU.Adjusted), by = c('site'='site_no'))

# plot:

p<-ggplot(df_Seg.EP, aes(x = EP, y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "Log-log\nCQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  new_scale_color()+
  geom_line(aes(x = Seg_EP, y = Seg_C_EP), color = 'yellow', size = 1.5)+
  scale_color_manual(name = "Log(C)~EP-Q\n Breakpoint Analysis", values = "yellow")+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free_y')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  scale_x_reverse()+
  geom_rect(data = df_Seg.EP%>%distinct(df_Seg.EP$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

save(p, file = 'Processed_Data/p1.Rdata')

# I like the idea of using the EP plots to help ID which
# sites are complex. For example, site 22 was deemed complex
# with log-log plots, but in log-EP space it looks similar to
# the other mobilization sites. 























#----------------------------------------#
####~~~  Reclassifying CQ type ~~~#### 
###~~~  using the EP-Q plots  ~~~### 
#----------------------------------------#

# rerun the df_Seg.2 workflow bt=ut changing the complex_sites vector:
temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.datalayers$Name)%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))
l.temp<-temp%>%
  split(., .$Name)
l.lm.CQ_slopes<-lapply(l.temp, \(i) lm(log_C~log_Q, data=i))
coef<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,1] ))), 'site_no')
pvals<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,4] ))), 'site_no')%>%
  rename(I.pval = 2, S.pval = 3)
m<-left_join(coef,pvals,by='site_no')
m<-mutate(m, Type = ifelse(S.pval>0.05, 'Stationary', ifelse(log_Q>0, 'Mobilization', 'Dilution')))
temp<-left_join(temp,m%>%select(site_no, Type),by=c('Name'='site_no'))
df_Seg.2<-filter(df_Seg, site %in% temp$Name)%>%
  left_join(.,m%>%select(site_no, Type),by=c('site'='site_no')) 
df_Seg.2<-df_Seg.2%>%
  mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))
hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)
complex_sites<-unique(df_Seg.2$site)[c(3,16,21,37)] 
df_Seg.2<-mutate(df_Seg.2, Type = ifelse(site %in% complex_sites, 'Complex', Type))
df_Seg.2<-left_join(df_Seg.2, df.datalayers%>%select(Name, USGS.LU.Adjusted, CAFO_count), by = c('site'='Name'))
df_Seg.2$CAFO_count[df_Seg.2$CAFO_count==0]<-NA

# plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

#







































































#### Delinating Watersheds (not run) ####

# # delinate the 137 NWIS sites:
# 
# l.SS_WS.NWIS<-lapply(seq_along(df.NWIS.TP_site_metadata$site_no), \(i) Delineate(df.NWIS.TP_site_metadata$dec_long_va[i], df.NWIS.TP_site_metadata$dec_lat_va[i]))
# 
# names(l.SS_WS.NWIS)<-df.NWIS.TP_site_metadata$site_no
# 
# # save(l.SS_WS.NWIS, file = 'C:/PhD/CQ/Downloaded_Data/l.SS_WS.NWIS.Rdata')
# 
# df.sf.NWIS<-fun.l.SS_WS.to.sfdf(l.SS_WS.NWIS)
# 
# # save(df.sf.NWIS, file = 'C:/PhD/CQ/Processed_Data/df.sf.NWIS.Rdata')
# 
# load('C:/PhD/CQ/Processed_Data/df.sf.NWIS.Rdata')
# 
# # calculate DA using vect:
# 
# vect.NWIS<-vect(df.sf.NWIS)
# 
# vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")
# 
# # add DA in mi2 to df.sf:
# 
# df.sf.NWIS$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102
# 
# # add NWIS tabulated area to df sf: first download the metadata for the sites using readNWISsite, then merge drain_area_va column to sfdf:
# 
# temp<-readNWISsite(df.sf.NWIS$Name)
# 
# df.sf.NWIS<-left_join(df.sf.NWIS, temp[,c(2,30)], by = c('Name'='site_no'))
# 
# # calcuate the percent error of the delination:
# 
# df.sf.NWIS$Delination_Error<-(df.sf.NWIS$area_sqmi-df.sf.NWIS$drain_area_va)/df.sf.NWIS$drain_area_va
# 
# # I want to look at the delimaitons that are on the cusp of not working:
# 
# temp<-df.sf.NWIS[,c(1,107)]%>%
#   filter(abs(Delination_Error)>.02)
# 
# # mapview(temp)
# 
# # after playing with the numbers, 2% seems like a good threshold for non-delineations:
# # subset df.sf to those sites that are considered to delineate correctly:
# 
# df.sf.NWIS.keep<-df.sf.NWIS%>%filter(abs(Delination_Error)<=.02)
# 
# # note I suspect that snapping the lat long to nhd would improve the delineation success notgoing to do that now









































############################################################
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
####~~~~~~~~~~~~~~~~~Watershed Attributes (NOT RUN)~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~(for correlation, eventually)~~~~~~~~~~~####
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
############################################################



####~~~~ Data Layers ~~~~####

#### CDL: ####

# note when I ran ths CDL code block I used the 103 sites thatwere apartofthe set, this is prior to filtering based onthe CDLdate... so I am just going to filter the resulting CDL df at the end. Just note that what is saved is of the 103 sites, but that the code doesnt reflect this, i.e. I just ran the filter anddidnt keep the code because df.NWIS.TP.keep would be down to 56 if I did it right 

# download:

rast.NWIS.CDL.2020 <- GetCDLData(aoi = df.sf.NWIS, year = "2020", type = "b", tol_time = 1000)
rast.NWIS.CDL.2008 <- GetCDLData(aoi = df.sf.NWIS, year = "2008", type = "b", tol_time = 1000)

# convert to rast:

rast.NWIS.CDL.2020<-rast(rast.NWIS.CDL.2020)
rast.NWIS.CDL.2008<-rast(rast.NWIS.CDL.2008)

# check CRS of both years:

crs(rast.NWIS.CDL.2020)
crs(rast.NWIS.CDL.2008)

# plot

# plot(rast.NWIS.CDL.2020)
# plot(rast.NWIS.CDL.2008)

# reproject to sample watershed vector data to match raster data:

vect.NWIS<-vect(df.sf.NWIS) # need to first ininalize vect since sometimes reading in rdata file (terra issue)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(rast.NWIS.CDL.2008))

# extract frequency tables for each sample watershed

l.NWIS.CDL <- terra::extract(rast.NWIS.CDL.2020, vect.NWIS.proj, table,ID=FALSE)[[1]]
l.NWIS.CDL.2008 <- terra::extract(rast.NWIS.CDL.2008, vect.NWIS.proj, table,ID=FALSE)[[1]]

# save:

# save(l.NWIS.CDL, file='Processed_Data/l.NWIS.CDL.Rdata')
# save(l.NWIS.CDL.2008, file='Processed_Data/l.NWIS.CDL.2008.Rdata')

load('Processed_Data/l.NWIS.CDL.Rdata')
load('Processed_Data/l.NWIS.CDL.2008.Rdata')

# convert resulting list of tables to list of dfs

l.NWIS.CDL<-lapply(l.NWIS.CDL, as.data.frame)
l.NWIS.CDL.2008<-lapply(l.NWIS.CDL.2008, as.data.frame)

# the next step is to join the CDL key with these dfs, but first:

# aggregate CDL: to do this:
# looking at the CDL legend (linkdata), I determined the following:

Ag<-c(1:6,10:14,21:39,41:61,66:72,74:77,204:214,216:227,229:250,254)
Pasture<-c(176)
Forest<-c(63,141:143)
Developed<-c(82,121:124)
Water<-c(83,111)
Wetlands_all<-c(87,190,195)
Other<-c(64,65,88,112,131,152) # Shrub<-c(64,152) Barren<-c()

l <- tibble::lst(Ag,Pasture,Forest,Developed,Water,Wetlands_all,Other)

reclass_CDL<-data.frame(lapply(l, `length<-`, max(lengths(l))))%>%
  pivot_longer(cols = everything(), values_to = 'MasterCat',names_to = 'Crop')%>%
  drop_na(MasterCat)

# left join each df in the list to the CDL legend key, as well as calcuate the pland:
  
l.NWIS.CDL<-lapply(l.NWIS.CDL, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., reclass_CDL, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes
l.NWIS.CDL.2008<-lapply(l.NWIS.CDL.2008, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., reclass_CDL, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes

# set names of list:

names(l.NWIS.CDL)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter
names(l.NWIS.CDL.2008)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter

# combine list into single df:

df.NWIS.CDL<-bind_rows(l.NWIS.CDL, .id = 'Name')
df.NWIS.CDL.2008<-bind_rows(l.NWIS.CDL.2008, .id = 'Name')

# remove potential for one of the MasterCats in the CDL to be empty, which is messing with the pivot_wider below:

df.NWIS.CDL<-filter(df.NWIS.CDL, Crop != '')
df.NWIS.CDL.2008<-filter(df.NWIS.CDL.2008, Crop != '')

# pivot wider:

# if using the CDL linkdata, use this:

# df.NWIS.CDL<- pivot_wider(df.NWIS.CDL[,-2], names_from = Crop, values_from = Freq)
# df.NWIS.CDL.2008<- pivot_wider(df.NWIS.CDL.2008[,-2], names_from = Crop, values_from = Freq)

# if using the reclass_CDL data, use this:

df.NWIS.CDL<-df.NWIS.CDL[,-2]%>%
  group_by(Name, Crop)%>%
  summarise(Freq=sum(Freq, na.rm = T))%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

df.NWIS.CDL.2008<-df.NWIS.CDL.2008[,-2]%>%
  group_by(Name, Crop)%>%
  summarise(Freq=sum(Freq, na.rm = T))%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

# check to see if add up to 100%:

sort(rowSums(df.NWIS.CDL[,-1], na.rm = T))
sort(rowSums(df.NWIS.CDL.2008[,-1], na.rm = T))

# one of the CDL 2020 is low:

which.min(rowSums(df.NWIS.CDL[,-1], na.rm = T))

# this is a long island site, so doesnt matter

# other than that site looks good. 

# Done with CDL

# actually I want to come up with a workflow to get all years of CDL for the watersheds,
# butthat isgoing totake a long timeto run...

# but first I'm going to look at the differences in the CDL between 2008 and 2020


























#### NLCD ####

# the CDL and NLCD differ in that the CDL combined Pasture and grassland while the nLCD does not
# thus, even though the CDL can be agrgated to et close to the NLCD, I want to also have the nLCD
# to the difference.

# note: the 53 sites post filtering based onthe CDLdate are used in this workflow

# I am going to run this workflow for NLCD 2019 and NLCD 2001 to see if any sites had major changes in land use:

# download: to do this:
# NLCD is not working when trying to download based on the entire polygon df, so going to use lapply to download individually:

l.rast.NWIS.NLCD.2019 <- lapply(seq_along(df.sf.NWIS.keep$Name), \(i) get_nlcd(template = st_cast(df.sf.NWIS.keep, "MULTIPOLYGON")[i,], label = as.character(i), year = 2019))
l.rast.NWIS.NLCD.2001 <- lapply(seq_along(df.sf.NWIS.keep$Name), \(i) get_nlcd(template = st_cast(df.sf.NWIS.keep, "MULTIPOLYGON")[i,], label = as.character(i), year = 2001))

# save(l.rast.NWIS.NLCD.2019, file='Downloaded_Data/l.rast.NWIS.NLCD.2019.Rdata')
# save(l.rast.NWIS.NLCD.2001, file='Downloaded_Data/l.rast.NWIS.NLCD.2001.Rdata')

# convert to SpatRasters:

l.rast.NWIS.NLCD.2019<-lapply(l.rast.NWIS.NLCD.2019, rast)
l.rast.NWIS.NLCD.2001<-lapply(l.rast.NWIS.NLCD.2001, rast)

# plot

# plot(rast.NWIS.NLCD.2019)

# see if 2001 and 2019 crs are the same:

crs(l.rast.NWIS.NLCD.2019[[2]])
crs(l.rast.NWIS.NLCD.2001[[2]])

# reproject to sample watershed vector data to match raster data:

vect.NWIS<-vect(df.sf.NWIS.keep) # need to first ininalize vect since sometimes reading in rdata file (terra issue)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(l.rast.NWIS.NLCD.2019[[1]]))

# extract frequency tables for each sample watershed

system.time({l.NWIS.NLCD.2019 <- lapply(seq_along(l.rast.NWIS.NLCD.2019), \(i) terra::extract(l.rast.NWIS.NLCD.2019[[i]], vect.NWIS.proj[i], ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2)))})
system.time({l.NWIS.NLCD.2001 <- lapply(seq_along(l.rast.NWIS.NLCD.2001), \(i) terra::extract(l.rast.NWIS.NLCD.2001[[i]], vect.NWIS.proj[i], ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2)))})

# save(l.NWIS.NLCD.2019, file = 'Processed_Data/l.NWIS.NLCD.2019.Rdata')
# save(l.NWIS.NLCD.2001, file = 'Processed_Data/l.NWIS.NLCD.2001.Rdata')

load("Processed_Data/l.NWIS.NLCD.2019.Rdata")
load("Processed_Data/l.NWIS.NLCD.2001.Rdata")

# reclassify: to do this:

# adjust the NLCD reclassify legend to match the CDL reclassify df made above:

legend.NWIS<-legend; legend.NWIS$Class3[c(13,17)]<-'Pasture';legend.NWIS$Class3[14]<-'Other';legend.NWIS$Class3[1]<-'Water';legend.NWIS$Class3[c(19,20)]<-'Wetlands_all'

sort(unique(legend.NWIS$Class3))==sort(unique(reclass_CDL$Crop))

# looks good

# reclassify the NLCD using this new legend and clean up the dataframe from the next step

l.NWIS.NLCD.2019<-lapply(l.NWIS.NLCD.2019, \(i) left_join(as.data.frame(i), legend.NWIS%>%select(Class, Class3), by = 'Class')%>%mutate(Class = Class3)%>%select(-Class3))
l.NWIS.NLCD.2001<-lapply(l.NWIS.NLCD.2001, \(i) left_join(as.data.frame(i), legend.NWIS%>%select(Class, Class3), by = 'Class')%>%mutate(Class = Class3)%>%select(-Class3))

# pivot_wider the df in the lists and add a Name column for the site:

l.NWIS.NLCD.2019<-lapply(seq_along(l.NWIS.NLCD.2019), \(i) l.NWIS.NLCD.2019[[i]]%>%group_by(Class)%>%summarise(Freq = sum(Freq))%>%pivot_wider(names_from = Class, values_from = Freq)%>%mutate(Name = df.sf.NWIS.keep$Name[i], .before = 1)%>%as.data.frame(.))
l.NWIS.NLCD.2001<-lapply(seq_along(l.NWIS.NLCD.2001), \(i) l.NWIS.NLCD.2001[[i]]%>%group_by(Class)%>%summarise(Freq = sum(Freq))%>%pivot_wider(names_from = Class, values_from = Freq)%>%mutate(Name = df.sf.NWIS.keep$Name[i], .before = 1)%>%as.data.frame(.))

# bind the lists into a single dataframe
# note some of the sites have different length dtaframes because they didnt have all the same number of NLCD classes. When binding rows this will give a dataframe of the maximum length and put NAs for sites where there wasn't a column: 

df.NWIS.NLCD.2019<-bind_rows(l.NWIS.NLCD.2019)
df.NWIS.NLCD.2001<-bind_rows(l.NWIS.NLCD.2001)

# create a list of these two dataframes:

l.NWIS.NLCD<-list(df.NWIS.NLCD.2019,df.NWIS.NLCD.2001)%>%purrr::set_names(c('2019','2001'))

# save thisprocessed list:

# save(l.NWIS.NLCD, file='Processed_Data/l.NWIS.NLCD.Rdata')

load('Processed_Data/l.NWIS.NLCD.Rdata')

#

























#### NED ####

# download:

DEM.NWIS<-get_ned(df.sf.NWIS, label = '2') # already SpatRaster!

writeRaster(DEM.NWIS, file='Downloaded_Data/DEM.NWIS.tif', overwrite=TRUE)

# plot(DEM.NWIS)

# extract elevation metrics over each sample watershed: to do this:
# build a function with multiple functions:

f <- function(x, na.rm = T) {
  c(mean=mean(x, na.rm = na.rm),
    range=max(x, na.rm = na.rm)-min(x, na.rm = na.rm),
    sd=sd(x, na.rm = na.rm)
  )
}

# reproject NWIS basins to DEM crs:

vect.NWIS<-vect(df.sf.NWIS)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(DEM.NWIS))

# extract the metrics over each watershed using the function above:

df.NWIS.DEM <- as.data.frame(terra::extract(DEM.NWIS, vect.NWIS.proj, f))

# set the names of the df:

names(df.NWIS.DEM)<-c('Name', 'Elev_Avg', 'Elev_Range', 'Elev_SD')

# set the names of the sites:

df.NWIS.DEM$Name<-df.sf.NWIS$Name

# finally save the df:

# save(df.NWIS.DEM, file= 'Processed_Data/df.NWIS.DEM.Rdata')

load('Processed_Data/df.NWIS.DEM.Rdata')

# done with DEM

























#### Climate (not run) ####

# steps:

# Loop through the three metrics and metric specific operations.
# for each metric, perform:

# 1) download all gridded climate data for a metric (e.g.: precip)
# 2) extract daily mean of the metric for each subbasin into df
# 3) reformat (pivot longer),convert the date, and calcualte annual statistics based on the metric (e.g. precip is annual totals)

# the definition of:
# pr = total annual precio
# tmmn = minimum of all annual daily minimum temperature 
# tmmx = maximum of all annual daily maximum temperature

vars<-c('pr', 'tmmn', 'tmmx')

vars_funs<-c(`sum`, `min`, `max`)

# set up a df for appending into each loop:

df.NWIS.Climate<-data.frame(ID = 1:dim(df.NWIS.CDL)[1])

i<-1

# loop through the climate variables in vars:

for (i in seq_along(vars)){
  
  # download a Spatraster for the entire bounding box of all the watersheds of the climate variable data for the loops iteration
  # each layer of the spatraster is a day for the climate variable in thedaterange provided inthe function call:
  
  climate<-getGridMET(vect.NWIS, varname = vars[i], startDate = "2016-01-01", endDate = "2019-12-31")[[1]]
  
  # reproject the watersheds to the crs of the climate variables raster:
  
  vect.NWIS.proj<-terra::project(vect.NWIS, crs(climate))
  
  # extract the daily mean of the climate variables data over the watersheds:
  # thisreturns a dataframe thefirstcolumn is the site (justnumbered 1-56 right now) and the other columns are the daily values of the climate vairable (taking the mean ofthe raster cells over the watershed): 
  
  climate <- terra::extract(climate, vect.NWIS.proj, mean)
  
  # reformat, convert dates, and calcualte annual stats:
  
  # pivot Longer: each of the columns for the daily values of the climate variable starts withthevariable name and an underscore, so using that pattern to make a long dataframe 
  # format the dates: the date is contained withinthe column that was priviously pivoted into a single column. so reformating this column and turing into a date:
  # grouping by the site and year, summarize the daily values into an annualvalue based onthe operation for the climatevariable. for pricep, sum is used, for min temp min is used, and formax temp max is used!:
  # rename the reulting year(Date) column to a nicer name:
  # renaming the values ofthe year column: add the climate variable to the year,this pairs withthe next andfinal step:
  # pivot wider to get unique columns for the years of the cimlatevariable (here its 2016-2019):
  
  climate1<-climate%>%
    pivot_longer(cols= starts_with(paste0(vars[i],"_")), names_to = 'Date',values_to = "Value")%>%
    mutate(Date=as.Date(str_replace(Date, paste0(vars[i],"_"), "")))%>%
    group_by(year(Date),ID)%>%
    summarize(Annual = vars_funs[[i]](Value))%>%
    rename(Year=1)%>%
    mutate(Year=paste0(vars[i],Year))%>%
    pivot_wider(names_from = Year, values_from = Annual)
  
  # merge the climate variables data to a common dataframe (initalized prior to loop).
  # New columns are added for each site:
  
  df.NWIS.Climate<-left_join(df.NWIS.Climate, climate1, by = 'ID')
  
}

# rename the sites:

df.NWIS.Climate[,1]<-df.NWIS.DEM$Name

# rename the sites:

names(df.NWIS.Climate)[1]<-'Name'

# save df:

save(df.NWIS.Climate, file = 'Processed_Data/df.NWIS.Climate.Rdata')

#























##### Determine if CAFO is in watershed ####

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

# determine how many CAFOs intersect each watershed:

df.sf.NWIS.keep$CAFO_count <- lengths(st_intersects(df.sf.NWIS.keep, CAFOs))

# make a map colored by CAFO count:

# mapview(df.sf.NWIS.keep, zcol = 'CAFO_count')

#































































####~~~~ Sufficentcy of Gauges 2 (NOT RUN) ~~~~####

# below is the code used when I was looking at the sufficency of gauges 2 sites
# regardless of the year of TP data. But filtring on 2001 above, I dont really care how many years of 
# data they have (right now at least...)

# # import gages 2 and filter to NYS TP CQ sites:
# 
# l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm_EDITED_again.xlsx")
# df.G2<-reduce(l.G2, full_join, by = "STAID")
# df.G2<-select(df.G2, STAID | where(is.numeric))
# df.G2<-df.G2%>%filter(STAID %in% TP_sites$site_no)
# 
# # add the n_sample_rank: to do this:
# # rerun the temp dataframe (since I keep writing over it in different sections:
# 
# temp<-df.NWIS.TP%>%group_by(site_no)%>%
#   summarise(n=n())%>%
#   arrange(desc(n))%>%
#   mutate(n_sample_rank=rank(-n, ties.method='first'))
# 
# # left join to add n_sample rank column:
# 
# df.G2<-left_join(df.G2, temp, by = c('STAID'='site_no'))
# 
# # I want look at the number of years of data for each of these:
# 
# G2_n_years<-filter(df.NWIS.TP_CQ, site_no %in% df.G2$STAID)%>%
#   group_by(site_no)%>%
#   summarize(min_date=min(year(sample_dt)),
#             max_date=max(year(sample_dt)),
#             n_years=max(year(sample_dt))-min(year(sample_dt)))%>%
#   arrange(desc(n_years))
# 
# # I want to look at the time series of the sites with little years of data: to do this:
# 
# # filter to sits with lessthan 10 years of data:
# 
# x<-filter(G2_n_years, n_years<10)
# 
# # plot:
# 
# p<-filter(df.NWIS.TP_CQ, site_no %in% x$site_no)%>%
#   # filter(., n_sample_rank==196)%>%
#   ggplot(., aes(x = as.Date(sample_dt), y = result_va))+
#   geom_point()+
#   facet_wrap('n_sample_rank', scales ='free')
# 
# # from this plot the following sites I say have good data even though they have less than 10 years of samples:
# 
# below_10_keep<-c(46,70,72,118,119,120,122,125,126,134)
# 
# # now come up with a final list of sites: to do this:
# 
# # merge the n_samplerank and n_years dfs and filter based on n_years and the numbers of n_sample_ranks in below_10_keep:
# 
# G2_keep<-left_join(G2_n_years, temp, by = 'site_no')%>%
#   filter(., n_years>=10 | n_sample_rank %in% below_10_keep)
# 
# # make a plot:
# 
# p1<-filter(df.NWIS.TP_CQ, site_no %in% G2_keep$site_no)%>%
#   # filter(., n_sample_rank==196)%>%
#   ggplot(., aes(x = sample_dt, y = result_va))+
#   facet_wrap('n_sample_rank', scales ='free')+
#   geom_point()
# 
# # now filter the gauges2 to these sites:
# 
# df.G2<-filter(df.G2, STAID %in% G2_keep$site_no)























































































#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                     ####~~~~ MLR regular ~~~~####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# set up list of dataframes for different selection methods used below:

# create a list of dfs, one for each cQ parameter, with the values of the spearman correlations for each predictor:

l.cor.MLR<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(desc(abs(Spearman_Correlation)))%>%  
           # slice_head(n = 7)%>%
           # bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           # distinct(term, .keep_all = T)%>%
           # filter(!between(Spearman_Correlation, -0.25,.25))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           # mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           # filter(sig_0.05 == 'sig')%>%
           as.data.frame())


x<-l.cor.MLR[[1]]

# create a list of dfs from df.OLS.Sens, which contains the values of the predictors, subseting using the names in each l.cor[[i]]$term:

l.cor.MLR.full<-lapply(1:5, \(i) df.OLS_Sens%>%
                    select(i+1, l.cor.MLR[[i]]$term)%>%
                    as.data.frame()%>%
                      rename(term = 1))

names(l.cor.MLR.full)<-names(df.OLS_Sens)[2:6]

y<-l.cor.MLR.full[[1]]

y.round<-round(y, 2)

# create another list but keep the site names for use at end of loop in site outliers:

l.cor.MLR.full.w_name<-lapply(1:5, \(i) df.OLS_Sens%>%
                         select(1, i+1, l.cor.MLR[[i]]$term)%>%
                         as.data.frame())

names(l.cor.MLR.full.w_name)<-names(df.OLS_Sens)[2:6]

z<-l.cor.MLR.full.w_name[[1]]

# Trying no feature selection:

# loop through the 5 CQ parameters:

i<-1

for (i in 1:5){
  
  df<-l.cor.MLR.full[[i]]
  
  library(caret)
  
  set.seed(123)
  
  # Set up repeated k-fold cross-validation
  
  train.control <- trainControl(method = "cv", number = 10)
  
  # Train the model
  
  step.model <- train(term ~., data = df,
                      method = "leapForward", 
                      tuneGrid = data.frame(nvmax = 1:5),
                      trControl = train.control
  )
  

  # look at model results:
  
  step.model$results
  step.model$bestTune
  summary(step.model$finalModel)
  coef(step.model$finalModel, 4)
  
  # get df of just the predictors in this model:
  
  pred<-names(coef(step.model$finalModel, 4))[-1]
  
  df.2<-df%>%select(term, pred)
  
  # make lm:
  
  m<-lm(term~., data=df.2)
  
  summary(m)
  
  # look at univarite
  
  plots:



#

# loop:

i<-1

for (i in 1:5){
  
  # when i<-1, OLS.I
  
  # regardless of the method to pick the predictors for the MLR model, no more than 4 wil be included as per the run p < n/10
  
  #################################################
  # Spearman correlations::
  #################################################

  # build a lm using the top 4 correlates from the spearman correlaiton analysis:
  
  data.temp<-l.cor.MLR.full[[i]][,1:5]
  m.spear<-lm(term~., data = data.temp)
  
  summary(m.spear) # Adj R2 = 0.2121, pval = 0.01154, none of the terms signficant
  
  plot(m.spear) # 12, 22
  
  terms.m.spear<-names(m.spear$model)[-1]
  
  #################################################
  # Random Forest: 
  # changes everytime slightly even if seed is same, 
  # I think this is called unstable
  #################################################
  
  set.seed(1)
  
  # build RF model:
  x<-randomForest(term~., data = l.cor.MLR.full[[i]], na.action=na.omit)
  # arrange by RF variable importance 
  y<-tibble::rownames_to_column(as.data.frame(x$importance), "Attribute")%>%
    arrange(desc(IncNodePurity))
  # build lm of top 4:
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, y$Attribute[1:4])
  m.rf<-lm(term~., data = data.temp)
  summary(m.rf) # adj R2 = 0.2494, pval = 0.005165, only RIP100_FOREST signficant
  plot(m.rf) # 22, 12
  
  terms.m.rf<-names(m.rf$model)[-1]
  
  #################################################
  # stepAIC:
  #################################################
  
  #### v1: full model:
  
  # stepAIC workflow:
  fit0 <- lm(term~1,data=l.cor.MLR.full[[i]])
  up.model <- paste("~", paste(colnames(l.cor.MLR.full[[i]][,-1]), collapse=" + "))
  m.step.all <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.all) # 
  
  #### v2: reduce predictor set prior to building full model:
  
  # scale variables because some have small ranges:
  
  data.temp<-l.cor.MLR.full[[i]][,-1]
  
  ### v2.1: check which predictors have the highest coefficent of variability among the 40 observations:
  
  v.mean<-sapply(l.cor.MLR.full[[i]][,-1], mean, na.rm=T)
  v.sd<-sapply(l.cor.MLR.full[[i]][,-1], sd, na.rm=T)
  v.cv<-v.mean/v.sd
  df.cv<-tibble::rownames_to_column(as.data.frame(v.cv), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.cv$Attribute[1:39]
  # use stepAIC to build a 'full' model using the top n-1:
  fit0 <- lm(term~1,data=l.cor.MLR.full[[i]])
  up.model <- paste("~", paste(colnames(l.cor.MLR.full[[i]]%>%select(df.cv$Attribute[1:39])), collapse=" + "))
  m.step.39.cv <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.cv) # adj R2 = 0.7047, pval = 7.4e-08
  # this is pretty good but I still want p<n/10 
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.cv$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, var.temp)
  m.leaps.step.cv.39to8to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.cv.39to8to4)
  rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1])
  m.leaps.step.cv.39to8to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.cv.39to8to4) # adj r2 = 0.5395, pval = 2.2e-06
  # plot(m.leaps.step.cv.39to8to4) # 39, 21, 23 (21 in cooks)
  
  terms.m.step.cv<-names(m.leaps.step.cv.39to8to4$model)[-1]
  
  ### v2.2: check which predictors have the highest variability among the 40 observations:
  
  v.var<-sapply(l.cor.MLR.full[[i]][,-1], var, na.rm=T)
  df.var<-tibble::rownames_to_column(as.data.frame(v.var), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.var$Attribute[1:39]
  # build the full model using the top n-1:
  fit0 <- lm(term~1,data=l.cor.MLR.full[[i]])
  up.model <- paste("~", paste(colnames(l.cor.MLR.full[[i]]%>%select(df.var$Attribute[1:39])), collapse=" + "))
  m.step.39.var <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.var) # adj R2 = 0.9259, pval = 1.01e-08
  # this is amazing but we shouldn't use a model with 22 predictors
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.var$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, var.temp)
  m.leaps.step.var.39to22to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.var.39to22to4)
  rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1])
  m.leaps.step.var.39to22to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.var.39to22to4) # adj r2 = 0.467, pval = 2.5e-05
  # plot(m.leaps.step.var.39to22to4) # 23, 25, 10, 21 (21 in cooks)
  
  terms.m.step.var<-names(m.leaps.step.var.39to22to4$model)[-1]
  
  #################################################
  # LEAPS:
  # (straight up with all 306 varables):
  #################################################
  
  # leaps: (doesnt run, too big)
  
  # fit.reg<-leaps::regsubsets(term~.,data=l.cor.MLR.full[[i]],nbest=4, really.big=T)
  
  #################################################
  # LASSO:
  # I dont know the difference between the two lasso methods I have below
  #################################################
  
  #### lasso v1: 
  
  # Define predictor and response variables
  y <- l.cor.MLR.full[[i]]$term
  x <- data.matrix(l.cor.MLR.full[[i]][, -1])
  #fit lasso regression model using k-fold cross-validation
  cv_model <- cv.glmnet(x, y, alpha = 1)
  best_lambda <- cv_model$lambda.min
  #display optimal lambda value
  best_lambda
  #view plot of test MSE's vs. lambda values
  plot(cv_model)
  #view coefficients of best model
  best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  coef(best_model)
  # get non zero coeffs:
  Shat.v1 <- rownames(coef(best_model))[which(coef(best_model) != 0)]
  # since only 3 given, build lm with these!:
  data.temp<-select(l.cor.MLR.full[[i]], term, Shat.v1[-1])
  m.lasso.v1<-lm(term~., data=data.temp)
  summary(m.lasso.v1) # adj r2 = 0.376, pval = 0.00016
  plot(m.lasso.v1) # 23, 15, 21, 28 (28 in cooks)
  
  terms.m.lasso.v1<-names(m.lasso.v1$model)[-1]
  
  #### lasso v2:
  
  #Define predictor and response variables
  ytrain <- l.cor.MLR.full[[i]]$term
  xtrain <- data.matrix(l.cor.MLR.full[[i]][, -1])
  # fit model:
  fit.lasso.glmnet <-glmnet(x=xtrain,y=ytrain,alpha=1) 
  # fit cv (?)
  cv.lasso.glmnet <-cv.glmnet(x=xtrain,y=ytrain,alpha=1) 
  plot(cv.lasso.glmnet)
  # save beta coefficents (?)
  beta.lasso <- coef(fit.lasso.glmnet, s = cv.lasso.glmnet$lambda.min)
  beta.lasso
  # get non zero coeffs:
  Shat.v2 <- rownames(beta.lasso)[which(beta.lasso != 0)]
  # since only 4 given, build lm with these!:
  data.temp<-select(l.cor.MLR.full[[i]], term, Shat.v2[-1])
  m.lasso.v2<-lm(term~., data=data.temp)
  summary(m.lasso.v2) # adj r2 = 0.43, pval = 7.8e-05
  plot(m.lasso.v2) # 21, 39, 23, 28 (28 in cooks)
  
  terms.m.lasso.v2<-names(m.lasso.v2$model)[-1]

  #################################################
  # Correlations:
  #################################################

  # calculate correlation matrix
  correlationMatrix <- cor(l.cor.MLR.full[[i]][,-1])
  # summarize the correlation matrix
  print(correlationMatrix)
  # ggcorrplot(correlationMatrix)
  # find attributes that are highly corrected (ideally >0.75)
  highlyCorrelated <- caret::findCorrelation(correlationMatrix, cutoff=0.5)
  # print indexes of highly correlated attributes
  print(highlyCorrelated)
  unique(highlyCorrelated)
  # remove highly correlated variables:
  data.temp<-l.cor.MLR.full[[i]][,-c(highlyCorrelated+1)]
  # stepAIC workflow:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp[,-1]), collapse=" + "))
  m.step.cor <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.cor) # adj r2=.4423,pval=0.0006
  # leaps workflow:
  var.temp<-strsplit(as.character(m.step.cor$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, var.temp)
  m.leaps.step.cor<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.cor)
  rownames(as.data.frame(coef(m.leaps.step.cor, 4)))[-1]
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.cor, 4)))[-1])
  m.leaps.step.cor<-lm(term~., data=data.temp)
  summary(m.leaps.step.cor) # adj r2 = 0.2911, pval = 0.0027
  plot(m.leaps.step.cor) # 23, 1, 21 
  
  terms.m.cor<-names(m.leaps.step.cor$model)[-1]
  
  #################################################
  # Feature Importance
  #################################################
  
  control <- caret::trainControl(method="repeatedcv", number=10, repeats=3)
  # train the model
  model <- caret::train(term~., data=l.cor.MLR.full[[i]], method="leapSeq", preProcess="scale", trControl=control)
  # estimate variable importance
  importance <- varImp(model, scale=FALSE)
  # summarize importance
  print(importance)
  # plot importance
  plot(importance)
  # get top 4 most important:
  x<-tibble::rownames_to_column(importance$importance, "Attribute")%>%arrange(desc(Overall))
  # build lm with top 4:
  data.temp<-select(l.cor.MLR.full[[i]], term, x$Attribute[1:4])
  m.feat_import.4<-lm(term~., data= data.temp)
  summary(m.feat_import.4) # adj r2 = 0.2229, pval=0.01148
  plot(m.feat_import.4) # 23, 21, 11
  
  terms.m.feat_import<-names(m.feat_import.4$model)[-1]
  
  #################################################
  # Feature Selection:
  #################################################
  
  # define the control using a random forest selection function
  control <- caret::rfeControl(functions=rfFuncs, method="cv", number=10)
  # run the RFE algorithm
  results <- caret::rfe(l.cor.MLR.full[[i]][,-1], l.cor.MLR.full[[i]][,1], sizes=c(1:4), rfeControl=control)
  # summarize the results
  print(results)
  # list the chosen features
  predictors(results)
  # plot the results
  plot(results, type=c("g", "o"))
  # build lm with top 4:
  x<-predictors(results)[1:4]
  # build lm with top 4:
  data.temp<-select(l.cor.MLR.full[[i]], term, x)
  m.feat_select.4<-lm(term~., data= data.temp)
  summary(m.feat_select.4) # adj r2 = 0.244, pval=0.0074
  plot(m.feat_select.4) # 23, 21, 27
  
  terms.m.feat_select<-names(m.feat_select.4$model)[-1]
  
  #################################################
  # PCA:
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  #################################################
  
  data.temp<-l.cor.MLR.full[[i]][,-1]
  res.pca <- PCA(data.temp, graph = FALSE)
  print(res.pca)
  eig.val <- get_eigenvalue(res.pca)
  eig.val # sum(eig.val[,1]) # The sum of all the eigenvalues give a total variance of p (305 here).
  fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
  var <- get_pca_var(res.pca)
  var
  fviz_pca_var(res.pca, col.var = "black")
  library("corrplot")
  corrplot(var$cos2, is.corr=FALSE)
  fviz_cos2(res.pca, choice = "var", axes = 1:2)
  fviz_pca_var(res.pca, col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) # Avoid text overlapping)
  fviz_pca_var(res.pca, alpha.var = "cos2")
  corrplot(var$contrib, is.corr=FALSE) 
  # Contributions of variables to PC1
  fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
  # Contributions of variables to PC2
  fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
  # The total contribution to PC1 and PC2 is obtained with the following R code:
  fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)
  # The most important (or, contributing) variables can be highlighted on the correlation plot as follow:
  fviz_pca_var(res.pca, col.var = "contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
  # Create a random continuous variable of length 10
  set.seed(123)
  my.cont.var <- rnorm(305)
  # Color variables by the continuous variable
  fviz_pca_var(res.pca, col.var = my.cont.var,gradient.cols = c("blue", "yellow", "red"),legend.title = "Cont.Var")
  # Create a grouping variable using kmeans
  # Create 3 groups of variables (centers = 3)
  set.seed(123)
  res.km <- kmeans(var$coord, centers = 3, nstart = 25)
  grp <- as.factor(res.km$cluster)
  # Color variables by groups
  fviz_pca_var(res.pca, col.var = grp, palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),legend.title = "Cluster")
  # graph of individuals
  fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE) # Avoid text overlapping (slow if many points))
  # You can also change the point size according the cos2 of the corresponding individuals:
  fviz_pca_ind(res.pca, pointsize = "cos2", pointshape = 21, fill = "#E7B800",repel = TRUE) # Avoid text overlapping (slow if many points)
  # Total contribution on PC1 and PC2
  fviz_contrib(res.pca, choice = "ind", axes = 1:2)
  # biplot:
  fviz_pca_biplot(res.pca, repel = TRUE,col.var = "#2E9FDF", col.ind = "#696969"  )
  
  # Im just going to stop with PCA for now because I dont know how to interpret results:
  
  #################################################
  # Looking at regression outliers:
  #################################################
  
  # The following observations (1-40) came up in the plot(lm) calls throughout thecode above for i=1 (OLS.I):
  # n times (second number)
  
  x<-c(1, 1, 10, 1 , 11, 2, 15, 1, 21, 9, 23, 9, 25, 1, 27, 1, 28, 2, 37, 2, 39, 2)
  x<-as.data.frame(do.call(cbind, split(x, c("Site", "n")))) %>% separate_rows("n", sep = ",")%>%arrange(n)
  
  # more than any others, sites 21 and 23 kept coming up as outliers
  # look at these sites. to do this:
  
  look_at_sites<-l.cor.MLR.full.w_name[[i]][c(21,23),]$Name 
  
  map.1<-df.sf.NWIS%>%filter(Name %in% look_at_sites)
    
  mapview(map.1)
  
  # if I look at the facet plot for CQ type and land use, nothing stands out 
  # about these problem sites. I mean, two are ag and two are undeveloped, but I 
  # dont trust the land use classifications. Also two classes out of 4 isnt that impressive. 
  
  #################################################
  # Looking at frequcny of terms across the models:
  #################################################
  
  x<-grep("terms",names(.GlobalEnv),value=TRUE)
  x<-do.call("list",mget(Pattern1))
  x<-unlist(x)
  x<-as.data.frame(table(x))%>%arrange(desc(Freq))
  
  #################################################
  # Building MLR with handpicked predictors, the ones that capture disturbance:
  #################################################
  
  # the following gauges 2 attributes are kept:
  
  pred_to_keep<-c("term", 
                  "HYDRO_DISTURB_INDX", 
                  "pre1990_STOR", 
                  names(l.G2$Landscape_Pat), 
                  names(l.G2$LC06_Basin),
                  names(l.G2$LC_Crops),
                  'HGA', 'HGB', 'HGC', 'HGD',
                  names(l.G2$Topo)[!grepl("ASPECT|MEAN", names(l.G2$Topo))]
                  )
  
  # I orginally also included in pred_to_keep:
  
  # names(select(l.G2$Climate, ends_with('BASIN'))), names(l.G2$Topo), "BAS_COMPACTNESS", "DRAIN_SQKM", "BFI_AVE","TOPWET"
  
  # removing duplicates of STADID:
  
  pred_to_keep <- pred_to_keep[!(pred_to_keep %in% "STAID")]
  
  # removing predictors not in the correlaiton df columns (some didnt make it):
  
  pred_to_keep <- pred_to_keep[(pred_to_keep %in% names(l.cor.MLR.full[[i]]))]
  
  # filter the correlation df to these predictors

  data.temp<-l.cor.MLR.full[[i]]%>%select(pred_to_keep)

  # the two stepAIC parts from above had really good adj R2 so starting with those to see
  # how model building with these predictor set does:
  
  ### v2.1: check which predictors have the highest coefficent of variability among the 40 observations:
  
  v.mean<-sapply(data.temp[,-1], mean, na.rm=T)
  v.sd<-sapply(data.temp[,-1], sd, na.rm=T)
  v.cv<-v.mean/v.sd
  df.cv<-tibble::rownames_to_column(as.data.frame(v.cv), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.cv$Attribute[1:39]
  # use stepAIC to build a 'full' model using the top n-1:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp%>%select(df.cv$Attribute[1:39])), collapse=" + "))
  m.step.39.cv <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.cv) # adj R2 = 0.7909, pval = 2.076e-05
  # this is pretty good but I still want p<n/10 
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.cv$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-data.temp%>%select(term, var.temp)
  m.leaps.step.cv.39to8to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.cv.39to8to4)
  rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1] 
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1])
  m.leaps.step.cv.39to8to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.cv.39to8to4) # adj r2 = 0.3339, pval = 0.0009889
  # plot(m.leaps.step.cv.39to8to4) # 39, 21, 23 (21 in cooks)
  
  terms.m.step.cv<-names(m.leaps.step.cv.39to8to4$model)[-1]
  
  ### v2.2: check which predictors have the highest variability among the 40 observations:
  
  # need to recreate data.temp:
  
  data.temp<-l.cor.MLR.full[[i]]%>%select(pred_to_keep)
  
  v.var<-sapply(data.temp[,-1], var, na.rm=T)
  df.var<-tibble::rownames_to_column(as.data.frame(v.var), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.var$Attribute[1:39]
  # build the full model using the top n-1:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp%>%select(df.var$Attribute[1:39])), collapse=" + "))
  m.step.39.var <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.var) # adj R2 = 0.9259, pval = 1.01e-08
  # this is amazing but we shouldn't use a model with 22 predictors
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.var$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-data.temp%>%select(term, var.temp)
  m.leaps.step.var.39to22to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.var.39to22to4)
  rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1])
  m.leaps.step.var.39to22to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.var.39to22to4) # adj r2 = 0.3967 , pval = 0.0001962
  # plot(m.leaps.step.var.39to22to4) # 23, 25, 10, 21 (21 in cooks)
  
  terms.m.step.var<-names(m.leaps.step.var.39to22to4$model)[-1]
  
  #################################################
  # standardize predictors
  #################################################
  
  # just running it on the reduced dataset for disturbance predictors:
  
  # need to recreate data.temp and standarize:
  
  data.temp<-l.cor.MLR.full[[i]]%>%select(pred_to_keep)%>%mutate(across(where(is.numeric), ~c(scale(.))))
  
  # now rerun stepAIC+leaps feature selection workflows o thisdataframe:
  
  ### v2.1: check which predictors have the highest coefficent of variability among the 40 observations:
  
  v.mean<-sapply(data.temp[,-1], mean, na.rm=T)
  v.sd<-sapply(data.temp[,-1], sd, na.rm=T)
  v.cv<-v.mean/v.sd
  df.cv<-tibble::rownames_to_column(as.data.frame(v.cv), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.cv$Attribute[1:39]
  # use stepAIC to build a 'full' model using the top n-1:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp%>%select(df.cv$Attribute[1:39])), collapse=" + "))
  m.step.39.cv <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.cv) # adj R2 = 0.7909, pval = 2.076e-05
  # this is pretty good but I still want p<n/10 
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.cv$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-data.temp%>%select(term, var.temp)
  m.leaps.step.cv.39to8to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.cv.39to8to4)
  rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1] 
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1])
  m.leaps.step.cv.39to8to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.cv.39to8to4) # adj r2 = 0.3339, pval = 0.0009889
  # plot(m.leaps.step.cv.39to8to4) # 39, 21, 23 (21 in cooks)
  
  terms.m.step.cv<-names(m.leaps.step.cv.39to8to4$model)[-1]
  
  ### v2.2: check which predictors have the highest variability among the 40 observations:
  
  # need to recreate data.temp:
  
  data.temp<-l.cor.MLR.full[[i]]%>%select(pred_to_keep)
  
  v.var<-sapply(data.temp[,-1], var, na.rm=T)
  df.var<-tibble::rownames_to_column(as.data.frame(v.var), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.var$Attribute[1:39]
  # build the full model using the top n-1:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp%>%select(df.var$Attribute[1:39])), collapse=" + "))
  m.step.39.var <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.var) # adj R2 = 0.9259, pval = 1.01e-08
  # this is amazing but we shouldn't use a model with 22 predictors
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.var$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-data.temp%>%select(term, var.temp)
  m.leaps.step.var.39to22to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.var.39to22to4)
  rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1])
  m.leaps.step.var.39to22to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.var.39to22to4) # adj r2 = 0.3967 , pval = 0.0001962
  # plot(m.leaps.step.var.39to22to4) # 23, 25, 10, 21 (21 in cooks)
  
  terms.m.step.var<-names(m.leaps.step.var.39to22to4$model)[-1]
  
  
  
  
  
  
#########
  
}


# summary table of models: to do this:

# create a list of model names from gloabl enviornment:

rm(m.simple)

l.m<-names(.GlobalEnv) %>% str_subset(pattern = "^m\\.")

# put the models by their model names into a list:

l.m<-do.call("list",mget(l.m))

# filter the list to those with 5 or less terms (1 for intercept):

l.m<-l.m[sapply(l.m, \(i) ifelse(length(coef(i))<6, T, F))]

tab_model(l.m, dv.labels = names(l.m), file="MLR.tab_model.OLS.I.html")

#


############################################################################
                    ####~~~~ MLR standardized ~~~~####
############################################################################

# set up list of dataframes for different selection methods used below:

# create the l..cor list as above but now only filter out the 
# non-significant attributes:

l.cor.MLR<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(Spearman_Correlation)%>%  
           # slice_head(n = 7)%>%
           # bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           # distinct(term, .keep_all = T)%>%
           # filter(!between(Spearman_Correlation, -0.25,.25))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           # mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           # filter(sig_0.05 == 'sig')%>%
           as.data.frame())


# create a list of dataframes from df.OLS.Sens and subseting the attributes using the 
# names in each l.cor[[i]]$term for simple and full models:

l.cor.MLR.full<-lapply(1:4, \(i) df.OLS_Sens%>%
                         select(i+1, l.cor.MLR[[i]]$term)%>%
                         as.data.frame()%>%
                         rename(term = 1))

# standardize predictors (use same name so code infor loop doesnthaveto be changed)

l.cor.MLR.full<-lapply(l.cor.MLR.full, \(i) i%>%mutate(across(where(is.numeric), scale)))


# create another list but keep the site names for use at end of loop in site outliers:

l.cor.MLR.full.w_name<-lapply(1:4, \(i) df.OLS_Sens%>%
                                select(1, i+1, l.cor.MLR[[i]]$term)%>%
                                as.data.frame())

x<-l.cor.MLR.full.w_name[[4]]

# 

# loop:

i<-1

for (i in 1:4){
  
  # regardless of the method to pick the predictors for the MLR model, no more than 4 wil be included as per the run p < n/10
  
  #################################################
  # Spearman correlations::
  #################################################
  
  # build a lm using the top 4 correlates from the spearman correlaiton analysis:
  
  data.temp<-l.cor.MLR.full[[i]][,1:5]
  m.spear<-lm(term~., data = data.temp)
  
  summary(m.spear) # Adj R2 = 0.2074, pval = 0.0156
  
  # plot(m.spear) # 21, 37,23
  
  terms.m.spear<-names(m.spear$model)[-1]
  
  #################################################
  # Random Forest: 
  # changes everytime slightly even if seed is same, 
  # I think this is called unstable
  #################################################
  
  set.seed(1)
  
  # build RF model:
  x<-randomForest(term~., data = l.cor.MLR.full[[i]], na.action=na.omit)
  # arrange by RF variable importance 
  y<-tibble::rownames_to_column(as.data.frame(x$importance), "Attribute")%>%
    arrange(desc(IncNodePurity))
  # build lm of top 4:
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, y$Attribute[1:4])
  m.rf<-lm(term~., data = data.temp)
  summary(m.rf) # adj R2 = 0.2444, pval = 0.0073
  plot(m.rf) # 11, 21, 23, 37 (in cooks)
  
  terms.m.rf<-names(m.rf$model)[-1]
  
  #################################################
  # stepAIC:
  #################################################
  
  #### v1: full model:
  
  # stepAIC workflow:
  fit0 <- lm(term~1,data=l.cor.MLR.full[[i]])
  up.model <- paste("~", paste(colnames(l.cor.MLR.full[[i]][,-1]), collapse=" + "))
  m.step.all <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.all) # doesnt work
  
  #### v2: reduce predictor set prior to building full model:
  
  # scale variables because some have small ranges:
  
  data.temp<-l.cor.MLR.full[[i]][,-1]
  
  ### v2.1: check which predictors have the highest coefficent of variability among the 40 observations:
  
  v.mean<-sapply(l.cor.MLR.full[[i]][,-1], mean, na.rm=T)
  v.sd<-sapply(l.cor.MLR.full[[i]][,-1], sd, na.rm=T)
  v.cv<-v.mean/v.sd
  df.cv<-tibble::rownames_to_column(as.data.frame(v.cv), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.cv$Attribute[1:39]
  # use stepAIC to build a 'full' model using the top n-1:
  fit0 <- lm(term~1,data=l.cor.MLR.full[[i]])
  up.model <- paste("~", paste(colnames(l.cor.MLR.full[[i]]%>%select(df.cv$Attribute[1:39])), collapse=" + "))
  m.step.39.cv <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.cv) # adj R2 = 0.75
  # this is pretty good but I still want p<n/10 
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.cv$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, var.temp)
  m.leaps.step.cv.39to8to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.cv.39to8to4)
  rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1])
  m.leaps.step.cv.39to8to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.cv.39to8to4) # adj r2 = 0.5395, pval = 2.2e-06
  plot(m.leaps.step.cv.39to8to4) # 39, 21, 23 (21 in cooks)
  
  terms.m.step.cv<-names(m.leaps.step.cv.39to8to4$model)[-1]
  
  ### v2.2: check which predictors have the highest variability among the 40 observations:
  
  v.var<-sapply(l.cor.MLR.full[[i]][,-1], var, na.rm=T)
  df.var<-tibble::rownames_to_column(as.data.frame(v.var), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.var$Attribute[1:39]
  # build the full model using the top n-1:
  fit0 <- lm(term~1,data=l.cor.MLR.full[[i]])
  up.model <- paste("~", paste(colnames(l.cor.MLR.full[[i]]%>%select(df.var$Attribute[1:39])), collapse=" + "))
  m.step.39.var <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.var) # adj R2 = 0.9259, pval = 1.01e-08
  # this is amazing but we shouldn't use a model with 22 predictors
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.var$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, var.temp)
  m.leaps.step.var.39to22to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.var.39to22to4)
  rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1])
  m.leaps.step.var.39to22to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.var.39to22to4) # adj r2 = 0.467, pval = 2.5e-05
  plot(m.leaps.step.var.39to22to4) # 23, 25, 10, 21 (21 in cooks)
  
  terms.m.step.var<-names(m.leaps.step.var.39to22to4$model)[-1]
  
  #################################################
  # LEAPS:
  # (straight up with all 306 varables):
  #################################################
  
  # leaps: (doesnt run, too big)
  
  # fit.reg<-leaps::regsubsets(term~.,data=l.cor.MLR.full[[i]],nbest=4, really.big=T)
  
  #################################################
  # LASSO:
  # I dont know the difference between the two lasso methods I have below
  #################################################
  
  #### lasso v1: 
  
  # Define predictor and response variables
  y <- l.cor.MLR.full[[i]]$term
  x <- data.matrix(l.cor.MLR.full[[i]][, -1])
  #fit lasso regression model using k-fold cross-validation
  cv_model <- cv.glmnet(x, y, alpha = 1)
  best_lambda <- cv_model$lambda.min
  #display optimal lambda value
  best_lambda
  #view plot of test MSE's vs. lambda values
  plot(cv_model)
  #view coefficients of best model
  best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  coef(best_model)
  # get non zero coeffs:
  Shat.v1 <- rownames(coef(best_model))[which(coef(best_model) != 0)]
  # since only 3 given, build lm with these!:
  data.temp<-select(l.cor.MLR.full[[i]], term, Shat.v1[-1])
  m.lasso.v1<-lm(term~., data=data.temp)
  summary(m.lasso.v1) # adj r2 = 0.376, pval = 0.00016
  plot(m.lasso.v1) # 23, 15, 21, 28 (28 in cooks)
  
  terms.m.lasso.v1<-names(m.lasso.v1$model)[-1]
  
  #### lasso v2:
  
  #Define predictor and response variables
  ytrain <- l.cor.MLR.full[[i]]$term
  xtrain <- data.matrix(l.cor.MLR.full[[i]][, -1])
  # fit model:
  fit.lasso.glmnet <-glmnet(x=xtrain,y=ytrain,alpha=1) 
  # fit cv (?)
  cv.lasso.glmnet <-cv.glmnet(x=xtrain,y=ytrain,alpha=1) 
  plot(cv.lasso.glmnet)
  # save beta coefficents (?)
  beta.lasso <- coef(fit.lasso.glmnet, s = cv.lasso.glmnet$lambda.min)
  beta.lasso
  # get non zero coeffs:
  Shat.v2 <- rownames(beta.lasso)[which(beta.lasso != 0)]
  # since only 4 given, build lm with these!:
  data.temp<-select(l.cor.MLR.full[[i]], term, Shat.v2[-1])
  m.lasso.v2<-lm(term~., data=data.temp)
  summary(m.lasso.v2) # adj r2 = 0.43, pval = 7.8e-05
  plot(m.lasso.v2) # 21, 39, 23, 28 (28 in cooks)
  
  terms.m.lasso.v2<-names(m.lasso.v2$model)[-1]
  
  #################################################
  # Correlations:
  #################################################
  
  # calculate correlation matrix
  correlationMatrix <- cor(l.cor.MLR.full[[i]][,-1])
  # summarize the correlation matrix
  print(correlationMatrix)
  # ggcorrplot(correlationMatrix)
  # find attributes that are highly corrected (ideally >0.75)
  highlyCorrelated <- caret::findCorrelation(correlationMatrix, cutoff=0.5)
  # print indexes of highly correlated attributes
  print(highlyCorrelated)
  unique(highlyCorrelated)
  # remove highly correlated variables:
  data.temp<-l.cor.MLR.full[[i]][,-c(highlyCorrelated+1)]
  # stepAIC workflow:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp[,-1]), collapse=" + "))
  m.step.cor <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.cor) # adj r2=.4423,pval=0.0006
  # leaps workflow:
  var.temp<-strsplit(as.character(m.step.cor$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-l.cor.MLR.full[[i]]%>%select(term, var.temp)
  m.leaps.step.cor<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.cor)
  rownames(as.data.frame(coef(m.leaps.step.cor, 4)))[-1]
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.cor, 4)))[-1])
  m.leaps.step.cor<-lm(term~., data=data.temp)
  summary(m.leaps.step.cor) # adj r2 = 0.2911, pval = 0.0027
  plot(m.leaps.step.cor) # 23, 1, 21 
  
  terms.m.cor<-names(m.leaps.step.cor$model)[-1]
  
  #################################################
  # Feature Importance
  #################################################
  
  control <- caret::trainControl(method="repeatedcv", number=10, repeats=3)
  # train the model
  model <- caret::train(term~., data=l.cor.MLR.full[[i]], method="leapSeq", preProcess="scale", trControl=control)
  # estimate variable importance
  importance <- varImp(model, scale=FALSE)
  # summarize importance
  print(importance)
  # plot importance
  plot(importance)
  # get top 4 most important:
  x<-tibble::rownames_to_column(importance$importance, "Attribute")%>%arrange(desc(Overall))
  # build lm with top 4:
  data.temp<-select(l.cor.MLR.full[[i]], term, x$Attribute[1:4])
  m.feat_import.4<-lm(term~., data= data.temp)
  summary(m.feat_import.4) # adj r2 = 0.2229, pval=0.01148
  plot(m.feat_import.4) # 23, 21, 11
  
  terms.m.feat_import<-names(m.feat_import.4$model)[-1]
  
  #################################################
  # Feature Selection:
  #################################################
  
  # define the control using a random forest selection function
  control <- caret::rfeControl(functions=rfFuncs, method="cv", number=10)
  # run the RFE algorithm
  results <- caret::rfe(l.cor.MLR.full[[i]][,-1], l.cor.MLR.full[[i]][,1], sizes=c(1:4), rfeControl=control)
  # summarize the results
  print(results)
  # list the chosen features
  predictors(results)
  # plot the results
  plot(results, type=c("g", "o"))
  # build lm with top 4:
  x<-predictors(results)[1:4]
  # build lm with top 4:
  data.temp<-select(l.cor.MLR.full[[i]], term, x)
  m.feat_select.4<-lm(term~., data= data.temp)
  summary(m.feat_select.4) # adj r2 = 0.244, pval=0.0074
  plot(m.feat_select.4) # 23, 21, 27
  
  terms.m.feat_select<-names(m.feat_select.4$model)[-1]
  
  #################################################
  # PCA:
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  #################################################
  
  data.temp<-l.cor.MLR.full[[i]][,-1]
  res.pca <- PCA(data.temp, graph = FALSE)
  print(res.pca)
  eig.val <- get_eigenvalue(res.pca)
  eig.val # sum(eig.val[,1]) # The sum of all the eigenvalues give a total variance of p (305 here).
  fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
  var <- get_pca_var(res.pca)
  var
  fviz_pca_var(res.pca, col.var = "black")
  library("corrplot")
  corrplot(var$cos2, is.corr=FALSE)
  fviz_cos2(res.pca, choice = "var", axes = 1:2)
  fviz_pca_var(res.pca, col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) # Avoid text overlapping)
  fviz_pca_var(res.pca, alpha.var = "cos2")
  corrplot(var$contrib, is.corr=FALSE) 
  # Contributions of variables to PC1
  fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
  # Contributions of variables to PC2
  fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
  # The total contribution to PC1 and PC2 is obtained with the following R code:
  fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)
  # The most important (or, contributing) variables can be highlighted on the correlation plot as follow:
  fviz_pca_var(res.pca, col.var = "contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
  # Create a random continuous variable of length 10
  set.seed(123)
  my.cont.var <- rnorm(305)
  # Color variables by the continuous variable
  fviz_pca_var(res.pca, col.var = my.cont.var,gradient.cols = c("blue", "yellow", "red"),legend.title = "Cont.Var")
  # Create a grouping variable using kmeans
  # Create 3 groups of variables (centers = 3)
  set.seed(123)
  res.km <- kmeans(var$coord, centers = 3, nstart = 25)
  grp <- as.factor(res.km$cluster)
  # Color variables by groups
  fviz_pca_var(res.pca, col.var = grp, palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),legend.title = "Cluster")
  # graph of individuals
  fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE) # Avoid text overlapping (slow if many points))
  # You can also change the point size according the cos2 of the corresponding individuals:
  fviz_pca_ind(res.pca, pointsize = "cos2", pointshape = 21, fill = "#E7B800",repel = TRUE) # Avoid text overlapping (slow if many points)
  # Total contribution on PC1 and PC2
  fviz_contrib(res.pca, choice = "ind", axes = 1:2)
  # biplot:
  fviz_pca_biplot(res.pca, repel = TRUE,col.var = "#2E9FDF", col.ind = "#696969"  )
  
  # Im just going to stop with PCA for now because I dont know how to interpret results:
  
  #################################################
  # Looking at regression outliers:
  #################################################
  
  # The following observations (1-40) came up in the plot(lm) calls throughout thecode above for i=1 (OLS.I):
  # n times (second number)
  
  x<-c(1, 1, 10, 1 , 11, 2, 15, 1, 21, 9, 23, 9, 25, 1, 27, 1, 28, 2, 37, 2, 39, 2)
  x<-as.data.frame(do.call(cbind, split(x, c("Site", "n")))) %>% separate_rows("n", sep = ",")%>%arrange(n)
  
  # more than any others, sites 21 and 23 kept coming up as outliers
  # look at these sites. to do this:
  
  look_at_sites<-l.cor.MLR.full.w_name[[i]][c(21,23),]$Name 
  
  map.1<-df.sf.NWIS%>%filter(Name %in% look_at_sites)
  
  mapview(map.1)
  
  # if I look at the facet plot for CQ type and land use, nothing stands out 
  # about these problem sites. I mean, two are ag and two are undeveloped, but I 
  # dont trust the land use classifications. Also two classes out of 4 isnt that impressive. 
  
  #################################################
  # Looking at frequcny of terms across the models:
  #################################################
  
  x<-grep("terms",names(.GlobalEnv),value=TRUE)
  x<-do.call("list",mget(Pattern1))
  x<-unlist(x)
  x<-as.data.frame(table(x))%>%arrange(desc(Freq))
  
  #################################################
  # Building MLR with handpicked predictors, the ones that capture disturbance:
  #################################################
  
  # the following gauges 2 attributes are kept:
  
  pred_to_keep<-c("term", "HYDRO_DISTURB_INDX", "BAS_COMPACTNESS", "DRAIN_SQKM", names(select(l.G2$Climate, ends_with('BASIN'))),
                  "BFI_AVE", "TOPWET", "pre1990_STOR", names(l.G2$Landscape_Pat), names(l.G2$LC06_Basin),names(l.G2$LC_Crops), names(l.G2$Topo))
  
  # removing duplicates of STADID:
  
  pred_to_keep <- pred_to_keep[!(pred_to_keep %in% "STAID")]
  
  # removing predictors not in the correlaiton df columns (some didnt make it):
  
  pred_to_keep <- pred_to_keep[(pred_to_keep %in% names(l.cor.MLR.full[[i]]))]
  
  # filter the correlation df to these predictors
  
  data.temp<-l.cor.MLR.full[[i]]%>%select(pred_to_keep)
  
  # the two stepAIC parts from above had really good adj R2 so starting with those to see
  # how model building with these predictor set does:
  
  ### v2.1: check which predictors have the highest coefficent of variability among the 40 observations:
  
  v.mean<-sapply(data.temp[,-1], mean, na.rm=T)
  v.sd<-sapply(data.temp[,-1], sd, na.rm=T)
  v.cv<-v.mean/v.sd
  df.cv<-tibble::rownames_to_column(as.data.frame(v.cv), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.cv$Attribute[1:39]
  # use stepAIC to build a 'full' model using the top n-1:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp%>%select(df.cv$Attribute[1:39])), collapse=" + "))
  m.step.39.cv <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.cv) # adj R2 = 0.7909, pval = 2.076e-05
  # this is pretty good but I still want p<n/10 
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.cv$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-data.temp%>%select(term, var.temp)
  m.leaps.step.cv.39to8to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.cv.39to8to4)
  rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.cv.39to8to4, 4)))[-1])
  m.leaps.step.cv.39to8to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.cv.39to8to4) # adj r2 = 0.3339, pval = 0.0009889
  plot(m.leaps.step.cv.39to8to4) # 39, 21, 23 (21 in cooks)
  
  terms.m.step.cv<-names(m.leaps.step.cv.39to8to4$model)[-1]
  
  ### v2.2: check which predictors have the highest variability among the 40 observations:
  
  # need to recreate data.temp:
  
  data.temp<-l.cor.MLR.full[[i]]%>%select(pred_to_keep)
  
  v.var<-sapply(data.temp[,-1], var, na.rm=T)
  df.var<-tibble::rownames_to_column(as.data.frame(v.var), "Attribute")%>%
    rename(var=2)%>%
    arrange(desc(var))
  # the following top n-1 predictors are:
  df.var$Attribute[1:39]
  # build the full model using the top n-1:
  fit0 <- lm(term~1,data=data.temp)
  up.model <- paste("~", paste(colnames(data.temp%>%select(df.var$Attribute[1:39])), collapse=" + "))
  m.step.39.var <- stepAIC(fit0,direction="both",scope=list(lower=fit0,upper=up.model),trace = FALSE)
  summary(m.step.39.var) # adj R2 = 0.9259, pval = 1.01e-08
  # this is amazing but we shouldn't use a model with 22 predictors
  # so use leaps to determine best 4:
  var.temp<-strsplit(as.character(m.step.39.var$terms)[3],split=' + ', fixed=TRUE)[[1]]
  data.temp<-data.temp%>%select(term, var.temp)
  m.leaps.step.var.39to22to4<-regsubsets(term~.,data=data.temp,nbest=1)
  summary(m.leaps.step.var.39to22to4)
  rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1]
  # now build an lm with these 4:
  data.temp<-select(data.temp, term, rownames(as.data.frame(coef(m.leaps.step.var.39to22to4, 4)))[-1])
  m.leaps.step.var.39to22to4<-lm(term~., data=data.temp)
  summary(m.leaps.step.var.39to22to4) # adj r2 = 0.3967 , pval = 0.0001962
  plot(m.leaps.step.var.39to22to4) # 23, 25, 10, 21 (21 in cooks)
  
  terms.m.step.var<-names(m.leaps.step.var.39to22to4$model)[-1]
  
  
  #
  
  
  
  
  
  #########
  
}


# summary table of models: to do this:

# create a list of model names from gloabl enviornment:

rm(m.simple)

l.m<-names(.GlobalEnv) %>% str_subset(pattern = "^m\\.")

# put the models by their model names into a list:

l.m<-do.call("list",mget(l.m))

# filter the list to those with 5 or less terms (1 for intercept):

l.m<-l.m[sapply(l.m, \(i) ifelse(length(coef(i))<6, T, F))]

tab_model(l.m, dv.labels = names(l.m), file="MLR.tab_model.OLS.I.html")

#













































































#### Categorizing land use ####

# USGS criteria:
# Agricultural sites have >50% agricultural land and ≤5% urban land;
# urban sites have >25% urban and ≤25% agricultural land; 
# undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# all other combinations of urban, agricultural, and undeveloped lands are classified as mixed

# combine Ag and Pasture into a single landuse for Ag:

df.datalayers<-mutate(df.datalayers, Ag2 = Ag+Pasture)

# set NA to zero

df.datalayers[is.na(df.datalayers)]<-0

# create the land use class column based on USGS critiera:
# adjusted thresholds:

df.datalayers<-df.datalayers%>%
  mutate(USGS.LU.Adjusted = 'Mixed')%>%
  mutate(USGS.LU.Adjusted = case_when(.default = 'Mixed',
                                      Ag2 > .30 & Developed <= .1 ~ 'Agriculture',
                                      Developed > .1 & Ag2 <= .3 ~ 'Urban',
                                      Developed <= .1 & Ag2 <= .1 ~ 'Undeveloped'))

# merge the OLS and Sens slopes and intercepts with this df:

df.datalayers<-left_join(df.datalayers, df.OLS_Sens[,1:5], by = 'Name')

# plot

p<-df.datalayers%>%
  pivot_longer(cols = 15:18, names_to = 'CQ_parameter', values_to = 'Value')%>%
  ggplot(., aes(x=USGS.LU.Adjusted, y=Value, color =USGS.LU.Adjusted ))+
    geom_boxplot(varwidth = TRUE, alpha=0.2)+
    # scale_x_discrete(labels=my_xlab)+
    facet_wrap('CQ_parameter', scales = 'free')+
    stat_compare_means(method = "anova", label.y = max(df.datalayers$Value))+      # Add global p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = "0.5") +
    ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')

p

#






















#### Categorizing land use (Not run) - experimenting with different data layers ####
# 
# # USGS criteria:
# # Agricultural sites have >50% agricultural land and ≤5% urban land;
# # urban sites have >25% urban and ≤25% agricultural land; 
# # undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# # all other combinations of urban, agricultural, and undeveloped lands are classified as mixed
# 
# # I have four different datasets to pull land use from:
# # NLCD 2001 and 2019
# # CDL 2008 and 2020
# 
# # The plan is to make two facet plots, one for orginal USGS thresholds and one for adjusted USGS thresholds:
# # facets will be the different CQ parameters with x axis being the different land use catgeory andyaxis beung the parameter value
# # for each land use, color will be used to determine which dataset the values came from
# 
# 
# # In a first pass using these thresholds the number of ag and urban sites wasvery low
# # I will play with these numbers to see what happens
# 
# # merge land use from CDL and NLCD to df.sf.NWIS.keep.2 and add a yearcolumn:
# 
# t1<-left_join(df.sf.NWIS.keep.2, df.NWIS.CDL%>%mutate(LU_source = 'CDL_2020'), by = 'Name')
# t2<-left_join(df.sf.NWIS.keep.2, df.NWIS.CDL.2008%>%mutate(LU_source = 'CDL_2008'), by = 'Name')
# t3<-left_join(df.sf.NWIS.keep.2, l.NWIS.NLCD$`2019`%>%mutate(LU_source = 'NLCD_2019'), by = 'Name')
# t4<-left_join(df.sf.NWIS.keep.2, l.NWIS.NLCD$`2001`%>%mutate(LU_source = 'NLCD_2001'), by = 'Name')
# 
# # merge the CDL and NLCD dfs:
# 
# df.LU<-bind_rows(t1,t2,t3,t4)
# 
# # combine Ag and Pasture into a single landuse for Ag:
# 
# df.LU<-mutate(df.LU, Ag = Ag+Pasture)
# 
# # set NA to zero
# 
# df.LU[is.na(df.LU)]<-0
# 
# # create the land use class column based on USGS critiera:
# 
# # orginal thresholds:
# 
# df.LU<-df.LU%>%
#   mutate(USGS.LU = 'Mixed')%>%
#   mutate(USGS.LU = case_when(.default = 'Mixed',
#                              Ag > .50 & Developed <= .05 ~ 'Agriculture',
#                              Developed > .25 & Ag <= .25 ~ 'Urban',
#                              Developed <= .05 & Ag <= .25 ~ 'Undeveloped')
#   )
# 
# # adjusted thresholds:
# 
# df.LU<-df.LU%>%
#   mutate(USGS.LU.Adjusted = 'Mixed')%>%
#   mutate(USGS.LU.Adjusted = case_when(.default = 'Mixed',
#                              Ag > .30 & Developed <= .1 ~ 'Agriculture',
#                              Developed > .1 & Ag <= .3 ~ 'Urban',
#                              Developed <= .1 & Ag <= .1 ~ 'Undeveloped'))
# 
# # merge the OLS and Sens slopes and intercepts with this df:
# 
# df.LU<-left_join(df.LU, df.OLS_Sens[,1:5], by = 'Name')
# 
# # pivot longer for geom_box + facet:
# 
# data<-df.LU%>%
#   pivot_longer(cols = 16:19, names_to = 'CQ_parameter', values_to = 'Value')%>%
#   mutate(USGS.LU=factor(USGS.LU))
# 
# # make ggplot:
# 
# # for orginal thresholds:
# 
# my_xlab <- paste(levels(factor(df.LU$USGS.LU)),"\n(N=",table(factor(df.LU$USGS.LU)),")",sep="")
# 
# ggplot(data, aes(x=LU_source, y=Value, color =USGS.LU ))+
#   geom_boxplot(varwidth = TRUE, alpha=0.2)+
#   # scale_x_discrete(labels=my_xlab)+
#   facet_wrap('CQ_parameter', scales = 'free')+
#   stat_compare_means(method = "anova", label.y = max(data$Value))+      # Add global p-value
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "0.5") +
#   ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')
# 
# # for adjusted thresholds:
# 
# my_xlab <- paste(levels(factor(df.LU$USGS.LU.Adjusted)),"\n(N=",table(factor(df.LU$USGS.LU.Adjusted)),")",sep="")
# 
# ggplot(data, aes(x=LU_source, y=Value, color =USGS.LU.Adjusted ))+
#   geom_boxplot(varwidth = TRUE, alpha=0.2)+
#   # scale_x_discrete(labels=my_xlab)+
#   facet_wrap('CQ_parameter', scales = 'free')+
#   stat_compare_means(method = "anova", label.y = max(data$Value))+      # Add global p-value
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "0.5") +
#   ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')
# 

























#### Grouping CQ curves (stationary, mobilization, dilutionary, complex) ####

# create a list of dataframes for each sites CQ observations:

temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.datalayers$Name)%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))

l.temp<-temp%>%
  split(., .$Name)

# create lm models for each site:

l.lm.CQ_slopes<-lapply(l.temp, \(i) lm(log_C~log_Q, data=i))

# save the model coef ad pvals:

coef<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,1] ))), 'site_no')

# save the pvalues 

pvals<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,4] ))), 'site_no')%>%
  rename(I.pval = 2, S.pval = 3)

# merge the two dfs:

m<-left_join(coef,pvals,by='site_no')

# add column for CQ type:

m<-mutate(m, Type = ifelse(S.pval>0.05, 'Stationary', ifelse(log_Q>0, 'Mobilization', 'Dilution')))

# merge labels with plotting df:

temp<-left_join(temp,m%>%select(site_no, Type),by=c('Name'='site_no'))

# create a df for CQplot of all sites with BP analysis:

df_Seg.2<-filter(df_Seg, site %in% temp$Name)%>%
 left_join(.,m%>%select(site_no, Type),by=c('site'='site_no')) 

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'purple', size = 2)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# looking at this plot I want to add a fourth CQ type for complex, if the slopes of the BP analysis look widely different. 
# I will start wit hcalcuating the angle between pre-post BP slope and see which sites get signled out:

# add a new column with the angle between the two lines:

df_Seg.2<-df_Seg.2%>%
  mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))

# create color pallete for the slope angle:

hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  new_scale_color() +
  geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
  geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
  scale_color_manual(name = "Slope Angle", values = hc)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p

# based on this plot, I would chose the following sites as complex:

complex_sites<-unique(df_Seg.2$site)[c(3,16,21,22,29,37)] 

# old complex site (when n=53): c(3,18,24,25,26,35,42,43)

df_Seg.2<-mutate(df_Seg.2, Type = ifelse(site %in% complex_sites, 'Complex', Type))

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# I want to add land use as color add number of CAFOs as text. to do this:

# merge the plotting df with the land use for NLCD 2001 and CAFO count and adjusted thresholds:

df_Seg.2<-left_join(df_Seg.2, df.datalayers%>%select(Name, USGS.LU.Adjusted, CAFO_count), by = c('site'='Name'))

# if CAFO count is zero set to NA:

df_Seg.2$CAFO_count[df_Seg.2$CAFO_count==0]<-NA

# plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

# plot of trouble sites from MLR analysis:


p<-filter(df_Seg.2, site %in% look_at_sites)%>%
  ggplot(., aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(site), scales = 'free')+
  # theme(
  #   strip.background = element_blank(),
  #   strip.text.x = element_blank()
  # )+
  geom_rect(data = .%>%distinct(.$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

# save(p1, file = 'Processed_Data/p1.Rdata') 
  
# add CQ type to the mapping df with polygon trasperncy based on number of samples:

map.CQ_Type<-df.sf.NWIS.keep%>%
  filter(Name %in% df.datalayers$Name)%>%
  left_join(.,distinct(df_Seg.2, site, .keep_all = T)%>%select(.,c(site, Type, n_sample_rank)),  by = c('Name'='site'))%>%
  select(Name, Type, n_sample_rank)%>%
  arrange(n_sample_rank)%>%
  mutate(NEW = (1/row_number())*40)
  # left_join(.,df.NWIS.TP%>%group_by(site_no)%>%summarise(n=n()), by = c('site'='site_no'))%>%
  # mutate(NEW = case_when(n < 50 ~ .05,
  #                        n >=50 & n < 100 ~ .25,
  #                        n >=100 & n < 500 ~ .5,
  #                        n >=500 ~ .9))
# map:

mapview(map.CQ_Type, zcol = 'Type', alpha.regions = 'NEW')

# map of trouble sites in MLR section:

x<-map.CQ_Type%>%filter(Name %in% look_at_sites)

mapview(x, zcol = 'Type')

######################################################
 #### ~~~~ XY plot of Slopes and Intercepts ~~~~ ####
####################################################

# plot of slopes and intercepts for each CQ type:

m<-m%>%mutate(Type = ifelse(site_no %in% complex_sites, 'Complex', Type))%>%
  left_join(., df.datalayers%>%select(Name, USGS.LU.Adjusted), by = c('site_no'='Name'))

ggplot(m, aes(x=`(Intercept)`, y=log_Q, color = Type))+
  geom_point(size = 2)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  facet_wrap('USGS.LU.Adjusted', scales = 'fixed')
  
# I want to do this with more sites























#### Land Use changes across watersheds (not run) ####

# # Look at changes in NLCD and CDL calculated Ag %land between 2008-2020 for each 53 watershed
# # to do this:
# 
# # pivot longerto getthe land use catgories in a single column:
# 
# df.LU.2<-pivot_longer(df.LU, cols = 6:12, names_to = 'LandUse', values_to = 'values')
# 
# # pivot wider to get the land use source (data layers) in their own columns:
# 
# df.LU.2<-pivot_wider(df.LU.2, names_from = 'LU_source', values_from = 'values')
# 
# # pivot longer to put CDL and NLCD in own columns:
# 
# df.LU.CDL<-pivot_longer(df.LU.2, cols = starts_with('CDL'), names_to = 'CDL', values_to = 'CDL_values')
# df.LU.NLCD<-pivot_longer(df.LU.2, cols = starts_with('NLCD'), names_to = 'NLCD', values_to = 'NLCD_values')
# 
# # group by land use and calculate the difference between the years:
# 
# df.LU.2<-df.LU.2%>%mutate(CDL_diff = abs(round(CDL_2020-CDL_2008,2)),NLCD_diff = abs(round(NLCD_2019-NLCD_2001,2)))%>%
#   arrange(desc(CDL_diff))
# 
























#### Seasonal ANanlysis ####

# is the CQ relaitonship different between seasons at these sites?

# I need to filter down to sites with over 100 samples for this, 20 samples isnt going tocut it:

df_Seg.3<-filter(df_Seg.2, n>100)

# add season column:

df_Seg.3$Season<-getSeason(df_Seg.3$Date)

# plot:

p1<-ggplot(df_Seg.3, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Season), size = 4)+
  # scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.title=element_text(size=14)
  )+
  geom_rect(data = df_Seg.3%>%distinct(df_Seg.3$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = Type), alpha = .15)+
  scale_fill_manual(name = "CQ Type", values = c("red", "blue","purple", "green"))

p1

# save(p1, file = 'Processed_Data/p1.Rdata') 

# 


























############################################
 #### ~~~~ Q Exceedence Proability ~~~~ ####
############################################




























####~~~~ Calcualte EP ~~~~####

# subset the daily flows to the 40 sites:

temp<-filter(df.NWIS.Q, site_no %in% df.datalayers$Name)

# split into lists by site:

temp<-split(temp, f = temp$site_no) 

# sort each dataframe in the list by the flow, add a column for m and n, convert back to a dataframe using bind_rows, and select only the needed columns for the left join (next step):

temp<-bind_rows(lapply(temp, \(i) i[order(i$X_00060_00003,decreasing = T),]%>%mutate(m = 1:n(), n = n())))%>%select(site_no, Date, m, n)

# append the value of M and n for each C-Q observation in the C-Q dataframe:
# need to rename the n column in the left dataframe as well:

df_Seg.2<-left_join(df_Seg.2%>%rename(n_samples = n), temp, by = c("site" = "site_no", "Date" = "Date"))

# calcualte the exceednce proability for Q for each C observation:

df_Seg.2$EP<-round(df_Seg.2$m/(df_Seg.2$n+1), 4)

# plot:

p1<-ggplot(df_Seg.2, aes(x = EP, y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "Log-Log\nCQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  # geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = .5, y =-2, label = CAFO_count), inherit.aes = FALSE, size=15)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  scale_x_reverse()+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p1

save(p1, file = 'Processed_Data/p1.Rdata')

load('Processed_Data/p1.Rdata')

p

#





















####~~~~ Estimate Break point analysis using EP-Q ~~~~####

#  rerun the BP analysis using EP-Q:

l_Seg.EP<-list()
davies.test.matrix.EP<-NULL
EP<-df.NWIS.TP_CQ%>%
  filter(site_no %in% df.datalayers$Name)%>%
  left_join(., df_Seg.2%>%select(site, Date, EP), by = c('site_no'='site', 'sample_dt'='Date'))%>%
  drop_na(result_va, Q_yield)%>%
  rename(site=site_no, Date = sample_dt)
temp.loop<-sort(unique(EP$site))

i<-1

for (i in seq_along(temp.loop)){
  
  tryCatch({
    
    # print the site name for loop debugging:
    
    print(i)
    print(temp.loop[i])
    
    # create a dataframe that will work with segmented. To do this: 
    # filter for the site the loop is in
    # add log transformed C and Q columns, as well as duplicated columns for renamed C and Q
    # filter for real log C and Q values so breakpoint analysis works smoothly:
    
    df<-EP%>%
      filter(site == temp.loop[i])%>%
      mutate(C = result_va)%>%
      mutate(log_C = log(C))%>%
      filter(is.finite(log_C))%>%
      filter(is.finite(EP))
    
    # build a single slope lm for log C and Q. Tis model is also used inthe breakpoint analysis inthenext step:
    
    m<-lm(log_C~EP, df)
    
    # perform breakpoint regression:
    
    m_seg<-segmented(obj = m, npsi = 1)
    
    # perform davies test for constant linear predictor:
    # the results are saved as a string with the the site name and true/false:
    
    x<-paste(temp.loop[i], '-', davies.test(m)$p.val<0.05)
    
    # add the results of davies test to the matrix made prior to this for loop:
    
    davies.test.matrix.EP<-c(davies.test.matrix.EP,x)
    
    # save the breakpoints
    
    bp<-m_seg$psi[1]
    
    # save the slopes: To do this:
    # a conditional statement is needed since sometimes the segmented function wont fit a two slope model at all and will return a object that doesnt work with the slope function used here:
    
    if(length(class(m_seg))==2){
      s<-as.data.frame(slope(m_seg))
    } else{
      s<-NA
    }
    
    # get the intercepts (again conditional statement is needed):
    
    if(length(class(m_seg))==2){
      inter<-as.data.frame(intercept(m_seg))
    } else{
      inter<-NA
    }
    
    
    # get the model fitted data, put in a dataframe, and reformat this dataframe to export out of the loop
    
    if(length(class(m_seg))==2){
      fit <- data.frame(site = df$site, Date = df$Date, Seg_C_EP = fitted(m_seg), I1_EP = inter$Est.[1], I2_EP = inter$Est.[2], Slope1_EP = s[1,1], Slope2_EP = s[2,1], BP_EP = bp)
      result_df<-left_join(df, fit, by = c('site', 'Date'))
    } else{
      fit <- data.frame(site = df$site, Date = df$Date, Seg_C_EP = fitted(m_seg), I1_EP = inter$Est.[1], I2_EP = NA, Slope1_EP = NA, Slope2_EP = NA, BP_EP = NA)
      result_df<-left_join(df, fit, by = c('site', 'Date'))
    }
    
    l_Seg.EP[[i]]<-result_df
    
  }, 
  
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

davies.test.matrix.EP
df.davies.EP<-as.data.frame(davies.test.matrix.EP)%>%
  separate_wider_delim(1, "-", names = c("site", "BP_EP_yes"))%>%
  mutate(across(c(1,2), trimws))
df_Seg.EP<-bind_rows(l_Seg.EP)%>%
  mutate(Seg_EP=EP)
df_Seg.EP<-df_Seg.EP%>%
  left_join(., df.davies.EP, by = 'site')%>%
  mutate(across(c(Seg_C_EP, Seg_EP), ~replace(., BP_EP_yes == 'FALSE', NA)))

# merge the land use, CQ type, data to this df:

df_Seg.EP<-left_join(df_Seg.EP, m%>%select(site_no, Type, USGS.LU.Adjusted), by = c('site'='site_no'))

# plot:

p<-ggplot(df_Seg.EP, aes(x = EP, y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "Log-log\nCQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  new_scale_color()+
  geom_line(aes(x = Seg_EP, y = Seg_C_EP), color = 'yellow', size = 1.5)+
  scale_color_manual(name = "Log(C)~EP-Q\n Breakpoint Analysis", values = "yellow")+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free_y')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  scale_x_reverse()+
  geom_rect(data = df_Seg.EP%>%distinct(df_Seg.EP$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

save(p, file = 'Processed_Data/p1.Rdata')

# I like the idea of using the EP plots to help ID which
# sites are complex. For example, site 22 was deemed complex
# with log-log plots, but in log-EP space it looks similar to
# the other mobilization sites. 























#----------------------------------------#
  ####~~~  Reclassifying CQ type ~~~#### 
   ###~~~  using the EP-Q plots  ~~~### 
#----------------------------------------#

# rerun the df_Seg.2 workflow bt=ut changing the complex_sites vector:
temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.datalayers$Name)%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))
l.temp<-temp%>%
  split(., .$Name)
l.lm.CQ_slopes<-lapply(l.temp, \(i) lm(log_C~log_Q, data=i))
coef<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,1] ))), 'site_no')
pvals<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,4] ))), 'site_no')%>%
  rename(I.pval = 2, S.pval = 3)
m<-left_join(coef,pvals,by='site_no')
m<-mutate(m, Type = ifelse(S.pval>0.05, 'Stationary', ifelse(log_Q>0, 'Mobilization', 'Dilution')))
temp<-left_join(temp,m%>%select(site_no, Type),by=c('Name'='site_no'))
df_Seg.2<-filter(df_Seg, site %in% temp$Name)%>%
  left_join(.,m%>%select(site_no, Type),by=c('site'='site_no')) 
df_Seg.2<-df_Seg.2%>%
  mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))
hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)
complex_sites<-unique(df_Seg.2$site)[c(3,16,21,37)] 
df_Seg.2<-mutate(df_Seg.2, Type = ifelse(site %in% complex_sites, 'Complex', Type))
df_Seg.2<-left_join(df_Seg.2, df.datalayers%>%select(Name, USGS.LU.Adjusted, CAFO_count), by = c('site'='Name'))
df_Seg.2$CAFO_count[df_Seg.2$CAFO_count==0]<-NA

# plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

#

























































# finally save image to workspace:

# save.image(file = 'Processed_Data/NWIS.Rdata')

load('Processed_Data/NWIS.Rdata')




