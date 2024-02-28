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

# read in site metadata and the df.TP_CQ df to get the sites to filter on:

load('Processed_Data/TP.Rdata')

df.NWIS.TP_site_metadata<-read.csv("Raw_Data/df.NWIS.TP_site_metadata.csv", colClasses = c(site_no = "character"))%>%
  filter(site_no %in% df.TP_CQ$site_no)

#### Delinating Watersheds for 42 sites ####

# delinate the sites:

# l.SS_WS.NWIS<-lapply(seq_along(df.NWIS.TP_site_metadata$site_no), \(i) Delineate(df.NWIS.TP_site_metadata$dec_long_va[i], df.NWIS.TP_site_metadata$dec_lat_va[i]))

# names(l.SS_WS.NWIS)<-df.NWIS.TP_site_metadata$site_no

# save(l.SS_WS.NWIS, file = 'C:/PhD/CQ/Downloaded_Data/NWIS.TP.SS_WS.Rdata')

# df.sf.NWIS<-fun.l.SS_WS.to.sfdf(l.SS_WS.NWIS)

# save(df.sf.NWIS, file = 'C:/PhD/CQ/Processed_Data/NWIS.TP.SS_WS_sf.Rdata')

load('C:/PhD/CQ/Processed_Data/NWIS.TP.SS_WS_sf.Rdata')

# which sites didnt work:

# first the ones that didnt make it out of going from SS_WS to df.sf:
  
setdiff(df.NWIS.TP_site_metadata$site_no, df.sf.NWIS$Name)

# then I used the streamstats webpage to get the watersheds for these
# reading them in:

WS_01357500<-st_read('Raw_Data/01357500_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '01357500')
WS_01362200<-st_read('Raw_Data/01362200_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '01362200')
WS_04213500<-st_read('Raw_Data/04213500_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '04213500')

# now combine these and add to df.sf.NWIS:

temp<-bind_rows(WS_01357500, WS_01362200, WS_04213500)

mapview(temp)

df.sf.NWIS<-bind_rows(df.sf.NWIS, temp)

# calculate DA using vect:

vect.NWIS<-vect(df.sf.NWIS)

vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")

# add DA in mi2 to df.sf:

df.sf.NWIS$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102

# add NWIS tabulated area to df sf: first download the metadata for the sites using readNWISsite, then merge drain_area_va column to sfdf:

temp<-readNWISsite(df.sf.NWIS$Name)

df.sf.NWIS<-left_join(df.sf.NWIS, temp[,c(2,30)], by = c('Name'='site_no'))

# calcuate the percent error of the delination:

df.sf.NWIS$Delination_Error<-(df.sf.NWIS$area_sqmi-df.sf.NWIS$drain_area_va)/df.sf.NWIS$drain_area_va

# I want to look at the delimaitons that are on the cusp of not working:

temp<-df.sf.NWIS%>%
  filter(abs(Delination_Error)>.02)

mapview(temp)

temp$Name

# then there are these, the oneswith high delineation error.
# going back to streamstats to get these:

WS_01362500<-st_read('Raw_Data/01362500_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '01362500')
WS_04231600<-st_read('Raw_Data/04231600_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '04231600')
WS_04249000<-st_read('Raw_Data/04249000_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '04249000')

# now combine these and add to df.sf.NWIS. to do this:

# reload in df.sf.NWIS:

load('C:/PhD/CQ/Processed_Data/NWIS.TP.SS_WS_sf.Rdata')

# remove these three:

df.sf.NWIS<-filter(df.sf.NWIS, !Name %in% temp$Name)

# add these three plus the three from before into temp:

temp<-bind_rows(WS_01362500, WS_04231600, WS_04249000, WS_01357500, WS_01362200, WS_04213500)

# now run through the workflow above again:

df.sf.NWIS<-bind_rows(df.sf.NWIS, temp)
vect.NWIS<-vect(df.sf.NWIS)
vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")
df.sf.NWIS$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102
temp<-readNWISsite(df.sf.NWIS$Name)
df.sf.NWIS<-left_join(df.sf.NWIS, temp[,c(2,30)], by = c('Name'='site_no'))
df.sf.NWIS$Delination_Error<-(df.sf.NWIS$area_sqmi-df.sf.NWIS$drain_area_va)/df.sf.NWIS$drain_area_va

# now check the delineation errors:

temp<-df.sf.NWIS%>%
  filter(abs(Delination_Error)>.02)

# good, there are none!

# export watershed shapefile sf df:

# save(df.sf.NWIS, file = 'Processed_Data/NWIS_Watershed_Shapefiles.Rdata')





























#### Delinating Watersheds for second pass ####

# the second pass sites are:

load('Processed_Data/TP.68.Rdata')

# note that 68 is meaning less..i changed the filter because I was wrong the first time

# this is the metadata for these sites:

df.NWIS.TP_site_metadata<-read.csv("Raw_Data/df.NWIS.TP_site_metadata.csv", colClasses = c(site_no = "character"))%>%
  filter(site_no %in% df.TP_CQ.68$site_no)

# there are 42 sites from the previous section that have good delineations
# thus there are 68-42=26 sites that need delineations:
# using the same workflow as above to start getting the watershed shapefiles for these 26 sites:
# note* there are 5 sites in this list that have NA for draiange area so removing those:

temp<-df.NWIS.TP_site_metadata%>%filter(!site_no %in% df.sf.NWIS$Name)%>%drop_na(drain_area_va)

# run streamstats downloads:

# l.SS_WS.NWIS<-lapply(seq_along(temp$site_no), \(i) Delineate(temp$dec_long_va[i], temp$dec_lat_va[i]))

# names(l.SS_WS.NWIS)<-temp$site_no

# save butadd '21' to file name so I know which ones are from this time around (i.e. second pass of this workflow):

# save(l.SS_WS.NWIS, file = 'Downloaded_Data/NWIS.TP.SS_WS.21.Rdata')
# load('Downloaded_Data/NWIS.TP.SS_WS.21.Rdata')

# convert to sf.df:

# df.sf.NWIS.21<-fun.l.SS_WS.to.sfdf(l.SS_WS.NWIS)

# save:

# save(df.sf.NWIS.21, file = 'Processed_Data/NWIS.TP.SS_WS_sf.21.Rdata')
load('Processed_Data/NWIS.TP.SS_WS_sf.21.Rdata')

# which sites didnt work:

# first the ones that didnt make it out of going from SS_WS to df.sf:

setdiff(temp$site_no, df.sf.NWIS.21$Name) # > "04226000" "04264331"

# there are two sites that did not work in the transition to sf.df:
# "04226000" is good to go, but "04264331" is the st. lawrence river, so excluding that one:

# then I used the streamstats webpage to get the watersheds for these
# reading them in:

WS_04226000<-st_read('Raw_Data/04226000_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '04226000')

# now combine these and add to df.sf.NWIS:

temp<-bind_rows(WS_04226000)

mapview(temp)

df.sf.NWIS.21<-bind_rows(df.sf.NWIS.21, temp)

# calculate DA using vect:

vect.NWIS<-vect(df.sf.NWIS.21)

vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")

# add DA in mi2 to df.sf:

df.sf.NWIS.21$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102

# add NWIS tabulated area to df sf: first download the metadata for the sites using readNWISsite, then merge drain_area_va column to sfdf:

temp<-readNWISsite(df.sf.NWIS.21$Name)

df.sf.NWIS.21<-left_join(df.sf.NWIS.21, temp[,c(2,30)], by = c('Name'='site_no'))

# calcuate the percent error of the delination:

df.sf.NWIS.21$Delination_Error<-(df.sf.NWIS.21$area_sqmi-df.sf.NWIS.21$drain_area_va)/df.sf.NWIS.21$drain_area_va

# I want to look at the delimaitons that are on the cusp of not working:

temp<-df.sf.NWIS.21%>%
  filter(abs(Delination_Error)>.02) 

mapview(temp)

# 1st one: 0423204920:
# the SS API gave 7.4 and it is listed as 9. I looked at the USGS page for the 
# gauge and the shape looks the same. I am going to leave it.

# 2nd one: 0424015305

# it is only 5% off, the shape looks right as per comparing to USGS gauge page. So we are all good!




























#### Merging watershed shapefile sf.df's ####

# now I can merge the 42 and 20 shapefiles:

df.sf.NWIS.62<-bind_rows(df.sf.NWIS, df.sf.NWIS.21)

# note* this results in 62 entires
# Inthe table in the manuscript, I though i would have 68,
# but 5 had NA for drainage area (and their gauge names looked weird)
# and one was the st lawrence river
# this leaves 62 = 42 from first pass (sites in G2) and 20 from second pass (sites not in G2)

# now saving this df:

# save(df.sf.NWIS.62, file = 'Processed_Data/NWIS_Watershed_Shapefiles.62.Rdata')

#





















#### finding WWTP in watersheds ####

# read in shapefile of NYS WWTP:

wwtp<-st_read('Raw_Data/NYS_Wastewater_Facility_shapefile/Wastewater_Facility.shp')

# create map of watersheds and WWTP:

mapview(df.sf.NWIS)+mapview(wwtp)

# find intersection of WWTP and watersheds:

x<-st_intersection(wwtp, df.sf.NWIS)

mapview(df.sf.NWIS)+mapview(x)

# looks good
# write out to csv:

write.csv(x[,c(1:6)], row.names = FALSE, file = 'Processed_Data/WWTP_in_NWIS_WS.csv')


























#### finding WWTP in watersheds (for the additional 20 site not in G2) ####

# create map of watersheds and WWTP:

mapview(df.sf.NWIS.21)+mapview(wwtp)

# find intersection of WWTP and watersheds:

x<-st_intersection(wwtp, df.sf.NWIS.21)

mapview(df.sf.NWIS.21)+mapview(x)

# looks good
# write out to csv:

write.csv(x[,c(1:6)], row.names = FALSE, file = 'Processed_Data/WWTP_in_NWIS_WS.21.csv')

