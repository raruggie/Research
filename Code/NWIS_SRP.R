# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

# just load the packages section in NWIS.R

####################### Functions #######################

source("Code/Ryan_functions.R")

####################### Goal of code #######################

# 1) Process the NWIS database for the CQ analysis
# 2) Run through CQ metrics/watershed attribute correlation

####################### Workflow #######################

# SRP pcode:

pcode<-'00660' 

ncode<-'SRP'
















#### NWIS Query ####

# this code adapted from NYS_site_ranking.R

# Download the metadata for sites with daily flow data. To do this:
# alread did this in Code/NWIS.R:

df.NWIS.Q_sites <- read.csv("Raw_Data/df.NWIS.Q_sites.csv", colClasses = c(site_no = "character"))

# clean up this dataframe

# there are duplicates of site numbers with different numbers of daily flow observations
# so group by the site no numbers and just keep the observation with the highest number of flow observations
# I realize now that this is proably not needed since I am grouping by the site number so the resulting number used in the
# raw flow data download will have what it has, but I am just goin to keep this (aka i could just use distinct?)
# also rename some columns and convert them to numeric:

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

# From the sites with daily flow data, determine which sites also have over 99 discrete samples . To do this:

# use function I created (sourced) to get a dataframe of sites for just SRP with #samples = 20 threshold:

# df.NWIS.SRP_sites<-fun.df.Pair_consit_flow(pcode, df.NWIS.Q_sites, n_samples = 20, state = 'NY')

# write.csv(df.NWIS.SRP_sites, "Raw_Data/df.NWIS.SRP_sites.csv", row.names=FALSE)

df.NWIS.SRP_sites<-read.csv("Raw_Data/df.NWIS.SRP_sites.csv", colClasses = c(site_no = "character"))

# download the raw daily flow data for these sites. To do this:
# I already downloaded the daily flow data for the TP and TN sites
# which ever sites do notoverlap with SRP will be downloaded here. To do this:

# read in the flow data for the TP and TN sites:

# df.NWIS.Q.combined<-read.csv("Raw_Data/df.NWIS.Q.combined.csv", colClasses = c(site_no = "character"))%>%mutate(Date = as.Date(Date))

# determine which of the SRP sites dont overlap:

# df.overlap<-df.NWIS.SRP_sites%>%filter(!site_no %in% df.NWIS.Q.combined$site_no)

# download the raw flow data for the sites still needed:

# df.NWIS.Q.SRP_needed<-readNWISdv(siteNumbers = df.overlap$site_no, parameterCd = '00060', startDate = "", endDate = "", statCd = "00003")

# unique(df.NWIS.Q.SRP_needed$site_no) # only get 1 of 19 that need... sucks but idk

# merge with the SRP sites flow data:

# df.NWIS.Q.combined<-bind_rows(df.NWIS.Q.combined,df.NWIS.Q.SRP_needed)

# write out this df to file: (be carefull!! writing over!!)

# write.csv(df.NWIS.Q.combined, "Raw_Data/df.NWIS.Q.combined.csv", row.names=FALSE)

# read it back in as 'df.NWIS.Q.for_SRP':

df.NWIS.Q.for_SRP<-read.csv("Raw_Data/df.NWIS.Q.combined.csv", colClasses = c(site_no = "character"))%>%mutate(Date = as.Date(Date))

# filter to just the SRP sites:

df.NWIS.Q.for_SRP<-df.NWIS.Q.for_SRP%>%filter(site_no %in% df.NWIS.SRP_sites$site_no)

# download the raw discrete SRP sample data:

# df.NWIS.SRP<-readNWISqw(siteNumbers = df.NWIS.SRP_sites$site_no, parameterCd = pcode)

# write.csv(df.NWIS.SRP, "Raw_Data/df.NWIS.SRP.csv", row.names=FALSE)

df.NWIS.SRP<-read.csv("Raw_Data/df.NWIS.SRP.csv", colClasses = c(site_no = "character"))%>%mutate(sample_dt = as.Date(sample_dt))

#























#### Building CQ df ####

# join the SRP with the flow df:

df.NWIS.SRP_CQ<-left_join(df.NWIS.SRP, df.NWIS.Q.for_SRP, by=c("site_no"="site_no", "sample_dt"="Date"))

# take average of multiple samples on the same day at the same site:

df.NWIS.SRP_CQ<-df.NWIS.SRP_CQ%>%
  group_by(site_no, sample_dt)%>%
  summarise_at(vars(result_va, X_00060_00003), funs(mean(., na.rm=TRUE)))%>%
  ungroup()

# remove sites where there are not >= 20 CQ pairs:

df.NWIS.SRP_CQ<-df.NWIS.SRP_CQ%>%
  drop_na(result_va, X_00060_00003)%>%
  group_by(site_no)%>%
  mutate(n=n())%>%
  ungroup()%>%
  filter(n>=20)%>%
  select(-n) # remove n column because going to get added back in next steps

# add column for number of paired observations by rank: to do this:

# create column for n_sample_rank:

temp<-df.NWIS.SRP%>%group_by(site_no)%>%
  summarise(n=n())%>%
  arrange(desc(n))%>%
  mutate(n_sample_rank=rank(-n, ties.method='first'))

# merge this df with df.NWIS.TP_CQ and arrange by the new column:

df.NWIS.SRP_CQ<-left_join(df.NWIS.SRP_CQ,temp, by='site_no')%>%
  arrange(n_sample_rank)

# next step is to determine which sites have draiange area listed
# if no draiange area is listed, there is no way to determine if delineation
# is correct. to do this:

# download the drainage areas from site metadata using readNWISsite

# df.NWIS.SRP_site_metadata<-readNWISsite(siteNumbers = unique(df.NWIS.SRP_CQ$site_no))

# write.csv(df.NWIS.SRP_site_metadata, "Raw_Data/df.NWIS.SRP_site_metadata.csv", row.names=FALSE)

df.NWIS.SRP_site_metadata<-read.csv("Raw_Data/df.NWIS.SRP_site_metadata.csv", colClasses = c(site_no = "character"))

# then select just the site number and DA column:

df.DA<-df.NWIS.SRP_site_metadata%>%
  select(site_no, drain_area_va)

# finally merge with df.NWIS.TP_CQ and create a new Q column with area normalized flows (not worrying about units right now): 
# Note: Some sites returned NA on draiange areas in readNWISsite:

df.NWIS.SRP_CQ<-left_join(df.NWIS.SRP_CQ, df.DA, by = 'site_no')%>%
  mutate(Q_yield = X_00060_00003/drain_area_va)%>%
  drop_na(drain_area_va)

# length(unique(df.NWIS.SRP_CQ$site_no))  97 sites


























#### Apply filters to sites ####

# number of sites with >=20 paired CQ observations:

length(unique(df.NWIS.SRP_CQ$site_no))

# 1) sites with samples after 2001:

# note that you need to rerun the 
# number of sample filter again because some sites 
# have data pre and post 2001

temp1<-df.NWIS.SRP_CQ%>%filter(sample_dt >= as.Date('2001-01-01'))%>%
  group_by(site_no)%>%
  mutate(n=n())%>%
  filter(n>=20)

length(unique(temp1$site_no))

# 59 sites

# 2) sites not on long island

# use latitude to filter:

temp2<-temp1%>%filter(site_no %in% filter(df.NWIS.SRP_site_metadata, dec_lat_va >40.9364)$site_no)

length(unique(temp2$site_no))

# 58 sites

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

# 40 sites

# one last test: redoing the filters but with no removing date restriction:

# number of sites in gauges 2:

x1<-filter(df.NWIS.SRP_CQ, site_no %in%df.G2$STAID)

length(unique(x1$site_no))

# 69 sites

# number of sites in gauges 2 and not on LI:

x2<-x1%>%filter(site_no %in% filter(df.NWIS.SRP_site_metadata, dec_lat_va >40.9364)$site_no)

length(unique(x2$site_no))

# 58 sites

























#### ! Finalizing df.TN_CQ to move on in code ####

# first pass 

df.SRP_CQ<-temp3

# second pass

df.SRP_CQ.SP<-temp2

# going to write over df.TN_CQ for this code to work smoothly with the these sites

df.SRP_CQ<-df.SRP_CQ.SP

#



























#### Restrictions on Data ####

# 1) remove the 4-5 largest sites:

x <- df.DA %>% 
  rename(Name = 1,Group = 2) %>% 
  filter(Name %in% df.SRP_CQ$site_no)

fun.map.DA(x$Name, x)

x <- x %>% 
  slice_min(Group, n=-4)

fun.map.DA(x$Name, x)

df.SRP_CQ <- df.SRP_CQ %>% filter(site_no %in% x$Name) 

# # 1) 50 samples:
# 
# df <- df.TN_CQ %>% filter(n>=50)
# 
# length(unique(df$site_no))
# 
# # 11?!?! 
# 
# # 2) over 20% of samples are above the 90% flow percentile:
# 
# # determine the 90% flow for each site:
# 
# df.Q.range <- df.NWIS.Q.for_TN %>% 
#   filter(site_no %in% df$site_no) %>% 
#   group_by(site_no) %>% 
#   dplyr::reframe(CI = quantile(X_00060_00003,c(.90), na.rm = T))
# 
# # merge this with the CQ df:
# 
# df.Q.range <- left_join(df, df.Q.range, by = 'site_no') 
# 
# # filter to samples above the 90% threshold column, add column, and filter:
# 
# df.Q.range <- df.Q.range%>% 
#   filter(X_00060_00003 >= CI) %>% 
#   mutate(n_90 = n()) %>% 
#   mutate(perc_90 = n_90/n) %>% 
#   distinct(site_no, perc_90) %>% 
#   arrange(perc_90) %>% 
#   filter(perc_90 >= .20)
# 
# fun.map.DA(df.Q.range$site_no)
# 
# df.TN_CQ <- df.TN_CQ %>% filter(site_no %in% df.Q.range$site_no)







































#### Fitting Breakpoints to CQ curves #### 

# use funciton in lapply to get list of df of segmented results:

l_Seg<-lapply(seq_along(unique(df.SRP_CQ$site_no)), \(i) fun.CQ.BP(i, df.SRP_CQ))

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

#






















#### Setting up Response Variables ####

#### ~ Set up df for this section: ####

df.RP <- df.SRP_CQ%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003, n_sample_rank)%>%
  mutate(sample_dt = as.Date(sample_dt), log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))

l.RP<-split(df.RP, f=df.RP$Name)

df.Q<-df.NWIS.Q.for_SRP%>% # this gives df of complete hydrographs for each site and includes columns to group by calender day later on
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

l.DA<-df.NWIS.SRP_site_metadata %>% filter(site_no %in% names(l.Q))%>%select(site_no, drain_area_va)%>%split(., .$site_no)%>%lapply(., \(i) i$drain_area_va)

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

df.median.Q<-df.RP%>%group_by(Name)%>%summarise(median_Q=median(X_00060_00003))

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
names(l.lm.pre)==df.median.Q$Name

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
names(l.lm.pre)==df.median.Q$Name

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

#### ~ 5) Flow Weighted Concentrations ####

# two approaches here:

# 1) using the individual CQ observations (Mean annual flow weighted concentration, MAFWC)
# 2) predicting C for eacch day in the record (Flow weighted Average concentration, FWAC)

# 1) MAFWC

# estimate load for each CQ observations (i.e. C*Q), 
# then divide each load observation by the average annual flow for the site
# then take the mean of these weighted concentrations:

# estimate average annual flow for each site:

df.AAF <- df.Q %>% group_by(site_no, year) %>% summarize(AAF = mean(X_00060_00003, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(site_no) %>% summarize(AAF = mean(AAF, na.rm = T))

# calculate loads for observed CQ pairs:

df.loads <- df.RP %>% mutate(load = C*Q)

# merge AAF for each site with df.load:

df.loads <- left_join(df.loads, df.AAF, by =c('Name'= 'site_no'))

# divide each load observation by AAF:

df.loads <- df.loads %>% mutate(MAFWC = load/AAF)

# take mean weighted C for each site:

df.MAFWC <- df.loads %>% group_by(Name) %>% summarize(MAFWC = mean(MAFWC, na.rm = T))

# 2) FWAC

# using l.pred.C from hydrosep
# combine l.pred.C into single df:

df.predC <- bind_rows(l.pred.C, .id = 'Name')

# calcualte load for everyday in record:

df.load <- df.predC %>% mutate(load = X_00060_00003*pred.C)

# group by and summarize to estimate (annual) mean daily load for each site:

df.load <- df.load %>% group_by(Name) %>% summarize(load = mean(load, na.rm = T))

# merge AAF for each site with df.load:

df.load <- left_join(df.load, df.AAF, by = c('Name'= 'site_no'))

# divide loads by average annual flow:

df.FWAC <- df.load %>% mutate(FWAC = load/AAF, .keep = 'unused')

#### ~ Putting response variables together into single df: ####

df.Response <- left_join(df.OLS, df.CV, by = 'Name')%>%
  left_join(.,df.medC, by = 'Name') %>% 
  left_join(., df.MAFWC, by = 'Name') %>% 
  left_join(., df.FWAC, by = 'Name')

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
                           grepl("simple", name) ~ "simple",
                           grepl("MAFWC", name) ~ "MAFWC",
                           grepl("FWAC", name) ~ "FWAC"
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

v.resp.cols <- 2:20

v.pred.cols <- 21:ncol(df.setup)

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
  corrr::focus(2:20)%>%
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

# I also want to look at the relaitionship between mean c and MAFWC:

x <- left_join(df.setup, df.DA, by = c('Name'='site_no')) %>% 
  select(Name,meanC, MAFWC, FWAC, drain_area_va)

p1 <- ggplot(x, aes(x = meanC, y = MAFWC))+
  geom_point(aes(size = drain_area_va))+
  geom_line(aes(x = meanC, y = meanC), color = 'red') +
  geom_smooth(method = 'lm')+
  scale_size_continuous(breaks=c(25, 250, 500, 1000, 1250, 2000, 5000))

p2 <- ggplot(x, aes(x = meanC, y = FWAC))+
  geom_point(aes(size = drain_area_va))+
  geom_line(aes(x = meanC, y = meanC), color = 'red') +
  geom_smooth(method = 'lm')+
  scale_size_continuous(breaks=c(25, 250, 500, 1000, 1250, 2000, 5000))

p3 <- ggplot(x, aes(x = MAFWC, y = FWAC))+
  geom_point(aes(size = drain_area_va))+
  geom_line(aes(x = MAFWC, y = MAFWC), color = 'red') +
  geom_smooth(method = 'lm')+
  scale_size_continuous(breaks=c(25, 250, 500, 1000, 1250, 2000, 5000))

ggarrange(p1,p2,p3, common.legend = T)

#



























#### PCA ####

# ID and then map sites that are outliers:

out.pos <- lapply(df.setup[,v.resp.cols], FindOutliers) %>% unlist %>% unique
outlier.sites <- df.setup$Name[out.pos]

fun.map.DA(outlier.sites)

# set up a df for PCA:

df.PCA <- df.setup
rownames(df.PCA) <- df.PCA[,1] # set the row names as the site names

# determine the top predictors based on variance:

preds.var <- df.PCA[,v.pred.cols] %>% pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  dplyr::summarize(var = var(value)) %>% arrange(desc(var)) %>% ungroup() %>% top_n(16)

# also do variances based on standarized 0-1 variables:

preds.var.01 <- df.PCA[,v.pred.cols] %>%
  lapply(., range01) %>%
  bind_cols() %>%
  pivot_longer(cols = everything()) %>%
  group_by(name) %>%
  dplyr::summarize(var = var(value)) %>% arrange(desc(var)) %>% ungroup() %>% top_n(16)

# lets look at the differences between the predictor sets for regular and standardized (0-1) top variance:

(x <- preds.var$name)
(y <- preds.var.01$name)

setdiff(x,y) # "mean.daily.Q" "Elev_Avg"     "Elev_Median"  "Sum.WWTP.Q"   "HSG_C" 

setdiff(names(df.setup[,v.pred.cols]),y)

setdiff(y,x) # "R_CROPSNLCD06" "Corn" "R_RIP100_PLANT" "RIP.CSA.100" "BFI"  

# ready for PCA

# there are lots of iterations to be had here:

# 1) all sites with all predictors
# 2) all sites with top half of predictors based on real variance
# 3) all sites with top half of predictors based on 0-1 variance
# 4) remove outlier sites with all predictors
# 5) remove outlier sites with top half of predictors based on real variance
# 6) remove outlier sites with top half of predictors based on 0-1 variance

pc1 <- prcomp(df.PCA[,v.pred.cols], center = TRUE, scale. = TRUE)
pc2 <- prcomp(df.PCA[,v.pred.cols] %>% select(preds.var$name), center = TRUE, scale. = TRUE)
pc3 <- prcomp(df.PCA[,v.pred.cols] %>% select(preds.var.01$name), center = TRUE, scale. = TRUE)
pc4 <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name), center = TRUE, scale. = TRUE)
pc5 <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name) %>% select(preds.var$name), center = TRUE, scale. = TRUE)
pc6 <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name) %>% select(preds.var.01$name), center = TRUE, scale. = TRUE)

# make individuals biplots:

l.pca <- list(pc1,pc2,pc3,pc4,pc5,pc6)

l.pca.plot.biplot.obs <- lapply(l.pca, fun.PCA.biplot.individuals)

# look at individuals biplots for each PCA done above,
# as well as the cummulative proportion of explained variance at the second PC
# as well as the PC at which 95% of the variability is explained:

l.pca.plot.biplot.obs[[1]]
as.data.frame(summary(pc1)$importance)$PC2[3]
min(which(summary(pc1)$importance[3,]>.95))

l.pca.plot.biplot.obs[[2]]
as.data.frame(summary(pc2)$importance)$PC2[3]
min(which(summary(pc2)$importance[3,]>.95))

l.pca.plot.biplot.obs[[3]]
as.data.frame(summary(pc3)$importance)$PC2[3]
min(which(summary(pc3)$importance[3,]>.95))

l.pca.plot.biplot.obs[[4]]
as.data.frame(summary(pc4)$importance)$PC2[3]
min(which(summary(pc4)$importance[3,]>.95))

# l.pca.plot.biplot.obs[[5]]
as.data.frame(summary(pc5)$importance)$PC2[3]
min(which(summary(pc5)$importance[3,]>.95))

l.pca.plot.biplot.obs[[6]]
as.data.frame(summary(pc6)$importance)$PC2[3]
min(which(summary(pc6)$importance[3,]>.95))


# i want to make a map based on the three groups I am seeing in pc3 (all sites, restrctions on predictors):

# the first group will be greater than 0 for PC1 and greater than -0.5 for pc2
# the second group will be greater than 0 for PC1 and less than -0.5 for PC2
# the third group will be all sites less than 0 for pc1

l.pca.plot.biplot.obs[[2]]

fviz_pca_var(pc2,
             axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


df.pc3 <- as.data.frame(pc3$x)

g1<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2>-0.5,])
g2<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2<=-0.5,])
g3<- rownames(df.pc3[df.pc3$PC1<0,])

v.g <- c(g1,g2,g3)
v.g.names <- c(rep('g1', length(g1)),rep('g2', length(g2)),rep('g3', length(g3)))

df.g <- data.frame(v.g, v.g.names) %>% rename(Name = 1, Group = 2)

fun.map.DA(df.g$Name, zcol.key = df.g)

# what if I remove cattargus and genese river, because they are so weird/big and rerun pc3?
# this is basically pc6 with only 2 outlier sites now:

outlier.sites <- c("04231600","04213500")

pc3a <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name) %>% select(preds.var.01$name), center = TRUE, scale. = TRUE)

as.data.frame(summary(pc3a)$importance)$PC2[3]
min(which(summary(pc3a)$importance[3,]>.95))

fun.PCA.biplot.individuals(pc3a)

# the groupings look much different

df.pc3 <- as.data.frame(pc3a$x)

g1<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2>-1,])
g2<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2<=-1,])
g3<- rownames(df.pc3[df.pc3$PC1<0,])

v.g <- c(g1,g2,g3)
v.g.names <- c(rep('g1', length(g1)),rep('g2', length(g2)),rep('g3', length(g3)))

df.g <- data.frame(v.g, v.g.names) %>% rename(Name = 1, Group = 2)

fun.map.DA(df.g$Name, zcol.key = df.g)

# I dont see much difference between 1 and 2 (I dont know why just three are in 2)
# I suspect it could be something with near stream CSA,
# so I want to look at that map:

df.rip.csa <- df.setup %>% 
  filter(Name %in% df.g$Name) %>% 
  select(Name, R_PLANTNLCD06) %>% 
  rename(Group = 2)

fun.map.DA(df.rip.csa$Name, df.rip.csa)

# so those three stand out because they have much higher watershed percent agriculture

# the issue now becomes we only have three agricultural sites

# lets look at map of all 62 sites (minus the largest ones) or whatever based on agriculture:

outlier.sites <- c(outlier.sites,  c('04260500', '04231600', '01357500', '04249000'))

df.map <- df.datalayers.62 %>% 
  filter(!Name %in% outlier.sites) %>% 
  select(Name, R_PLANTNLCD06) %>% 
  rename(Group = 2)

fun.map.DA(df.map$Name, df.map)

# it doesnt look there are many more sites to add back in that have similar ag to the three I kept...

# lets look at the biplot of predictor variables:

fviz_pca_var(pc3a,
             axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# if we look at the contributions, forest is strongest in the negative direction for Dim1, and plant is strongest in the positive direction
# so can we say that Dim1 is a good contender variable for high correlation with positie nutrient loads?

# lets look at the actual ranks of importance:

df.importance <- as.data.frame(pc3a$rotation) %>% arrange(PC1, desc = F)

# the PCA for TN isnt too important because I am going to use these sites anyway...
# but I'm still going to use the plots for the manuscript





























#### Set up list to export for hypothesis testing ####

# set up list of dataframes for response variables and predictors

# OLS slope
# OLS slope stormflow samples
# MAFWC
# FWAC
# AANY.hydrosep
# medC

resp.vars <- c("OLS.Slope.1s", "OLS.Slope.2s_hydrosep_post", "MAFWC", 'FWAC', "OLS.AANY.2s.hydrosep.method2", "medC")

# One last check of the site restrictions:
# I want to make sure that at least 25 samples make up the stormflow segmenet:

test <- df.RP %>% 
  filter(Name %in% df.g$Name) %>% 
  filter(Q_type == 'Stormflow') %>% 
  group_by(Name) %>% 
  summarise(n_storm=n()) %>% 
  arrange(n_storm)

test1 <- df.SRP_CQ %>% 
  distinct(site_no, n) %>% 
  rename(Name = 1) %>% 
  left_join(., test, by = 'Name') %>% 
  mutate(frac = n_storm/n) %>% 
  arrange(frac)

ggplot(df.RP, aes(x = log_Q, y = log_C, color = Q_type)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap(~Name, scales = 'free')

# looks good

# set up dataframes for the 4 response variables::

l.resp.allpred <- lapply(resp.vars, \(i) df.PCA %>% filter(Name %in% df.g$Name) %>% select(i, v.pred.cols) %>% rename(term = 1)) %>% purrr::set_names(resp.vars)

# save:

save(l.resp.allpred, file = 'Processed_Data/SRP.l.resp.allpred.Rdata')

#






























#### Calculating Average Annual Consituent Yield ####
# (for use in correlaitons and MLR)

# step 1. calcualte average annual hydrograph for each site. to do this:

# filter the daily flow records to just the sites used:

l.avg_ann_hydro<-df.NWIS.Q.for_SRP%>%filter(site_no %in% df.SRP_CQ$site_no)

# split the daily flow records into list of dfs

l.avg_ann_hydro<-split(l.avg_ann_hydro, f = l.avg_ann_hydro$site_no)

# strip out year from date (just converts it to 2023) for each dataframe, then group by and summarize to get mean annual flow for each day of calendar year:

l.avg_ann_hydro<-lapply(l.avg_ann_hydro, \(i) i%>%mutate(Date = as.Date(Date))%>%mutate(year = year(Date), Date = as.Date(format(Date, "%m-%d"), "%m-%d"))%>%group_by(Date)%>%summarize(mean_Q = mean(X_00060_00003, na.rm= T)))

# step 2. predict C for each of day of the year using the mean flow rate for each site in the list:todothis:

# create a list of linear model objects for each sites CQ relaitonship:

l.C_daily<-df.SRP_CQ%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))%>%
  split(., f =.$Name)%>%
  lapply(., \(i) lm(log_C ~ log_Q, i))

# make sure the order of the sites in the flow and lm lists are the same:

names(l.C_daily)==names(l.avg_ann_hydro)

# use this list to create list of bias correction factors:

l.BCF<-lapply(l.C_daily, \(i) exp((summary(i)$sigma^2)/2))

# predict the concentration for each day for each site: this predicts the concnetration in log space, so transforming into real space here as well using BCF:

l.C_daily<-lapply(names(l.C_daily), \(i) l.BCF[[i]]*exp(as.numeric(predict(l.C_daily[[i]], data.frame(log_Q = log(l.avg_ann_hydro[[i]]$mean_Q) ) ) ) ) )%>%purrr::set_names(names(l.avg_ann_hydro))

# step 3. multiply predicted average daily C by daily Q to get daily load, sum up, then convert to yield. To do this:

# save just the daily Q:

l.Q.temp<-lapply(l.avg_ann_hydro, \(i) i$mean_Q)

# set up list of drainage areas:

l.DA<-readNWISsite(siteNumbers = names(l.C_daily))%>%select(site_no, drain_area_va)%>%split(., .$site_no)%>%lapply(., \(i) i$drain_area_va)

# multiply daily C by daily Q and convert units:
# mg/L * ft^3/sec * 28.3168 L/ft3 * 86400 sec /day * 1g/1000mg * 1kg/1000g * 1/mi^2 * 1mi^2/258.999 ha = kg/ha/day

unit_conversion=28.3168*86400*(1/1000)*(1/1000)*(1/258.999)

l.Yield<-lapply(names(l.C_daily), \(i) l.C_daily[[i]]*l.Q.temp[[i]]*(1/l.DA[[i]])*unit_conversion)%>%purrr::set_names(names(l.avg_ann_hydro))

# sum up for the year and convert to df

df.Yield <- lapply(l.Yield, sum)%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'Yield') # units of kg/ha/year

# compare with TP df.Yield:

# rename df.Yield:

df.Yield_SRP<-df.Yield%>%rename(Yield_SRP = Yield)

# load TP df.Yield

load('Processed_Data/df.Yield_TP.Rdata')

# rename TP df.Yield:

df.Yield_compare<-df.Yield%>%rename(Yield_TP = Yield)

# write over df.Yield with SRP data and rename column (for use later in this file, I want df.Yield to be like it was before comparing with TP)

df.Yield<-df.Yield_SRP%>%rename(Yield = Yield_SRP)

# create df.Yield_Compare:

df.Yield_compare<-left_join(df.Yield_compare, df.Yield_SRP, by = 'Name')%>%mutate(TPminusSRP = Yield_TP-Yield_SRP)

#




















#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Correlations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# reduce the predictor set down to those we believe are defensible:

# I also went in to the variable description excel sheet and isolated these variables as well:

temp1<-read.csv('Processed_Data/G2_pred_to_keep_2.csv')[,1]

# combining these with the vairbales in entire sheets as well:

pred_to_keep<-c("HYDRO_DISTURB_INDX", 
                "pre1990_STOR", 
                names(l.G2$Landscape_Pat), 
                names(l.G2$LC06_Basin),
                names(l.G2$LC_Crops),
                temp1,
                'HGA', 'HGB', 'HGC', 'HGD',
                "ELEV_MEDIAN_M_BASIN",
                "ELEV_STD_M_BASIN",
                "RRMEDIAN"
)

# removing 'STADID':

pred_to_keep <- pred_to_keep[!(pred_to_keep %in% "STAID")]

# removing the duplicates of DEV and PDEN, HIRES

pred_to_keep <- pred_to_keep[!(pred_to_keep %in% c("DEVOPENNLCD06","DEVLOWNLCD06","DEVMEDNLCD06","DEVHINLCD06", "PDEN_DAY_LANDSCAN_2007", "PDEN_NIGHT_LANDSCAN_2007", "HIRES_LENTIC_NUM"))]

# filtering df.G2 to this predictor list:

df.G2.reduced<-df.G2%>%select(c('STAID', pred_to_keep))

# take the average of the RIP and MAIN for each land use:

df.G2.reduced<-df.G2.reduced%>%filter(STAID %in% df.SRP_CQ$site_no)%>%
  rowwise() %>%
  mutate(MAIN_RIP_DEV_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("DEV")), na.rm = TRUE),
         MAIN_RIP_FOREST_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("FOREST")), na.rm = TRUE),
         MAIN_RIP_PLANT_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("PLANT")), na.rm = TRUE),
         .keep = 'unused'
  )

# adding up the other crops:

df.G2.reduced<-df.G2.reduced%>%
  mutate(sum_of_misc_crops = sum(c_across(ends_with(c("COTTON", "RICE", "SORGHUM", "SUNFLOWERS", "PEANUTS", "BARLEY", "DURUM_WHEAT", "OATS", "DRY_BEANS", "POTATOES", "ORANGES", 'OTHER_CROPS', 'IDLE', "WWHT_SOY_DBL_CROP"))), na.rm = TRUE), .keep = 'unused')

names(df.G2.reduced)

# combined the reduced predictors df with the already filtered CQ data:

df.Correlation<-df.SRP_CQ%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))%>%
  left_join(., df.G2.reduced, by = c('Name'='STAID'))

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

df.OLS_Sens<-left_join(df.OLS_Sens, df.G2.reduced, by = c('Name'='STAID'))

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

# I also want to try standardizing the predictors:

df.cor.0to1<-df.OLS_Sens %>% 
  mutate(across(7:last_col(), ~ (.-min(.))/(max(.)-min(.))))%>% # could add term: .names = "{.col}_standarized")
  correlate(method = 'spearman') %>%
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S, Yield))%>%
  pivot_longer(cols= c(2:6), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# looking at the regualr and 0 to 1 standarized spearman correlaitons, they are the same

df.cor$Spearman_Correlation==df.cor.0to1$Spearman_Correlation

# so not going to proceed with 0 to 1

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

OLS<-l.cor[[5]]%>%arrange(desc(Spearman_Correlation))

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

# just for SRP: go back and remove CDL_WWHT_SOY_DBL_CROP because of 

# create a function to use lapply to plot all 4 univariate sets:

# number<-5

make_plot<-function(number){
  
  OLS<-l.cor[[number]]%>%arrange(desc(Spearman_Correlation))
  
  x<-df.OLS_Sens%>%
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
  
  x
  
}

# use function in lapply (clear plot list first):

# lapply(c(1:5), \(i) make_plot(i))

#












































#### Categorizing land use ####

# USGS criteria:
# Agricultural sites have >50% agricultural land and ≤5% urban land;
# urban sites have >25% urban and ≤25% agricultural land; 
# undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# all other combinations of urban, agricultural, and undeveloped lands are classified as mixed

# initate a dataframe with the landuse values from gauges 2:

df.NLCD06<-distinct(df.OLS_Sens, Name, .keep_all = T)

# I confirmed that pasture + crops = plant:

(df.NLCD06$PASTURENLCD06+df.NLCD06$CROPSNLCD06)==df.NLCD06$PLANTNLCD06

# and that the ones that start with DEV sum up to the first DEV one.

# create the land use class column based on USGS critiera:
# adjusted thresholds:

df.NLCD06<-df.NLCD06%>%
  mutate(USGS.LU.Adjusted = case_when(.default = 'Mixed',
                                      PLANTNLCD06 > 30 & DEVNLCD06 <= 10 ~ 'Agriculture',
                                      DEVNLCD06 > 10 & PLANTNLCD06 <= 30 ~ 'Urban',
                                      DEVNLCD06+PLANTNLCD06 <= 10 ~ 'Undeveloped'))

# frequency table:

table(df.NLCD06$USGS.LU.Adjusted)

# boxplots of

p<-df.NLCD06%>%
  pivot_longer(cols = 2:5, names_to = 'CQ_parameter', values_to = 'Value')%>%
  ggplot(., aes(x=USGS.LU.Adjusted, y=Value, color =USGS.LU.Adjusted ))+
  geom_boxplot(varwidth = TRUE, alpha=0.2)+
  # scale_x_discrete(labels=my_xlab)+
  facet_wrap('CQ_parameter', scales = 'free')+
  stat_compare_means(method = "anova")+      # Add global p-value # , label.y = max(df.datalayers$Value) ???
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5") +
  ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')

# p

#






































#### Grouping CQ curves (stationary, mobilization, dilutionary, complex) ####

# create a list of dataframes for each sites CQ observations:

df.SRP_CQ<-df.SRP_CQ%>%
  rename(Name = site_no)%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))

l.SRP_CQ<-df.SRP_CQ%>%
  split(., .$Name)

# create lm models for each site:

l.SRP_CQ.lm<-lapply(l.SRP_CQ, \(i) lm(log_C~log_Q, data=i))

# save the model coef ad pvals:

SRP_CQ.coef<-as.data.frame(bind_rows(lapply(l.SRP_CQ.lm, \(i) summary(i)$coefficients[,1]), .id = 'site_no'))%>%
  rename(Intercept = 2, Slope = 3)

# save the pvalues 

SRP_CQ.pvals<-as.data.frame(bind_rows(lapply(l.SRP_CQ.lm, \(i) summary(i)$coefficients[,4]), .id = 'site_no'))%>%
  rename(Intercept.pval = 2, Slope.pval = 3)

# merge the two dfs:

df.lm<-left_join(SRP_CQ.coef,SRP_CQ.pvals,by='site_no')

# add column for CQ type:

df.lm<-mutate(df.lm, Type = ifelse(Slope.pval>0.05, 'Stationary', ifelse(Slope>0, 'Mobilization', 'Dilution')))

# merge CQ type labels to the df for plotting

df.SRP_CQ<-left_join(df.SRP_CQ,df.lm%>%select(site_no, Type),by=c('Name'='site_no'))

# create a df for CQplot of all sites with BP analysis:

df_Seg.2<-filter(df_Seg, site %in% df.lm$site_no)%>%
  left_join(.,df.lm%>%select(site_no, Type),by=c('site'='site_no')) 

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

# p

# looking at this plot I want to add a fourth CQ type for complex, if the slopes of the BP analysis look widely different. 
# I will start with calcuating the angle between pre-post BP slope and see which sites get signled out:

# add a new column with the angle between the two lines:

df_Seg.2<-df_Seg.2%>%
  mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))

# I want to add land use as color as well. Merge the plotting df with the land use df:

df_Seg.2<-left_join(df_Seg.2, df.NLCD06%>%select(Name, USGS.LU.Adjusted), by = c('site'='Name'))

# create color pallete for the slope angle:

hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  ylab('log(SRP)')+
  new_scale_color() +
  geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
  geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
  scale_color_manual(name = "Slope Angle", values = hc)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .35)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","yellow", "green"))


p

# export df_Seg.2:

# save(df_Seg.2, file = 'Processed_Data/SRP.df_Seg.2.Rdata')

# # based on this plot, I think I need to zoom in on each one:
# 
# # create a plotting function for each site:
# 
# p.fun<-function(df){
#   ggplot(df, aes(x = log(Q_real), y = log(C)))+
#     geom_point(aes(color = Type))+
#     scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
#     geom_smooth(method = 'lm')+
#     new_scale_color() +
#     geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
#     geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
#     scale_color_manual(name = "Slope Angle", values = hc)+
#     ggtitle(df_Seg.2$site[df_Seg.2$n_sample_rank==i])
# }
# 
# # use function in lappy to make lots of plots (clear plot list first)
# 
# lapply(sort(unique(df_Seg.2$n_sample_rank)), \(i) df_Seg.2%>%filter(n_sample_rank==i)%>%p.fun(.))
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
# 






























#### Steve's Idea: triangle CQ plot ####

# set up dataframe:

# add percent ag and urban to df.lm:

df.tri<-left_join(df.lm, df.NLCD06%>%select(Name, DEVNLCD06, PLANTNLCD06), by = c('site_no'='Name'))%>%select(site_no, Slope, Type, DEVNLCD06, PLANTNLCD06)%>%mutate(Slope = round(Slope, 2), DEVNLCD06=round(DEVNLCD06, 2), PLANTNLCD06=round(PLANTNLCD06, 2))

# create plot:

p<-ggplot(df.tri, aes(x=DEVNLCD06, y=PLANTNLCD06)) +
  geom_point(aes(shape=as.factor(Type), fill=abs(Slope)), size = 4) +
  scale_shape_manual(values = c("Mobilization" = 24, "Dilution" = 25, "Stationary" = 22),
                     guide = guide_legend(override.aes = list(fill = "pink")))+
  scale_fill_gradient(low = "yellow", high = "red")+
  xlab('Percent Developed')+
  ylab('Percent Agriculture')+
  ggtitle(paste(ncode, 'CQ type and slope magnitude as a function of percent Ag and Developed Land'))

p$labels$fill <- "Slope Magnitude"
p$labels$shape <- "CQ Type"

p

# want to recreate this plot in Code/FingerLakesPresentation.R so exporting df.tri:

# save(df.tri, file = 'Processed_Data/df.tri.SRP.Rdata')

#



























#### MLR ####

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
  coef(step.model$finalModel, as.numeric(step.model$bestTune))
  
  # get df of just the predictors in this model:
  
  pred<-names(coef(step.model$finalModel, as.numeric(step.model$bestTune)))[-1]
  
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
  
  p
  
  # look at model residuals:
  
  # plot(m)
  
}

# set m.list names:

names(m.list)<-names(l.cor.MLR.full)

tab_model(m.list, dv.labels = names(m.list), title = paste('Comparison of MLR models for',ncode, 'using forward selection implemented in caret::train'), file="temp.html")

# export m.list for use in Code/FingerLakesPresentation.R:

save(m.list, file = 'Processed_Data/m.list.SRP.Rdata')

#

