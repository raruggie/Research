# 

# Nitrate

NO3 <- read.csv('Downloaded_Data/WQP/NO3.csv')

unique(NO3$ActivityMediaSubdivisionName)

NO3 <- NO3 %>% filter(ActivityMediaSubdivisionName == "Surface Water")

unique(NO3$ResultSampleFractionText)

NO3 <- NO3 %>% filter(ResultSampleFractionText == "Total")

unique(NO3$USGSPCode)

unique(NO3$CharacteristicName)

NO3 <- NO3 %>% mutate(X1 = paste0(ResultSampleFractionText,'-',CharacteristicName,'-',ResultMeasure.MeasureUnitCode))

unique(NO3$X1)

NO3$site_no <- gsub("USGS-", '', NO3$MonitoringLocationIdentifier)

NO3 <- NO3 %>% drop_na(ResultMeasureValue)

NO3 <- NO3 %>% group_by(site_no) %>% mutate(n_samples = n())

sort(unique(NO3$n_samples))

# Phosphorus

# compare WQP P with NWIS query P:

pcode<-'00665'

WQP_P <- read.csv('Downloaded_Data/WQP/P.csv')

WQP_P <- WQP_P %>% filter(ActivityMediaSubdivisionName=="Surface Water")

unique(WQP_P$ResultSampleFractionText)

unique(WQP_P$CharacteristicName)

unique(WQP_P$ResultMeasure.MeasureUnitCode)

WQP_P <- WQP_P %>% mutate(X1 = paste0(ResultSampleFractionText,'-',CharacteristicName,'-',ResultMeasure.MeasureUnitCode))

unique(WQP_P$X1)

unique(WQP_P$USGSPCode)

# TN

TN <- read.csv('Downloaded_Data/WQP/TN.csv')

unique(TN$ActivityMediaSubdivisionName)

TN <- TN %>% filter(ActivityMediaSubdivisionName=="Surface Water")

TN <- TN %>% mutate(X1 = paste0(ResultSampleFractionText,'-',CharacteristicName,'-',ResultMeasure.MeasureUnitCode))

unique(TN$X1)

TN <- TN %>% filter(X1 %in% c("Total-Nitrogen-ug/L", "Total-Nitrogen-mg/L", "Total-Nitrogen-"))

# See if a site is in WQP:

site <- "04231000"

df.site <- NO3 %>% filter(site_no == site)

# there are dissolved samples for this site, but only 38/2 = 19...





