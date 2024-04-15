# 

NWIS<-read.csv("Raw_Data/df.NWIS.TP.csv", colClasses = c(site_no = "character"))
WQP<-readWQPqw(siteNumbers = paste0('USGS-', df.NWIS.TP_sites$site_no), parameterCd = pcode)

WQP2<-readWQPqw(siteNumbers = paste0('USGS-', df.NWIS.TP_sites$site_no), 'Total phosphorus')

x <- NWIS %>% group_by(site_no) %>% summarize(n=n())
y <- WQP %>% group_by(MonitoringLocationIdentifier) %>% summarise(nWQP = n()) %>% rename(site_no = 1) %>% mutate(site_no = gsub('USGS-', '', site_no))


z <- left_join(x,y, by = 'site_no')

# same 

dim(NWIS)
dim(WQP)













