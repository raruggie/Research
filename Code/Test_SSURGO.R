library(nhdplusTools)
library(sf)

test<-"04230500"

load(paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', test, '.Rdata'))

mapview(df.sf.soils, zcol = "hydgrpdcd")

i<-which(df.sf.NWIS$Name==test)

df.NWIS.TP_site_metadata<-read.csv("Raw_Data/df.NWIS.TP_site_metadata.csv", colClasses = c(site_no = "character"))%>%
  filter(site_no %in% df.sf.NWIS$Name)

df.site_metadata<-df.NWIS.TP_site_metadata[match(df.sf.NWIS$Name, df.NWIS.TP_site_metadata$site_no),]

long<-df.site_metadata$dec_long_va[i]
lat<-df.site_metadata$dec_lat_va[i]

# workflow for one site:

# run through nhdplusTools workflow:

start_point <- st_sfc(st_point(c(long, lat)), crs = 4269)
start_comid <- discover_nhdplus_id(start_point)

flowline <- navigate_nldi(list(featureSource = "comid", 
                               featureID = start_comid), 
                          mode = "upstreamTributaries", 
                          distance_km = 1000)

subset_file <- tempfile(fileext = ".gpkg")
subset <- subset_nhdplus(comids = as.integer(flowline$UT$nhdplus_comid),
                         # output_file = subset_file,
                         nhdplus_data = "download",
                         flowline_only = FALSE,
                         # overwrite = TRUE,
                         return_data = TRUE)

flowline <- subset$NHDFlowline_Network
catchment <- subset$CatchmentSP
waterbody <- subset$NHDWaterbody


mapview(df.sf.soils, zcol = "hydgrpdcd")+
  # mapview(flowline)+
  mapview(waterbody)
mapview(catchment)



















