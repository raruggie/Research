# function to turn streamstats polygons into vects
# use 100 for small_artefacts_threshold

SS_sf_to_vect<-function(sf_list, small_artefacts_threshold){
  
  # check validity
  # lapply(sf_list, st_is_valid, reason = TRUE) |> str()
  
  # make valid, check resulting classes
  l2 <- lapply(sf_list, st_make_valid)
  # lapply(l2, st_geometry) |> 
  #   lapply(class) |>
  #   str()
  
  # some shapes ended up as multipolygons, 
  # cast all to polygons and drop small artefacts (100m^2 area threshold)
  l3 <- lapply(l2, st_geometry) |>
    lapply(st_cast, "POLYGON") |> 
    lapply(\(g) g[units::drop_units(st_area(g)) > small_artefacts_threshold])
  # lapply(l3, class) |> str()
  
  # convert to SpatVector, check validity
  l_vect <- lapply(l3, vect)
  # lapply(l_vect, is.valid) |> str()
  
  # lets deal with that as well 
  l_vect <- lapply(l_vect, makeValid)
  
  # check resulting geomtypes
  # lapply(l_vect, geomtype) |> str()
  
}

# not a function but variables used for Reassigning NLCD legend to less classes;

legend<-pal_nlcd()%>%
  mutate(ID2=c(1,2,3,3,3,3,2,4,4,4,4,4,5,5,2,2,6,7,8,8))

legend_2<-data.frame(ID2 = unique(legend$ID2))%>%
  mutate(Class2 = c("Open Water", "Other", "Developed", "Forest", "Grassland", "Pasture/Hay", "Cultivated Crops", "Wetlands"),
         ID3 = c(1,2,3,4,5,6,6,7),
         Class3 = c("Open Water", "Other", "Developed", "Forest", "Grassland", "Ag", "Ag", "Wetlands"))

legend<-left_join(legend, legend_2, by = 'ID2')

legend$Class4<-c('undev','undev', rep('dev',4), rep('undev', 10), 'dev', 'dev', 'undev', 'undev')

legend$Class5<-c('non-Ag','non-Ag', rep('non-Ag',4), rep('non-Ag', 10), 'Ag', 'Ag', 'non-Ag', 'non-Ag')


# function to determine season from date:

getSeason <- function(DATES) {
  WS <- as.Date("2012-12-15", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2012-3-15",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2012-6-15",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2012-9-15",  format = "%Y-%m-%d") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Fall")))
}

# value for scaleing ggplot colors:
cc <- scales::seq_gradient_pal("purple", "yellow", "Lab")(seq(0,1,length.out=9))

# function for elapsed_months

elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

Delineate <- function(long, lat) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      message("Delineating Watershed:")
      
      delineateWatershed(xlocation = long, ylocation = lat, crs = 4326, includeparameters = 'true') 
      # The return value of `streamstats::delineateWatershed` is the actual value 
      # that will be returned in case there is no condition (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message("If no other red text than it worked!")
    }
  )    
  return(out)
}

fun.l.SS_WS.to.sfdf<-function(list_of_SS_WS){
  
  list_of_SS_WS<-list_of_SS_WS[sapply(list_of_SS_WS, function(x) class(x) == "watershed")]
  
  l.sp<-list_of_SS_WS
  
  for (i in seq_along(l.sp)){
    
    tryCatch({
      
      l.sp[[i]]<-toSp(watershed = l.sp[[i]], what = 'boundary')
      
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  # set the names of the list:
  
  names(l.sp)<-names(list_of_SS_WS)
  
  # need to remove the drainage areas that did not work in toSp (nothing we can do about losing these, idk why some dont come out of the function right):
  
  l.sp<-l.sp[!sapply(l.sp, function(x) class(x) == "watershed")]
  
  # convert from sp (old gis in r) to sf (new gis in r)
  
  l.sf<-lapply(l.sp, st_as_sf)
  
  # make valid:
  
  l.sf<-lapply(l.sf, st_make_valid)
  
  # need to convert the Shape_Area column to numeric for all dfsin the list or bind_rows wont work:
  
  l.sf<-lapply(l.sf, \(i) i%>%mutate(Shape_Area = as.numeric(Shape_Area)))
  
  # create a single sf df with all the sample site draiange areas:
  
  df.sf<-bind_rows(l.sf, .id = 'Name')%>%
    relocate(Name, .before= 1)
  
  return(df.sf)
}

Ryan_toSp <- function(DA) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      message("converting to sp:")
      
      toSp(DA) 
      # The return value of `streamstats::delineateWatershed` is the actual value 
      # that will be returned in case there is no condition (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message("If no other red text than it worked!")
    }
  )    
  return(out)
}


# Snap a point to the closest point on a line segment using sf
# https://stackoverflow.com/questions/51292952/snap-a-point-to-the-closest-point-on-a-line-segment-using-sf

st_snap_points = function(x, y, max_dist = 1000) {
  
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  return(out)
}

# a function is created to pull the data for each consituent this will be used in lapply to create a list with dataframes for  each consituent
# whatNWISdata is again used but this time the pcode for the different nutrients are inputs to the function created below:

# test variables for function building:

# parameter_code<-'00600' # pcode for TP
# 
# df.flow_query<-df.NWIS.Q_sites
# 
# state<-'NY'

fun.df.Pair_consit_flow<-function(parameter_code, df.flow_query, n_samples = 20, state = 'NY') {
  
  # download metadata for the sites in NY with data of the constituent of interest:
  
  df.consit_fun<-whatNWISdata(stateCd=state,parameterCd=parameter_code)
  
  # reduce duplicate site entries:
  
  df.consit_fun<-df.consit_fun%>% 
    group_by(site_no)%>% 
    slice(which.max(count_nu)) 
  
  # make sure the site_no columns are the same class:
  
  df.flow_query<-df.flow_query%>%mutate(site_no = as.character(site_no))
  
  # pair up this dataframe with the flow query dataframe:
  # select only asubsetof columns
  # filter to the sites returned in the flow query
  # then left join with the flow query df (could have just done a right join instead of filter+left_join)
  # filter for sites with over 99 samples
  # arrange the resulting df by the number of samples and the number of flow observations:
  
  df.consit_flow_fun<-df.consit_fun%>%
    select(c(2,3,5,6,22,23,24))%>% # "site_no"     "station_nm"  "dec_lat_va"  "dec_long_va" "begin_date"  "end_date"    "count_nu"
    filter(site_no %in% df.flow_query$site_no)%>%
    left_join(.,df.flow_query[,c(2,3,5,6,22,23,24)])%>%
    filter(count_nu > n_samples)%>%
    arrange(desc(count_nu), desc(nflowdays))
  
  # return this dataframe as the result of this function
  return(df.consit_flow_fun)
  

}

# read in allsheets of a excel workbook into a list (used it for gauges2)

library(readxl)    

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# function to take dataframe(s) of 21 day climate lags and turn into a single 1 row df of the processed climate variables:

fun.Process_climate<-function(df, site, date){
  #make single row dataframes for the processed climate variables:
  # precip
  df.rain<-data.frame(Name = rain_colnames, Value = round(cumsum(df$pr)[rain_cum_and_temp_deltas],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmax lag
  df.tmax.lag<-data.frame(Name = tmax_lag_colnames, Value = round(df$tmmx[temp_lags+1],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmax delta
  df.tmax.delta<-data.frame(Name = tmax_delta_colnames, Value = round(df$tmmx[rain_cum_and_temp_deltas+1]-df$tmmx[1],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmin lag
  df.tmin.lag<-data.frame(Name = tmin_lag_colnames, Value = df$tmmn[temp_lags+1])%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmin delta
  df.tmin.delta<-data.frame(Name = tmin_delta_colnames, Value = round(df$tmmn[rain_cum_and_temp_deltas+1]-df$tmmn[1],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # combine these dataframes and add a site and date column:
  df.date<-bind_cols(df.rain, df.tmax.lag, df.tmax.delta, df.tmin.lag, df.tmin.delta)%>%
    mutate(site_no = vect.NWIS$Name[i], Date = date, .before = 1)
  return(df.date)
}

# normalize from 0 to 1:

normalized<-function(y) {
  
  x<-y[!is.na(y)]
  
  x<-(x - min(x)) / (max(x) - min(x))
  
  y[!is.na(y)]<-x
  
  return(y)
}


add_columns <- function(df, columns){
  new <- rep(NA_character_, length(columns))
  names(new) <- columns
  mutate(df, !!!new)
}



# get watershed percent HSG for a single site:
# *use this function in lapply for multiple sites:
# if function doesnt work (maybe throws error that says file not found), try emptying Users/ryrug/AppData/Local/Temp folder
# 

# i<-12
# site_no<-df.sf.NWIS$Name[i]
# sf.df<-df.sf.NWIS

fun.SURRGO_HSG<-function(site_no, sf.df){
  
  print(site_no)
  
  # workflow for one site:
  
  # set function variables:
  
  template<-sf.df$geometry[sf.df$Name==site_no]
  
  # download soil data for site:
  # you may get: 
  # Error:There are no complete soil surveys in your study area.
  # so using tryCatch here:
  
  tryCatch({
    
    l.soils<-FedData::get_ssurgo(template = template, label = site_no)
    # look at map:
    
    # mapview(l.soils[[1]])+mapview(template)
    
    # will need to clip but can do that later...
    
    # the mukey in the spatial element can be used as a joining column for the HSG, which is located in l.soils$tabular$muaggatt:
    
    df.sf.soils<-left_join(l.soils[[1]], l.soils$tabular$muaggatt%>%select(mukey, hydgrpdcd)%>%mutate(MUKEY = as.character(mukey), .keep = 'unused'), by = 'MUKEY')
    
    save(df.sf.soils, file = paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', site_no, '.Rdata'))
    
    # load(paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', site_no, '.Rdata'))
    
    # lets filter the MU that are NA and replot:
    
    # df.sf.soils%>%
    #   # drop_na(hydgrpdcd)%>%
    #   mapview(., zcol = 'hydgrpdcd')+mapview(template)
    
    # map of HSG A and the rest:
    
    # df.sf.soils%>%
    #   # drop_na(hydgrpdcd)%>%
    #   mutate('hydgrpdcd'=ifelse(hydgrpdcd=='A', 'A', 'Else'))%>%
    #   mapview(., zcol = 'hydgrpdcd')
    
    # it looks like the NAs are water!!
    
    # convert to SpatVector to perform spatial analysis:
    
    vect.soils<-vect(df.sf.soils) # %>%drop_na(hydgrpdcd))
    
    vect.DA.template<-vect(template)
    
    # calcualte watershed area in unit ha:
    
    site_DA.ha<-round(expanse(vect.DA.template, unit="ha", transform=TRUE), 2)
    
    # crop soils vector to draiange area vector:
    
    vect.soils<-terra::crop(vect.soils, vect.DA.template)
    
    # look at map:
    
    plot(x= vect.soils, y='hydgrpdcd')
    
    # add column of each polygons area:
    
    vect.soils$area_ha<-expanse(x= vect.soils, unit = 'ha')
    
    # make a dataframe of the areas and HSG, and calcualte the percent of each HSG:
    
    df.HSG<-data.frame(Area = vect.soils$area_ha, HSG = vect.soils$hydgrpdcd)%>%
      drop_na(HSG)%>%
      group_by(HSG)%>%
      summarize(Watershed_Percent = round(sum(Area)/site_DA.ha, 4))
    
    return(df.HSG)
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# same function as above but for use when sites soil map has already been downloaded:

fun.SURRGO_HSG.already_downloaded<-function(site_no, sf.df){
  
  print(site_no)
  print(which(sf.df$Name==site_no))
  
  # workflow for one site:
  
  # set function variables:
  
  template<-sf.df$geometry[sf.df$Name==site_no]
  
  # download soil data for site:
  
  # l.soils<-FedData::get_ssurgo(template = template, label = site_no) 
  
  # look at map:
  
  # mapview(l.soils[[1]])+mapview(template)
  
  # will need to clip but can do that later...
  
  # the mukey in the spatial element can be used as a joining column for the HSG, which is located in l.soils$tabular$muaggatt:
  
  # df.sf.soils<-left_join(l.soils[[1]], l.soils$tabular$muaggatt%>%select(mukey, hydgrpdcd)%>%mutate(MUKEY = as.character(mukey), .keep = 'unused'), by = 'MUKEY')
  
  # save(df.sf.soils, file = paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', site_no, '.Rdata'))
  
  load(paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', site_no, '.Rdata'))
  
  # lets filter the MU that are NA and replot:
  
  # df.sf.soils%>%
  #   # drop_na(hydgrpdcd)%>%
  #   mapview(., zcol = 'hydgrpdcd')+mapview(template)
  
  # map of HSG A and the rest:
  
  # df.sf.soils%>%
  #   # drop_na(hydgrpdcd)%>%
  #   mutate('hydgrpdcd'=ifelse(hydgrpdcd=='A', 'A', 'Else'))%>%
  #   mapview(., zcol = 'hydgrpdcd')
  
  # it looks like the NAs are water!!
  
  # convert to SpatVector to perform spatial analysis:
  
  vect.soils<-vect(df.sf.soils) # %>%drop_na(hydgrpdcd))
  
  vect.DA.template<-vect(template)
  
  # calcualte watershed area in unit ha:
  
  site_DA.ha<-round(expanse(vect.DA.template, unit="ha", transform=TRUE), 2)
  
  # crop soils vector to draiange area vector:
  
  vect.soils<-terra::crop(vect.soils, vect.DA.template)
  
  # look at map:
  
  # plot(x= vect.soils, y='hydgrpdcd')
  
  # add column of each polygons area:
  
  vect.soils$area_ha<-expanse(x= vect.soils, unit = 'ha')
  
  # make a dataframe of the areas and HSG, and calcualte the percent of each HSG:
  
  df.HSG<-data.frame(Area = vect.soils$area_ha, HSG = vect.soils$hydgrpdcd)%>%
    drop_na(HSG)%>%
    group_by(HSG)%>%
    summarize(Watershed_Percent = round(sum(Area)/site_DA.ha, 4))
  
  return(df.HSG)
}


# function to calculate 100 and 800 meter riparian buffer percent land use:

i<-19

# mapview(df.sf.NWIS[i,])

fun.NHD.RIP.buffers<-function(i, df.site_metadata, save.flowline = TRUE){
  
  print(i)
  
  # set up function inputs:
  
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
                           flowline_only = TRUE, # if false you will get catchment + waterbody
                           # overwrite = TRUE,
                           return_data = TRUE)

  flowline <- subset$NHDFlowline_Network
  # catchment <- subset$CatchmentSP
  # waterbody <- subset$NHDWaterbody
  
  ## Or cando this by reading in from file:
  
  # flowline <- sf::read_sf(subset_file, "NHDFlowline_Network")
  # catchment <- sf::read_sf(subset_file, "CatchmentSP")
  # waterbody <- sf::read_sf(subset_file, "NHDWaterbody")
  
  # save sites flowline to file:
  
  if(save.flowline){
    
    save(flowline, file = paste('Downloaded_Data/NHD_flowlines/flowline_', df.sf.NWIS$Name[i], '.Rdata'))
    
    print('saved flowline as sf')
  } else {print('flow saving turned off')}
  
  ## now plot:
  
  # plot(sf::st_geometry(flowline), col = "blue")
  # plot(start_point, cex = 1.5, lwd = 2, col = "red", add = TRUE)
  # plot(sf::st_geometry(catchment), add = TRUE)
  # plot(sf::st_geometry(waterbody), col = rgb(0, 0, 1, alpha = 0.5), add = TRUE)

  # this is so cool!!
  
  # play around with flowline df:
  
  # mapview(flowline, zcol = "ftype")
  
  # my goal with the NHD is to use it to calcualte the landuse in the buffer,
  # so I really only need the StreamRiver ftype. 
  
  flowline1<-flowline%>%filter(ftype == 'StreamRiver')
  
  # create 100 and 800 meter buffer. use st_union to remove overlap:
  
  buffer.100<-st_buffer(flowline1, 100)%>%st_union()
  
  buffer.800<-st_buffer(flowline1, 800)%>%st_union()
  
  # Look at map:
  
  # mapview(buffer.100)+mapview(flowline1)
  
  # looks good
  
  # calculate the watershed 100 and 800 meter percent in major landuse types. to do this:
  
  # Load in the NLCD for the watershed:
  
  rast.NLCD<-rast(names.2006[i])
  
  # reproject to buffer vector data to match raster data:
  
  vect.buffer.100<-vect(buffer.100) 
  vect.buffer.800<-vect(buffer.800)
  
  vect.buffer.100.proj<-terra::project(vect.buffer.100, crs(rast.NLCD))
  vect.buffer.800.proj<-terra::project(vect.buffer.800, crs(rast.NLCD))
  
  # extract frequency tables:
  
  df.buffer.freq.100 <- terra::extract(rast.NLCD, vect.buffer.100.proj, ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2))
  df.buffer.freq.800 <- terra::extract(rast.NLCD, vect.buffer.800.proj, ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2))
  
  # makesure each NLCD class is represented. todo this: left join with legend:
  
  df.buffer.freq.100<-left_join(legend%>%select(Class), df.buffer.freq.100, by = 'Class')%>%replace(is.na(.), 0)
  df.buffer.freq.800<-left_join(legend%>%select(Class), df.buffer.freq.800, by = 'Class')%>%replace(is.na(.), 0)
  
  # reclassify: the GAGES II predictors for NLCD land use are the sum of a few NLCD classes (see the kable table of the variable descriptions made above). To do this:
  
  # create and ordered vector based on legend (legend comes from Ryan_funcitons.R)
  
  Class3.for.G2<-c("WATERNLCD06", "SNOWICENLCD06", rep("DEVNLCD06", 4), "BARRENNLCD06", "DECIDNLCD06", "EVERGRNLCD06", "MIXEDFORNLCD06", NA, "SHRUBNLCD06", "GRASSNLCD06", NA, NA, NA, "PASTURENLCD06", "CROPSNLCD06", "WOODYWETNLCD06", "EMERGWETNLCD06")
  
  # add an identifier to this vector so the column names are slightly different than the GAGES II predictors:
  
  Class3.for.G2<-paste0('R_', Class3.for.G2)
  
  # create a df from this vector (for latter use):
  
  df.Class3<-data.frame(Class = unique(Class3.for.G2)[complete.cases(unique(Class3.for.G2))])
  
  # create new df from the legend dataframe:
  
  legend.for.G2<-legend%>%mutate(Class3 = Class3.for.G2)
  
  # reclassify the NLCD using this new legend and clean up the dataframe:
  
  df.buffer.freq.100<-left_join(df.buffer.freq.100, legend.for.G2%>%select(Class, Class3), by = 'Class')%>%
    mutate(Class = Class3)%>%
    select(-Class3)%>%
    group_by(Class)%>%
    summarise(Freq = sum(Freq))%>%
    pivot_wider(names_from = Class, values_from = Freq)%>%
    mutate(Name = df.sf.NWIS$Name[i], .before = 1)%>%
    mutate(R_FORESTNLCD06 = R_DECIDNLCD06+R_EVERGRNLCD06+R_MIXEDFORNLCD06,
           R_PLANTNLCD06 = R_PASTURENLCD06+R_CROPSNLCD06)%>%
    pivot_longer(cols = -Name)
  
  df.buffer.freq.800<-left_join(df.buffer.freq.800, legend.for.G2%>%select(Class, Class3), by = 'Class')%>%
    mutate(Class = Class3)%>%
    select(-Class3)%>%
    group_by(Class)%>%
    summarise(Freq = sum(Freq))%>%
    pivot_wider(names_from = Class, values_from = Freq)%>%
    mutate(Name = df.sf.NWIS$Name[i], .before = 1)%>%
    mutate(R_FORESTNLCD06 = R_DECIDNLCD06+R_EVERGRNLCD06+R_MIXEDFORNLCD06,
           R_PLANTNLCD06 = R_PASTURENLCD06+R_CROPSNLCD06)%>%
    pivot_longer(cols = -Name)
  
  # I really only need DEV, PLANT, and FOREST.
  # so filtering these, pivoting wider, and merging with G2:
  
  df.buffer.freq.100<-df.buffer.freq.100%>%
    filter(name %in% c('R_DEVNLCD06', 'R_PLANTNLCD06', 'R_FORESTNLCD06'))%>%
    pivot_wider(names_from = name, values_from = value)%>%
    rename_with(~paste0(., "_100"), starts_with("R_"))
  
  df.buffer.freq.800<-df.buffer.freq.800%>%
    filter(name %in% c('R_DEVNLCD06', 'R_PLANTNLCD06', 'R_FORESTNLCD06'))%>%
    pivot_wider(names_from = name, values_from = value)%>%
    rename_with(~paste0(., "_800"), starts_with("R_"))
  
  # create a df of the 100 and 800 meter ripiran predictors:
  
  df.RIP<-left_join(df.buffer.freq.100, df.buffer.freq.800, by = 'Name')
  
  # rename columns to G2, but paste a unique ID to it:
  
  names(df.RIP)<-c('Name', paste0('R_', names(df.G2.reduced%>%select(contains('RIP')))))
  
  return(df.RIP)
}

# function to calcualte FRAGUN_BASIN:

i<-1

# landscape<-l.rast.NWIS.NLCD.2001[[i]]

fun.FRAGUN_BASIN<-function(i, landscape){
  
  print(i)
  
  # check to make sure landscape raster can be used in landscapemetrics functions:
  
  print(check_landscape(landscape))
  
  # crop landscape raster to watershed boundary: to do this:
  
  # load in wtershed boundary as Spatvector:
  
  vect.boundary<-vect(df.sf.NWIS[i])
  
  # reproject vector to match raster:
  
  vect.boundary<-terra::project(vect.boundary, crs(landscape))
  
  # crop:
  
  landscape<-crop(landscape, vect.boundary, mask=TRUE)
  
  # calcualte core area percentage of landscape (class level):
  
  cpland<-lsm_c_cpland(landscape,directions = 8,consider_boundary = FALSE,edge_depth = 1)
  
  # merge legend Class 4 to determine undeveloped classes:
  
  cpland.undev<-left_join(cpland, legend%>%select(ID, Class4), by = c('class'= 'ID'))%>%filter(Class4 == 'undev')
  
  # take 100-sum of column to get FRAGUN BASIN:
  
  FRAGUN_BASIN<-100-sum(cpland.undev$value)
  
  return(FRAGUN_BASIN)
  
}

# function to determine percent CSA in watershed:

i<-5

# WA<-df.sf.NWIS[i,]

fun.CSA_watershed_percent<-function(i, WA){
  
  print(i)
  
  tryCatch({
    
    # mapview(WA)
    
    # Step 1: create slope map.to do this:
    
    # convert watershed boundary to spatvector:
    
    vect.boundary<-vect(WA)
    
    # reproject to match DEM raster crs:
    
    vect.boundary<-project(vect.boundary, crs(DEM.NWIS))
    
    # crop NED:
    
    rast.NED<-crop(DEM.NWIS, vect.boundary, mask=T)
    
    # create slope map:
    
    rast.slope<-terrain(rast.NED, "slope")
    
    # bin slope into three categories:
    
    m<-c(0,5,1,5,10,2,10,100,3)
    rclmat<-matrix(m, ncol = 3, byrow = T)
    rast.slope<-classify(rast.slope,rclmat,include.lowest=T)
    
    # rename raster levels:
    
    rast.slope<-as.factor(rast.slope)
    
    levels(rast.slope)<-data.frame(ID=1:3,Slope=c('<5%','5-10%','>10%'))
    
    # Step 2: soils:
    
    # load in HSG map:
    
    load(paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', df.sf.NWIS$Name[i], '.Rdata'))
    
    # convert to spatvector:
    
    vect.soils<-vect(df.sf.soils%>%drop_na(hydgrpdcd))
    
    # replace X/D HSG with just D:
    
    vect.soils$hydgrpdcd<-gsub(".*/.*", "D", vect.soils$hydgrpdcd)
    
    # rasterize:
    
    rast.soils<-rasterize(vect.soils, rast.slope, field = 'hydgrpdcd')
    
    # classify levels to only 3 (will be used later for good, fair, and poor drainage):
    
    m<-c(0,1,1,2,2,3,3,3)
    rclmat<-matrix(m, ncol = 2, byrow = T)
    rast.soils<-classify(rast.soils,rclmat)
    
    # set levels to a name:
    
    rast.soils<-as.factor(rast.soils)
    levels(rast.soils)<-data.frame(ID = 1:3, Drainage = c('Good', 'Fair', 'Poor'))
    
    # Step 3: Land use:
    
    # load landuse raster:
    
    rast.LU<-rast(names.2006[i]) # no need to crop for now...
    
    # set levels:
    
    levels(rast.LU)<-levels(rast.LU)[[1]]%>%drop_na(Class)%>%left_join(., legend%>%select(ID, Class4), by = 'ID')%>%select(-Class)%>%rename(Dev = Class4)
    
    # classify down to just two levels:
    
    m<-c(legend$ID, ifelse(legend$Class4=='dev', 1, 0))
    rclmat<-matrix(m, ncol = 2, byrow = F)
    rast.LU<-classify(rast.LU,rclmat)
    rast.LU<-as.factor(rast.LU)
    levels(rast.LU)<-data.frame(ID = 0:1, Dev = c('undev', 'dev'))
    
    # Step 4: combine into single map:
    
    # first check crs:
    
    # crs(rast.slope)
    # crs(rast.soils)
    # crs(rast.LU)
    
    # reproject rast.LU:
    
    rast.LU<-project(rast.LU, rast.soils)
    
    # now merge maps:
    
    CSA.map<-concats(rast.slope, rast.soils)
    
    CSA.map<-concats(CSA.map, rast.LU)
    
    # turn into binary CSA map:
    
    levels(CSA.map)<-data.frame(ID = 0:17, CSA = c(rep('non-CSA',5), 'CSA', rep('non-CSA', 3), 'CSA', 'non-CSA', 'CSA', 'non-CSA', 'CSA', 'non-CSA', 'CSA', 'non-CSA', 'CSA'))
    m<-c(0:17,ifelse(levels(CSA.map)[[1]]$CSA=='CSA', 1,0))
    rclmat<-matrix(m, ncol = 2, byrow = F)
    CSA.map<-classify(CSA.map,rclmat)
    CSA.map<-as.factor(CSA.map)
    levels(CSA.map)<-data.frame(ID = 0:1, CSA = c('non-CSA', 'CSA'))
    
    # save CSA map at this stage:
    
    writeRaster(CSA.map, file = paste0('Processed_Data/CSA_maps/CSA_map_', df.sf.NWIS$Name[i], '.tif'), overwrite=TRUE)
    
    # sum up CSA area for watershed:
    
    df.CSA<-freq(CSA.map)
    
    # determine percent watershed:
    
    df.CSA$perc <- round(100 * df.CSA$count / sum(df.CSA$count), 1)
    
    percent.CSA<-df.CSA$perc[df.CSA$value=='CSA']
    
    return(percent.CSA)
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
}

# function to determine percent riparian buffer CSA in watershed:

i<-8

# sf.df<-df.sf.NWIS.62

fun.Riparian_CSA<-function(i, sf.df){
  
  # set funciton variable:
  
  WA<-sf.df[i,]
  
  print(i)
  print(WA$Name)
  
  tryCatch({
    
    # mapview(WA)
    
    # Step 1: load in sites flowline:
    
    load(paste0('Downloaded_Data/NHD_flowlines/flowline_', sf.df$Name[i], '.Rdata'))
    
    # Step 2: load in sites CSA map:
    
    CSA.map<-rast(paste0('Processed_Data/CSA_maps/CSA_map_', sf.df$Name[i], '.tif'))
    
    # Step 3: extract CSA map over ripiran buffer: todothis:
    
    # follow buffer workflow in fun.NHD.RIP.buffers:
    
    flowline1<-flowline%>%filter(ftype == 'StreamRiver')
    buffer.100<-st_buffer(flowline1, 100)%>%st_union()
    buffer.800<-st_buffer(flowline1, 800)%>%st_union()
    vect.buffer.100<-vect(buffer.100) 
    vect.buffer.800<-vect(buffer.800)
    vect.buffer.100.proj<-terra::project(vect.buffer.100, crs(CSA.map))
    vect.buffer.800.proj<-terra::project(vect.buffer.800, crs(CSA.map))
    
    # mapview(WA)+mapview(flowline)
    
    # extract CSA.map over buffer:
    
    CSA.map.100<-terra::extract(CSA.map, vect.buffer.100.proj, ID=FALSE)%>%drop_na(CSA)
    CSA.map.800<-terra::extract(CSA.map, vect.buffer.800.proj, ID=FALSE)%>%drop_na(CSA)
    
    # calculate proportion of CSA pixels:
    
    v.100<-(sum(CSA.map.100$CSA=='CSA')/length(CSA.map.100$CSA))*100
    v.800<-(sum(CSA.map.800$CSA=='CSA')/length(CSA.map.800$CSA))*100

    # create df to return from function:
    
    df.RIP.CSA<-data.frame(Name = WA$Name, RIP.CSA.100 = v.100, RIP.CSA.800 = v.800)

    # return:
    
    return(df.RIP.CSA)
    
    #
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
}

# function for calcualting watershed percent Ag in CSA:

i<-2

# WA<-df.sf.NWIS[i,]

fun.Ag.CSA<-function(i, WA){
  
  print(i)
  print(df.sf.NWIS$Name[i])
  
  tryCatch({
    
    # workflow:
    
    # create categoical slope raster:
    
    vect.boundary<-vect(WA)
    vect.boundary<-project(vect.boundary, crs(DEM.NWIS))
    rast.NED<-crop(DEM.NWIS, vect.boundary, mask=T)
    rast.slope<-terrain(rast.NED, "slope")
    # plot(rast.slope)
    
    # create LU map and set levels to just non-Ag and Ag:
    
    rast.LU<-rast(names.2016[i])
    vect.DA<-vect(WA)
    vect.DA<-project(vect.DA, crs(rast.LU))
    rast.LU<-crop(rast.LU, vect.DA, mask=TRUE)
    rast.LU<-project(rast.LU, rast.slope)
    levels(rast.LU)<-levels(rast.LU)[[1]]%>%drop_na(Class)%>%left_join(., legend%>%select(ID, Class5), by = 'ID')%>%select(-Class)%>%rename(Ag = Class5)
    m<-c(legend$ID, ifelse(legend$Class5=='Ag', 1, 0)) 
    rclmat<-matrix(m, ncol = 2, byrow = F)
    rast.LU<-classify(rast.LU,rclmat)
    rast.LU<-as.factor(rast.LU)
    levels(rast.LU)<-data.frame(ID = 0:1, Ag = c('non-Ag', 'Ag'))
    # plot(rast.LU)
    
    # create HSG raster and classify to two levels, good (0) and poor (1) soils:
    
    load(paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', df.sf.NWIS$Name[i], '.Rdata'))
    vect.soils<-vect(df.sf.soils%>%drop_na(hydgrpdcd))
    vect.soils$hydgrpdcd<-gsub(".*/.*", "D", vect.soils$hydgrpdcd)
    rast.soils<-rasterize(vect.soils, rast.slope, field = 'hydgrpdcd')
    m<-c(0,0,1,0,2,1,3,1)
    rclmat<-matrix(m, ncol = 2, byrow = T)
    rast.soils<-classify(rast.soils,rclmat)
    rast.soils<-as.factor(rast.soils)
    # plot(rast.soils)
    
    # test: Create a raster of just slope where there is land use Ag:
    
    Ag <- ifel(rast.LU == 'Ag', 1, NA)
    # plot(Ag)
    slope.Ag <- mask(rast.slope, Ag)
    
    # plot:
    
    # plot(slope.Ag)
    # hist(slope.Ag)
    # max(slope.Ag)
    
    # it looks like most of the Ag land is below 5%, but goes 
    # as high as 24% (that could be bad data too...)
    
    # now calculate watershed percent Ag with over 5% slope + HSG C/D. to do this:
    
    # create steep slope, Ag, and poor draiange rasters but use 0 instead of NA in ifel:
    
    Ag<-as.factor(Ag <- ifel(rast.LU == 'Ag', 1, 0))
    # plot(Ag)
    Steep<-as.factor(ifel(rast.slope>5, 1,0))
    # plot(Steep)
    Soggy<-rast.soils
    # plot(Soggy)
    
    # merge rasters:
    
    Steep.Wet.Ag<-concats(Steep, Soggy)
    Steep.Wet.Ag<-concats(Steep.Wet.Ag, Ag)
    
    # plot:
    
    # plot(Steep.Wet.Ag)
    # levels(Steep.Wet.Ag)
    
    # calculate frequency table (1_1_1 is Ag CSA:)
    
    CSA.Ag<-freq(Steep.Wet.Ag)
    
    # determine percent watershed:
    
    CSA.Ag$perc <- round(100 * CSA.Ag$count / sum(CSA.Ag$count), 10)
    
    # extract out Ag CSA:
    
    Ag.CSA<-round(CSA.Ag$perc[CSA.Ag$value=='1_1_1'], 3)
    
    return(Ag.CSA)
    
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  
}

# function for calcuating breakpoint in CQ curve:

i<-1

# df.CQ<-df.TP_CQ


fun.CQ.BP<-function(i, df.CQ){
  
  site<-sort(unique(df.CQ$site_no))[i]
  
  # workflow:
  
  tryCatch({
    
    # print the site name for loop debugging:
    
    print(i)
    print(site)

    
    # create a dataframe that will work with segmented. To do this: 
    # filter for the site the loop is in
    # add log transformed C and Q columns, as well as duplicated columns for renamed C and Q
    # filter for real log C and Q values so breakpoint analysis works smoothly:
    
    df<-df.CQ%>%
      filter(site_no == site)%>%
      mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
      filter(is.finite(log_C))%>%
      filter(is.finite(log_Q))
    
    # build a single slope lm for log C and Q. Tis model is also used inthe breakpoint analysis inthenext step:
    
    m<-lm(log_C~log_Q, df)
    
    # perform breakpoint regression:
    
    m_seg<-segmented(obj = m, npsi = 1)
    
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
    
    # get the model fitted data and put in a dataframe:
    
    fit <- data.frame(Q = df$log_Q, Seg_C = fitted(m_seg))
    
    # reformat this dataframe to export out of the loop
    
    if(length(class(m_seg))==2){
      result_df<-fit%>%mutate(site = site, Date = df$sample_dt, n = df$n, Q_real = df$Q, C = df$C, I1 = inter$Est.[1], I2 = inter$Est.[2], Slope1 = s[1,1], Slope2 = s[2,1], BP = bp)
    } else{
      result_df<-fit%>%mutate(site = site, Date = df$sample_dt, n = df$n, Q_real = df$Q, C = df$C, I1 = NA, I2 = NA, Slope1 = NA, Slope2 = NA, BP = NA)
    }
    
    # set BP_yes column:
    
    result_df$BP_yes<-davies.test(m)$p.val<0.05
    
    # return:
    
    return(result_df)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# function for WRTDS:

fun.R.modelEstimation<-function(eList){
  
  tryCatch({
    
    eList<-modelEstimation(eList)
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  return(eList)
  
  
}

# funtion to extract lm from results of train w/ forward selection:

# model <- l.step.model[[2]]
# 
# x <- model$finalModel
# 
# x
# 
# as.numeric(model$bestTune)
# 
# coef(x, as.numeric(model$bestTune)-2)
# 
# df <- l.setup.MLR[[i]]

fun.train.MLR.to.lm <- function(model, df){
  
  # get df of just the predictors in this model:
  
  pred<-names(coef(model$finalModel, as.numeric(model$bestTune)))[-1]
  
  # remove potential backticks:
  
  pred<-gsub("`", "", pred)
  
  # select down df:
  
  df.2<-df%>%select(term, pred)
  
  # make lm:
  
  m<-lm(term~., data=df.2)
  
  # Return:
  
  return(m)
  
}

# function to extract significant model terms:

# model <- l.m.final[[i]]

fun.extract.sig.model.terms <- function(model){
  
  v.sig <- coef(summary(model))[,4][-1]
  
  v.sig.names <- names(v.sig[v.sig<0.05])
  
  return(v.sig.names)
  
}

# function to compare full lm, regsbsest, and stepAIC models:

i <- 6

# df <- l.setup.MLR[[i]]

fun.compare.linear.model.selections <- function(df, nbest = 1){
  
  # Step 1):
  
  # build full model:
  
  m.1 <- lm(term~., df)
  
  # extract signficant predictors name and coef:
  
  coef.pvals <- coef(summary(m.1))[,4][-1]
  v.terms.coef.names <- names(coef.pvals[coef.pvals<0.05])
  v.terms.coef <- coef(summary(m.1))[,1][-1][coef(summary(m.1))[,4][-1]<0.05]
  df.coef.1 <- data.frame(name=v.terms.coef.names,value=v.terms.coef)
  rownames(df.coef.1) <- NULL
  
  # extract adjR2 and vif:
  
  if(dim(df.coef.1)[1]!=0){
    df.coef.1$adjR2 <- summary(m.1)$adj.r.squared
    ifelse(length(m.1$coefficients)>2, df.vif <- car::vif(m.1) %>% as.data.frame() %>% tibble::rownames_to_column(., "name") %>% rename(vif=2) %>% mutate(name = gsub("`", "", name)), df.vif <- data.frame(name=names(m.1$coefficients[2]), vif=NA))
    # df.vif <- car::vif(m.1) %>% as.data.frame() %>% tibble::rownames_to_column(., "name") %>% rename(vif=2) %>% mutate(name = gsub("`", "", name))
    df.coef.1 <- left_join(df.coef.1, df.vif, by = 'name')
  }
  
  
  # Step 2)
  
  # build regsubset model:
  
  Best_Subset <- regsubsets(term~., data =df,
                            nbest = 1,      # 1 best model for each number of predictors
                            nvmax = 6,    # NULL for no limit on number of variables
                            force.in = NULL, 
                            force.out = NULL,
                            method = "exhaustive")
  
  # find best model:
  
  highest.adjR2 <- which.max(summary(Best_Subset)$adjr2)
  
  # extract predictors from best model:
  
  preds <- names(coef(Best_Subset, highest.adjR2))[-1]
  
  # remove potential backticks:
  
  preds<-gsub("`", "", preds)
  
  # build lm with best predictors:
  
  m.2 <- lm(term ~ ., data = df %>% select(term, preds))
  
  # extract signficant predictors name and coef:
  
  coef.pvals <- coef(summary(m.2))[,4][-1]
  v.terms.coef.names <- names(coef.pvals[coef.pvals<0.05])
  v.terms.coef <- coef(summary(m.2))[,1][-1][coef(summary(m.2))[,4][-1]<0.05]
  df.coef.2 <- data.frame(name=v.terms.coef.names,value=v.terms.coef)
  rownames(df.coef.2) <- NULL
  
  # extract adjR2:
  
  if(dim(df.coef.2)[1]!=0){
    df.coef.2$adjR2 <- summary(m.2)$adj.r.squared
    ifelse(length(m.2$coefficients)>2, df.vif <- car::vif(m.2) %>% as.data.frame() %>% tibble::rownames_to_column(., "name") %>% rename(vif=2) %>% mutate(name = gsub("`", "", name)), df.vif <- data.frame(name=names(m.2$coefficients[2]), vif=NA))
    # df.vif <- car::vif(m.2) %>% as.data.frame() %>% tibble::rownames_to_column(., "name") %>% rename(vif=2) %>% mutate(name = gsub("`", "", name))
    df.coef.2 <- left_join(df.coef.2, df.vif, by = 'name')
  }
  
  # Step 3)
  
  # build stepAIC model:
  
  m.3 <- stepAIC(m.1, direction = "both", trace = FALSE)
  
  # extract signficant predictors name and coef:
  
  coef.pvals <- coef(summary(m.3))[,4][-1]
  v.terms.coef.names <- names(coef.pvals[coef.pvals<0.05])
  v.terms.coef <- coef(summary(m.3))[,1][-1][coef(summary(m.3))[,4][-1]<0.05]
  df.coef.3 <- data.frame(name=v.terms.coef.names,value=v.terms.coef)
  rownames(df.coef.3) <- NULL
  
  # extract adjR2:
  
  if(dim(df.coef.3)[1]!=0){
    df.coef.3$adjR2 <- summary(m.3)$adj.r.squared
    ifelse(length(m.3$coefficients)>2, df.vif <- car::vif(m.3) %>% as.data.frame() %>% tibble::rownames_to_column(., "name") %>% rename(vif=2) %>% mutate(name = gsub("`", "", name)), df.vif <- data.frame(name=names(m.3$coefficients[2]), vif=NA))
    df.coef.3 <- left_join(df.coef.3, df.vif, by = 'name')
  }
  
  # Step 4) 
  
  # combine df.coef's into single df:
  
  df.coef <- list(df.coef.1, df.coef.2, df.coef.3) %>% 
    purrr::set_names(c('Full.lm', 'regsubset', 'stepAIC')) %>% 
    bind_rows(., .id='Type') # %>%
    # pivot_wider(names_from = name)
  
  # note if df.coef.1, 2, and 3 have no observations, you cannont arrage by adjR2:
  
  if(dim(df.coef)[1]!=0){
    df.coef <- df.coef %>%
      arrange(desc(adjR2))
  }
  
  return(df.coef) 
  
}

# trouble shoot function:

# for (i in seq_along(l.setup.MLR)){
# 
#   fun.compare.linear.model.selections(l.setup.MLR[[i]])
# }

# function to compute best MLR model
# note* uses function from above:

i <- 5

# df <- l.setup.MLR[[i]]

fun.build.best.lm <- function(df, nbest = 1){
  
  # initalize df.coef using function:
  
  df.compare <- fun.compare.linear.model.selections(df)
  
  # get signficant predictors from the best model:
  
  preds <- df.compare$name[df.compare$Type==df.compare$Type[which.max(df.compare$adjR2)]]
  
  # select term and these predictors only, and also add in the highly correlated one left in:
  
  df1 <- select(df, c(term, preds))
  df2 <- select(df, c(term, preds, highCorr[leave]))
  
  # build model with these predictor sets:
  
  m.1 <- lm(term~.,data=df1)
  m.2 <- lm(term~.,data=df2)
  m.3 <- lm(term~.,data=df)
  
  # extract model results:
  
  summary(m.1)
  summary(m.2)
  summary(m.3)
  
}

##Function to calculate outliers
FindOutliers <- function(data) {
  lowerq = quantile(data)[2]
  upperq = quantile(data)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  # we identify extreme outliers
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  result <- which(data > extreme.threshold.upper | data < extreme.threshold.lower)
}

# function to plot one or multiple NWIS sites by site no:

fun.map.DA <- function(v.site_no) {
  
  load('Processed_Data/NWIS_Watershed_Shapefiles.62.Rdata')
  temp<-df.sf.NWIS.62%>%filter(Name %in% v.site_no)
  mapview(temp)
  
}

# function to normalize data between zero and 1:

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# function to make plot of pca biplot by observation:

fun.PCA.biplot.individuals <- function(pc.obj){
  
  p <- fviz_pca_ind(pc.obj,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  
  return(p)
  
}

# function to perform test/train CV for MLR:
# mxkuhn@gmail.com

trC.repcv <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
trC.cv <- trainControl(method = "cv", number = 10)
trC.none <- trainControl(method = "none")
trC.LGOCV <- trainControl(method = "LGOCV", p = p)
trC.boot <- trainControl(method = "boot")
trC.repcv <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
trC.LOOCV <- trainControl(method = "LOOCV")

# df <- l.resp.allpred[[1]]

fun.LGOCV.MLR <- function(df, p = 0.8){
  
  set.seed(2) # 2 give r2 of >0.7 ?!?!
  
  # set up training samples:
  
  training.samples <- df[,1] %>% createDataPartition(p = p, list = FALSE)
  
  # partion data:
  
  train.data  <- df[training.samples, ]
  test.data <- df[-training.samples, ]
  
  # build models:
  
  model <- train(term ~., data = train.data, method = "leapSeq", trControl = trC.LGOCV)
  
  # use model to predict test data:
  
  model.predictions <- predict(model, test.data)
  
  # get model performance by comparing test predictions and test observations:
  
  postResample(pred = model.predictions, obs = test.data$term)

  # look at model terms:
  
  coef(model$finalModel, as.numeric(model$bestTune))
  
  # create df of results:
  
  
  
  # method 2: use trainControl:
  
  # Set up LGOCV cross-validation:
  
  train.control <- trainControl(method = "LGOCV", p = p)
  
  # build model:
  
  model <- train(term ~., 
                      data = df,
                      method = "leapForward",
                      trControl = train.control
  )
  
  model.predictions <- predict(model, df)
  
  # get model performanceby comparing test predictions and test observations:
  
  postResample(pred = model.predictions, obs = df$term)
  
  # look at model terms:
  
  coef(model$finalModel, as.numeric(model$bestTune))
  
  
  
  
  
  
  
  # 
  
  ?resamples
  
  
  
  
}
  