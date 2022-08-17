
######################################

library(ggplot2)
library(ggmap)
library(RgoogleMaps)
library(viridis)
library(ggsci)
library(sp)
library(rgdal)


crs_temp32 <- CRS("+proj=utm +zone=32 +ellps=WGS84")
crs_temp33 <- CRS("+proj=utm +zone=33 +ellps=WGS84")
crs_temp33_var <- CRS("+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs")
crs_temp_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

### Reading shapefile with municipality info; 

Map1 <- readOGR("R://GeoSpatialData/AdministrativeUnits/Norway_AdministrativeUnits/Processed/Kommuner_2018_WGS84.shp", 
                stringsAsFactors=FALSE, encoding="UTF-8", use_iconv=TRUE)
Map1 <- sp::spTransform(Map1, crs_temp_longlat)

Map2 <- readOGR("R://GeoSpatialData/AdministrativeUnits/Norway_AdministrativeUnits/Processed/Fylker_2019_ERTS89_33N.shp", 
                stringsAsFactors=FALSE, encoding="UTF-8", use_iconv=TRUE)
Map2 <- sp::spTransform(Map2, crs_temp_longlat)

## Fig: Norway; 

p2 <- ggplot() + 
  xlim(-20, 50)+
  ylim(40, 80) +
  borders("world", fill="grey90", colour = "chartreuse4") +
  borders(regions = "Norway(?!:Svalbard)", colour = "chartreuse4", fill = "chartreuse4") +
  theme_minimal()

p2


ggsave("figurer/Map_europe.png", p2, bg = "transparent")



######################################################################
# FIG 1: MAP of detections from line transect data set;

library(RODBC)
myconn <-odbcDriverConnect("driver={SQL Server}; server=ninsql07; database=Honsefugl; trusted_connection=true")

Obs_map <- as_tibble(sqlFetch(myconn, "Observasjon"))
close(myconn)

obs_map1 <- Obs_map %>% filter(Aar==2019 & LinjeAvstand<200) %>% select("ObservasjonId", "LinjeAvstand", "Latitude", "Longitude")



p3 <- ggplot() + 
  xlim(3, 33)+
  ylim(57, 72) +
  borders(regions = "Norway(?!:Svalbard)", colour = "gray50", fill = "chartreuse4") +
  scale_color_gradient2(midpoint = 0.5, low = "grey60", mid = "gold1",
                        high = "firebrick1", space = "Lab" )+
  geom_point(data=obs_map1, aes(x=Longitude, y=Latitude), col="orange", size=0.4) + 
  
  theme_minimal()

p3

ggsave("figurer/Map_line_transect.png", p3, bg = "transparent")


################################

library(rgbif)

d1 <- occ_search(scientificName = "Lagopus lagopus", hasCoordinate = TRUE, country = "NO", limit=10000)

d2 <- d1[[3]] %>% filter(year>1990) %>% select(key, decimalLongitude, decimalLatitude)

p4 <- ggplot() + 
  xlim(3, 33)+
  ylim(57, 72) +
  borders(regions = "Norway(?!:Svalbard)", colour = "gray50", fill = "chartreuse4") +
  scale_color_gradient2(midpoint = 0.5, low = "grey60", mid = "gold1",
                        high = "firebrick1", space = "Lab" )+
  geom_point(data=d2, aes(x=decimalLongitude, y=decimalLatitude), col="orange", size=1.2) + 
  
  theme_minimal()

p4

ggsave("figurer/Map_gbif.png", p4, bg = "transparent")


####################################################################

d4 <- data.frame(long=14.0, lat=64.4)

p4 <- ggplot() + 
  xlim(3, 33)+
  ylim(57, 72) +
  borders(regions = "Norway(?!:Svalbard)", colour = "gray50", fill = "chartreuse4") +
  scale_color_gradient2(midpoint = 0.5, low = "grey60", mid = "gold1",
                        high = "firebrick1", space = "Lab" )+
  geom_point(data=d4, aes(x=long, y=lat), col="orange", size=5) + 
  
  theme_minimal()
p4

ggsave("figurer/Map_Lierne.png", p4, bg = "transparent")


