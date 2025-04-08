library(tidyverse)
library(skimr)

d <- read_csv("./data/Black-backed jackal, Etosha National Park, Namibia.csv", col_names=TRUE)
head(d)
sum(is.na(d))
unique(d$visible)
d <- d %>%
  filter(year(timestamp) >= 2009 & year(timestamp) <= 2010) %>% # filter 2009 & 2010
  filter(`individual-local-identifier` %in% c("CM05","CM08","CM09","CM10","CM15",
                                              "CM44","CM47","CM62","CM69","CM83")) %>%
  dplyr::select(timestamp,`location-long`,`location-lat`,`individual-local-identifier`) %>% 
  arrange(`individual-local-identifier`, timestamp)
skim(d)

# GPS locations with duplicates removed
d <- d %>%
  distinct(timestamp,`location-long`,`location-lat`,`individual-local-identifier`, .keep_all = TRUE)

#########################
library(ggplot2)

ggplot(d, aes(x = `location-long`, y = `location-lat`, group = `individual-local-identifier`, color = `individual-local-identifier`)) +
  geom_path(linewidth = 0.8, alpha = 0.2) +
  geom_point(size = 0.1, alpha = 0.2) +
  coord_fixed() +
  labs(title = "Black-backed Jackal Trajectories",
       x = "Longitude", y = "Latitude", color = "Jackals") +
  theme_minimal()

#########################
# split data by id
d_list <- d %>% 
  group_split(`individual-local-identifier`) %>% 
  set_names(unique(d$`individual-local-identifier`))

# The temporal threshold used to define simultaneous fixes was 1,800 s.

# number of simultaneous fixes

count_simultaneous_fixes <- function(data1, data2, time_threshold = 1800) {
  # Ensure timestamps are in POSIXct format
  data1$timestamp <- as.POSIXct(data1$timestamp)
  data2$timestamp <- as.POSIXct(data2$timestamp)
  # Compute absolute time differences (in seconds) between all timestamp pairs
  time_diff_matrix <- outer(
    as.numeric(data1$timestamp),
    as.numeric(data2$timestamp),
    FUN = function(x, y) abs(x - y)
  )
  # Count the number of timestamp pairs within the time threshold
  simultaneous_fix_count <- sum(time_diff_matrix <= time_threshold)
  
  return(simultaneous_fix_count)
}


count_simultaneous_fixes(d_list$CM05, d_list$CM44)
count_simultaneous_fixes(d_list$CM08, d_list$CM09)
count_simultaneous_fixes(d_list$CM08, d_list$CM10)
count_simultaneous_fixes(d_list$CM09, d_list$CM10)
count_simultaneous_fixes(d_list$CM09, d_list$CM83)
count_simultaneous_fixes(d_list$CM15, d_list$CM47)
count_simultaneous_fixes(d_list$CM44, d_list$CM69)
count_simultaneous_fixes(d_list$CM44, d_list$CM83)
count_simultaneous_fixes(d_list$CM47, d_list$CM83)
count_simultaneous_fixes(d_list$CM62, d_list$CM83)

# mean proximity(m)

library(sf)
library(sp)

mean_proximity <- function(data1, data2, time_threshold = 1800) {
  # Ensure timestamps are in POSIXct format
  data1$timestamp <- as.POSIXct(data1$timestamp)
  data2$timestamp <- as.POSIXct(data2$timestamp)
  
  # Create time difference matrix (seconds)
  time_diff_matrix <- outer(
    as.numeric(data1$timestamp),
    as.numeric(data2$timestamp),
    FUN = function(x, y) abs(x - y)
  )
  
  # Create coordinate matrix (in meters)
  coords1 <- st_coordinates(st_as_sf(data1, coords = c("location-long", "location-lat"), crs = 4326))  # lon-lat
  coords2 <- st_coordinates(st_as_sf(data2, coords = c("location-long", "location-lat"), crs = 4326))
  
  # Distance matrix (great-circle distance, meters)
  distance_matrix <- spDists(coords1, coords2, longlat = TRUE) * 1000
  
  # Filter distances where time difference is within threshold
  valid_distances <- distance_matrix[time_diff_matrix <= time_threshold]
  
  if (length(valid_distances) == 0) {
    return(NA)  # no simultaneous fixes
  } else {
    return(mean(valid_distances))
  }
}

mean_proximity(d_list$CM05, d_list$CM44)
mean_proximity(d_list$CM08, d_list$CM09)
mean_proximity(d_list$CM08, d_list$CM10)
mean_proximity(d_list$CM09, d_list$CM10)
mean_proximity(d_list$CM09, d_list$CM83)
mean_proximity(d_list$CM15, d_list$CM47)
mean_proximity(d_list$CM44, d_list$CM69)
mean_proximity(d_list$CM44, d_list$CM83)
mean_proximity(d_list$CM47, d_list$CM83)
mean_proximity(d_list$CM62, d_list$CM83)


###############
library(dplyr)
library(geosphere)

calculate_dyad_summary <- function(data1, data2, name1, name2, time_threshold = 1800) {
  # POSIXct 변환
  data1$timestamp <- as.POSIXct(data1$timestamp)
  data2$timestamp <- as.POSIXct(data2$timestamp)
  
  # 두 데이터의 모든 조합 만들기
  df1 <- data1[rep(1:nrow(data1), each = nrow(data2)), ]
  df2 <- data2[rep(1:nrow(data2), times = nrow(data1)), ]
  
  # 시간 차 계산
  time_diff <- abs(as.numeric(df1$timestamp) - as.numeric(df2$timestamp))
  valid <- time_diff <= time_threshold
  
  if (sum(valid) == 0) {
    return(data.frame(
      Dyad = paste(name1, "and", name2),
      Simultaneous_Fixes = 0,
      Mean_Proximity = NA,
      Contacts_100 = 0,
      Contacts_500 = 0,
      Contacts_1000 = 0,
      Contacts_5000_5500 = 0
    ))
  }
  
  # 유효한 쌍 좌표 추출
  coords1 <- df1[valid, c("location-long", "location-lat")]
  coords2 <- df2[valid, c("location-long", "location-lat")]
  distances <- distHaversine(coords1, coords2)
  
  # 거리 기반 contact 개수
  contacts_100 <- sum(distances <= 100)
  contacts_500 <- sum(distances <= 500)
  contacts_1000 <- sum(distances <= 1000)
  contacts_5500 <- sum(distances >= 5000 & distances <= 5500)
  
  data.frame(
    Dyad = paste(name1, "and", name2),
    Simultaneous_Fixes = length(distances),
    Mean_Proximity = round(mean(distances), 1),
    Contacts_100 = contacts_100,
    Contacts_500 = contacts_500,
    Contacts_1000 = contacts_1000,
    Contacts_5000_5500 = contacts_5500
  )
}


dyad_list <- list(
  c("CM05", "CM44"),
  c("CM08", "CM09"),
  c("CM08", "CM10"),
  c("CM09", "CM10"),
  c("CM09", "CM83"),
  c("CM15", "CM47"),
  c("CM44", "CM69"),
  c("CM44", "CM83"),
  c("CM47", "CM83"),
  c("CM62", "CM83")
)

results <- lapply(dyad_list, function(pair) {
  calculate_dyad_summary(d_list[[pair[1]]], d_list[[pair[2]]], pair[1], pair[2])
})

final_table <- do.call(rbind, results)

final_table




















##

library(move2)
CM09 <- st_as_sf(d_list$CM09, coords = c("location-long", "location-lat"), crs = 4326)
CM83 <- st_as_sf(d_list$CM83, coords = c("location-long", "location-lat"), crs = 4326)
CM09 <- mt_as_move2(CM09, time_column = "timestamp", track_id_column = "individual-local-identifier")
CM83 <- mt_as_move2(CM83, time_column = "timestamp", track_id_column = "individual-local-identifier")

#####
BoDa<-mt_stack(CM09, CM83)
dyad<-BoDa
z<-conProcess(dyad, dc=10**10, tc=1800)
t<-as.data.frame(z)

b<-subset(t, contact_n==1) 

