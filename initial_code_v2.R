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

#############################################################################################################################
library(ggplot2)

ggplot(d, aes(x = `location-long`, y = `location-lat`, group = `individual-local-identifier`, color = `individual-local-identifier`)) +
  geom_path(linewidth = 0.8, alpha = 0.2) +
  geom_point(size = 0.1, alpha = 0.2) +
  coord_fixed() +
  labs(title = "Black-backed Jackal Trajectories",
       x = "Longitude", y = "Latitude", color = "Jackals") +
  theme_minimal()

# split the data by individual id
d_list <- d %>% 
  group_split(`individual-local-identifier`) %>% 
  set_names(unique(d$`individual-local-identifier`))

#############################################################################################################################
# Description of jackal dyad
# table 사진 첨부
# The temporal threshold used to define simultaneous fixes was 1,800 s.
# number of simultaneous fixes, mean proximity (m), and the number of contacts

# Make function for table 2.
table_2_function <- function(data1, data2, time_threshold = 1800) {
  library(sf)
  library(sp)
  
  ########## number of simultaneous fixes ##########
  # Compute absolute time differences (in seconds) between all timestamp pairs using matrix
  time_diff_matrix <- outer(
    as.numeric(data1$timestamp),
    as.numeric(data2$timestamp),
    FUN = function(x, y) abs(x - y)
  )
  # Count the number of timestamp pairs within the time threshold
  simultaneous_fix_count <- sum(time_diff_matrix <= time_threshold)
  
  ########## mean proximity (m) ##########
  # Create coordinate matrix
  coords1 <- st_coordinates(st_as_sf(data1, coords = c("location-long", "location-lat"), crs = 4326))
  coords2 <- st_coordinates(st_as_sf(data2, coords = c("location-long", "location-lat"), crs = 4326))
  # Distance matrix (meters)
  distance_matrix <- spDists(coords1, coords2, longlat = TRUE) * 1000
  # Filter distances where time difference is within threshold
  simultaneous_distances <- distance_matrix[time_diff_matrix <= time_threshold]
  
  ########## number of contacts at specified distance ##########
  # Count contacts at distance thresholds
  contact_100m <- sum(simultaneous_distances <= 100)
  contact_500m <- sum(simultaneous_distances <= 500)
  contact_1000m <- sum(simultaneous_distances <= 1000)
  contact_5000_5500m <- sum(simultaneous_distances >= 5000 & simultaneous_distances <= 5500)
  
  ########## output ##########
  return(data.frame(
    simultaneous_fix_count = simultaneous_fix_count,
    mean_proximity_m = round(mean(simultaneous_distances)),
    contact_100m = contact_100m,
    contact_500m = contact_500m,
    contact_1000m = contact_1000m,
    contact_5000_5500m = contact_5000_5500m
  ))
}

# results

# Dyads
dyads <- list(
  c("CM05", "CM44"), c("CM08", "CM09"), c("CM08", "CM10"), c("CM09", "CM10"), c("CM09", "CM83"),
  c("CM15", "CM47"), c("CM44", "CM69"), c("CM44", "CM83"), c("CM47", "CM83"), c("CM62", "CM83")
)
# use lapply 
results <- lapply(dyads, function(pair) {
  table_2_function(d_list[[pair[1]]], d_list[[pair[2]]])
})
# dataframe
df <- do.call(rbind, lapply(seq_along(results), function(i) {
  cbind(Dyad = paste(dyads[[i]][1], "and", dyads[[i]][2]), results[[i]])
}))
# kable
library(kableExtra)
df %>% 
  kable(align = "c") %>% 
  kable_classic_2(full_width = FALSE, html_font = "Cambria")

#############################################################################################################################
# Calculate Movement Parameters (3.2)
# table 사진 첨부

# individual's movement parameters
library(adehabitatLT)
# 05 예시
traj <- as.ltraj(
  xy = as.data.frame(d_list$CM05[, c("location-long", "location-lat")]),
  date = d_list$CM05$timestamp,
  id = d_list$CM05$`individual-local-identifier`
)

move_params <- ld(traj)

move_params <- move_params %>%
  mutate(
    speed = dist / as.numeric(dt),         # Speed (m/s)
    PI = cos(rel.angle),                   # Persistence Index (PI)
    Vp = speed * PI,                       # Persistence Velocity (Vp)
    Vt = speed * sin(rel.angle)            # Turning Velocity (Vt)
  ) %>%
  dplyr::select(x, y, date, id, dist, speed, abs.angle, rel.angle, PI, Vp, Vt)


# similarity of movement of two individuals in a dyad
library(wildlifeDI)
library(sf)
library(move2)
# 05 44 예시
matches <- expand.grid(1:nrow(d_list$CM05), 1:nrow(d_list$CM44))
matches <- matches[abs(as.numeric(difftime(d_list$CM05$timestamp[matches$Var1],
                                           d_list$CM44$timestamp[matches$Var2],
                                           units = "secs"))) <= 1800, ]
sim_data1 <- d_list$CM05[matches$Var1, ]
sim_data2 <- d_list$CM44[matches$Var2, ] # 있으나 마나 인듯
# 1. sf 변환
sim_data1_sf <- st_as_sf(sim_data1, coords = c("location-long", "location-lat"), crs = 4326)
sim_data2_sf <- st_as_sf(sim_data2, coords = c("location-long", "location-lat"), crs = 4326)
# 2. move2 변환
m1 <- mt_as_move2(sim_data1_sf, time_column = "timestamp", track_id_column = "individual-local-identifier")
m2 <- mt_as_move2(sim_data2_sf, time_column = "timestamp", track_id_column = "individual-local-identifier")
# 3. 병합 및 정렬
m <- mt_stack(m1, m2)
# 4. DI 계산
DI(dyad, tc=1800)
