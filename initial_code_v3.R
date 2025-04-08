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
  distinct(timestamp,`individual-local-identifier`, .keep_all = TRUE)

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
  distance_matrix <- spDists(coords1, coords2, longlat = TRUE) * 1000 # spDists function works with WGS84 ellipsoid
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
move_params_function <- function(data_1){
  # 먼저 wgs84 좌표계를 UTM 좌표계로 투영시켜줘야 미터 계산이 가능하다. 따라서 연구 대상지와 맞는 좌표계 선택하여 투영한다.
  library(sf)
  data_sf <- st_as_sf(data_1, coords = c("location-long", "location-lat"), crs = 4326)
  data_proj <- st_transform(data_sf, crs = 32733)  # UTM 33S (남반구)
  # 좌표 추출
  coords_proj <- st_coordinates(data_proj)
  # 원래 데이터에 x, y 열 붙이기
  data_1$x <- coords_proj[,1]
  data_1$y <- coords_proj[,2]
  
  # movement parameters를 계산한다.
  library(adehabitatLT)
  traj <- as.ltraj(
    xy = as.data.frame(data_1[, c("x", "y")]),
    date = data_1$timestamp,
    id = data_1$`individual-local-identifier`
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
  
  return(move_params)
}

CM05_move_params <- move_params_function(d_list$CM05)
CM08_move_params <- move_params_function(d_list$CM08)
CM09_move_params <- move_params_function(d_list$CM09)
CM10_move_params <- move_params_function(d_list$CM10)
CM15_move_params <- move_params_function(d_list$CM15)
CM44_move_params <- move_params_function(d_list$CM44)
CM47_move_params <- move_params_function(d_list$CM47)
CM62_move_params <- move_params_function(d_list$CM62)
CM69_move_params <- move_params_function(d_list$CM69)
CM83_move_params <- move_params_function(d_list$CM83)

#############################################################################################################################
# contact 열 만들기

contact_flag_function <- function(data1, data2, time_threshold=1800, distance_threshold) {
  library(dplyr)
  library(sf)
  library(sp)
  # 시간 차 행렬
  time_diff_matrix <- outer(data1$date, data2$date, FUN = function(x, y) abs(x - y))
  # 거리 행렬 계산 (단위: meter)
  coords_1 <- as.matrix(data1[, c("x", "y")]) # 이미 meter 단위 투영 좌표계라 그대로 사용
  coords_2 <- as.matrix(data2[, c("x", "y")])
  distance_matrix <- spDists(coords_1, coords_2, longlat = FALSE) # 이미 meter 단위 투영 좌표계라 longlat = FALSE
  
  # 유효한 interaction 인덱스 (TRUE인 위치)
  valid_interactions <- which((time_diff_matrix <= time_threshold) & (distance_matrix <= distance_threshold), arr.ind = TRUE)
  
  # 초기화
  data1$contact <- FALSE
  data2$contact <- FALSE
  
  # contact 표시 + 앞뒤 행도 contact로 만들기
  for (k in 1:nrow(valid_interactions)) {
    i <- valid_interactions[k, 1]
    j <- valid_interactions[k, 2]
    
    data1$contact[i] <- TRUE
    data2$contact[j] <- TRUE
    
    # 앞뒤 행 처리 (조건 체크 포함)
    if (i > 1) {data1$contact[i - 1] <- TRUE}
    if (i < nrow(data1)) {data1$contact[i + 1] <- TRUE}
    if (j > 1) {data2$contact[j - 1] <- TRUE}
    if (j < nrow(data2)) {data2$contact[j + 1] <- TRUE}
  }

  return(list(data1, data2))
}

jackal_contact_list <- contact_flag_function(CM05_move_params, CM44_move_params,distance_threshold = 100)
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM08_move_params, CM09_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM08_move_params, CM10_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM09_move_params, CM10_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM09_move_params, CM83_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM15_move_params, CM47_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM44_move_params, CM69_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM44_move_params, CM83_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM47_move_params, CM83_move_params,distance_threshold = 100))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM62_move_params, CM83_move_params,distance_threshold = 100))

#############################################################################################################################
# bootstrap

bootstrap_function <- function(data_jackle){
  library(circular)
  library(stats)
  set.seed(123)
  
  # 그룹 나누기
  contact_group <- data_jackle %>% filter(contact == TRUE)
  other_group <- data_jackle %>% filter(contact == FALSE)
  # 비교할 변수
  vars_watsons <- c("abs.angle", "rel.angle")  # Watson
  vars_ks <- c("speed", "dist", "PI", "Vp", "Vt")  # K-S
  # 부트스트랩 반복 수
  k <- 1000
  # 결과 저장 리스트
  bootstrap_stats <- list()
  bootstrap_pvals <- list()
  # 공통 sample size
  sample_size <- min(nrow(contact_group), nrow(other_group))
  
  # --- 1. Watson's U² 통계량 부트스트랩 ---
  for (var in vars_watsons) {
    stat_dist <- replicate(k, {
      x <- sample(na.omit(contact_group[[var]]), sample_size, replace = TRUE)
      y <- sample(na.omit(other_group[[var]]), sample_size, replace = TRUE)
      watson.two.test(as.circular(x, units = "radians"), as.circular(y, units = "radians"))$statistic
    })
    # 실제 통계량
    x_real <- na.omit(contact_group[[var]])
    y_real <- na.omit(other_group[[var]])
    real_stat <- watson.two.test(as.circular(x_real, units = "radians"), as.circular(y_real, units = "radians"))$statistic
    # p-value: 부트스트랩된 통계량이 실제보다 큰 비율
    p_val <- mean(stat_dist >= real_stat)
    
    bootstrap_stats[[var]] <- stat_dist
    bootstrap_pvals[[var]] <- p_val
  }
  
  # --- 2. K–S 통계량 부트스트랩 ---
  for (var in vars_ks) {
    stat_dist <- replicate(k, {
      x <- sample(na.omit(contact_group[[var]]), sample_size, replace = TRUE)
      y <- sample(na.omit(other_group[[var]]), sample_size, replace = TRUE)
      ks.test(x, y)$statistic
    })
    # 실제 통계량
    x_real <- na.omit(contact_group[[var]])
    y_real <- na.omit(other_group[[var]])
    real_stat <- ks.test(x_real, y_real)$statistic
    # p-value 계산
    p_val <- mean(stat_dist >= real_stat)
    
    bootstrap_stats[[var]] <- stat_dist
    bootstrap_pvals[[var]] <- p_val
  }
  
  return(bootstrap_pvals)
}

# ------result (100m)-------
# apply to 20 trajectories
bootstrap_pvals_list <- lapply(1:20, function(i) {
  bootstrap_function(jackal_contact_list[[i]])
})
# dataframe
bootstrap_pval_df <- do.call(rbind, lapply(seq_along(bootstrap_pvals_list), function(i) {
  data.frame(
    case = paste0("case_", i),
    variable = names(bootstrap_pvals_list[[i]]),
    p_value = unlist(bootstrap_pvals_list[[i]])
  )
}))
# wide format으로 변환
bootstrap_pval_df <- bootstrap_pval_df %>%
  pivot_wider(names_from = variable, values_from = p_value)
# 열별로 유의한 비율 계산
proportion_per_variable <- colMeans(bootstrap_pval_df[,-1] < 0.05, na.rm = TRUE)
# 새 행으로 추가
bootstrap_pval_df_100m <- rbind(
  bootstrap_pval_df,
  c(case = "Proportion_significant", proportion_per_variable)
)

#############################################################################################################################
# ------result (500m)-------
jackal_contact_list <- contact_flag_function(CM05_move_params, CM44_move_params,distance_threshold = 500)
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM08_move_params, CM09_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM08_move_params, CM10_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM09_move_params, CM10_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM09_move_params, CM83_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM15_move_params, CM47_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM44_move_params, CM69_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM44_move_params, CM83_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM47_move_params, CM83_move_params,distance_threshold = 500))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM62_move_params, CM83_move_params,distance_threshold = 500))
# apply to 20 trajectories
bootstrap_pvals_list <- lapply(1:20, function(i) {
  bootstrap_function(jackal_contact_list[[i]])
})
# dataframe
bootstrap_pval_df <- do.call(rbind, lapply(seq_along(bootstrap_pvals_list), function(i) {
  data.frame(
    case = paste0("case_", i),
    variable = names(bootstrap_pvals_list[[i]]),
    p_value = unlist(bootstrap_pvals_list[[i]])
  )
}))
# wide format으로 변환
bootstrap_pval_df <- bootstrap_pval_df %>%
  pivot_wider(names_from = variable, values_from = p_value)
# 열별로 유의한 비율 계산
proportion_per_variable <- colMeans(bootstrap_pval_df[,-1] < 0.05, na.rm = TRUE)
# 새 행으로 추가
bootstrap_pval_df_500m <- rbind(
  bootstrap_pval_df,
  c(case = "Proportion_significant", proportion_per_variable)
)

#############################################################################################################################
# ------result (1000m)-------
jackal_contact_list <- contact_flag_function(CM05_move_params, CM44_move_params,distance_threshold = 1000)
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM08_move_params, CM09_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM08_move_params, CM10_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM09_move_params, CM10_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM09_move_params, CM83_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM15_move_params, CM47_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM44_move_params, CM69_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM44_move_params, CM83_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM47_move_params, CM83_move_params,distance_threshold = 1000))
jackal_contact_list <- append(jackal_contact_list, contact_flag_function(CM62_move_params, CM83_move_params,distance_threshold = 1000))
# apply to 20 trajectories
bootstrap_pvals_list <- lapply(1:20, function(i) {
  bootstrap_function(jackal_contact_list[[i]])
})
# dataframe
bootstrap_pval_df <- do.call(rbind, lapply(seq_along(bootstrap_pvals_list), function(i) {
  data.frame(
    case = paste0("case_", i),
    variable = names(bootstrap_pvals_list[[i]]),
    p_value = unlist(bootstrap_pvals_list[[i]])
  )
}))
# wide format으로 변환
bootstrap_pval_df <- bootstrap_pval_df %>%
  pivot_wider(names_from = variable, values_from = p_value)
# 열별로 유의한 비율 계산
proportion_per_variable <- colMeans(bootstrap_pval_df[,-1] < 0.05, na.rm = TRUE)
# 새 행으로 추가
bootstrap_pval_df_1000m <- rbind(
  bootstrap_pval_df,
  c(case = "Proportion_significant", proportion_per_variable)
)

#############################################################################################################################
#contact & other 구분용 함수 (5000–5500m, 5500m 초과)
contact_distance_flag_function <- function(data1, data2, time_threshold = 1800, min_distance = 5000, max_distance = 5500) {
  library(dplyr)
  library(sp)
  # 시간 차 행렬 (초 단위 차이)
  time_diff_matrix <- outer(data1$date, data2$date, FUN = function(x, y) abs(x - y))
  # 거리 행렬 (이미 meter 단위 투영 좌표계라고 가정)
  coords_1 <- as.matrix(data1[, c("x", "y")])
  coords_2 <- as.matrix(data2[, c("x", "y")])
  distance_matrix <- spDists(coords_1, coords_2, longlat = FALSE)
  # 조건 1: contact (5000m 이상 5500m 이하 + 시간 조건)
  contact_condition <- (distance_matrix >= min_distance & distance_matrix <= max_distance) & (time_diff_matrix <= time_threshold)
  # 조건 2: other (5500m 초과 + 시간 조건)
  other_condition <- (distance_matrix > max_distance)
  
  # contact, other 초기화
  data1$contact <- NA
  data2$contact <- NA
  # --- contact 플래그 설정 ---
  valid_contacts <- which(contact_condition, arr.ind = TRUE)
  for (k in 1:nrow(valid_contacts)) {
    i <- valid_contacts[k, 1]
    j <- valid_contacts[k, 2]
    data1$contact[i] <- TRUE
    data2$contact[j] <- TRUE
    # 앞뒤 한 행씩도 contact로
    if (i > 1) data1$contact[i - 1] <- TRUE
    if (i < nrow(data1)) data1$contact[i + 1] <- TRUE
    if (j > 1) data2$contact[j - 1] <- TRUE
    if (j < nrow(data2)) data2$contact[j + 1] <- TRUE
  }
  # --- other 플래그 설정 ---
  other_rows_1 <- which(is.na(data1$contact) & apply(other_condition, 1, any))
  other_rows_2 <- which(is.na(data2$contact) & apply(other_condition, 2, any))
  
  data1$contact[other_rows_1] <- FALSE
  data2$contact[other_rows_2] <- FALSE
  
  return(list(data1, data2))
}

# ------result (5000 ~ 5500m)-------
jackal_contact_list <- contact_distance_flag_function(CM05_move_params, CM44_move_params)
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM08_move_params, CM09_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM08_move_params, CM10_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM09_move_params, CM10_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM09_move_params, CM83_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM15_move_params, CM47_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM44_move_params, CM69_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM44_move_params, CM83_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM47_move_params, CM83_move_params))
jackal_contact_list <- append(jackal_contact_list, contact_distance_flag_function(CM62_move_params, CM83_move_params))
# apply to 20 trajectories
bootstrap_pvals_list <- lapply(1:20, function(i) {
  bootstrap_function(jackal_contact_list[[i]])
})
# dataframe
bootstrap_pval_df <- do.call(rbind, lapply(seq_along(bootstrap_pvals_list), function(i) {
  data.frame(
    case = paste0("case_", i),
    variable = names(bootstrap_pvals_list[[i]]),
    p_value = unlist(bootstrap_pvals_list[[i]])
  )
}))
# wide format으로 변환
bootstrap_pval_df <- bootstrap_pval_df %>%
  pivot_wider(names_from = variable, values_from = p_value)
# 열별로 유의한 비율 계산
proportion_per_variable <- colMeans(bootstrap_pval_df[,-1] < 0.05, na.rm = TRUE)
# 새 행으로 추가
bootstrap_pval_df_5000_5500m <- rbind(
  bootstrap_pval_df,
  c(case = "Proportion_significant", proportion_per_variable)
)
