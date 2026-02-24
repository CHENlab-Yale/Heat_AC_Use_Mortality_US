rm(list=ls())

library(stats); library(splines); library(gam); library(dlnm); library(gnm); library(ggplot2)
library(dplyr)

### 0. basic settings
set.seed(0219)
# read the Chicago dataset from dlnm package
Chicago <- chicagoNMMAPS
# save path
save.path <- "~/"


### 1. format Chicago into individual-level long table
str(Chicago)
## (1) exposure data
exposure.data <- Chicago %>% 
  dplyr::select(date, temp, dptp, pm10, o3) %>% 
  rename(tavg = temp, dp=dptp)
saveRDS(exposure.data, file = paste0(save.path, "Data/Data.0.Chicago_exposure_ts.rds"))

## (2) health data
health.data <- Chicago %>% 
  dplyr::select(date, month, cvd) %>% 
  rename(case.month = month)
sum(health.data$cvd) # 260314
health.data.long <- health.data[rep(1:nrow(health.data), times=health.data$cvd),]
dim(health.data.long) # 260314 ===> correct
health.data.long <- health.data.long %>% 
  dplyr::select(-cvd) %>% 
  mutate(case.status = 1) %>% 
  rename(case.date = date) %>% 
  mutate(ID = 1:nrow(health.data.long)) %>% 
  relocate(ID, .before="case.date")
rownames(health.data.long) <- NULL


### 2. select controls for the cases
health.raw <- health.data.long
health.data.all <- health.data.long %>% mutate(date = case.date, 
                                               date.month = lubridate::month(date))
for (lag.week in 1:4) {
  control1 <- health.raw %>% 
    mutate(case.status = 0) %>% 
    mutate(date = case.date - lag.week * 7,
           date.month = lubridate::month(date)) %>% 
    filter(date.month == case.month)
  control2 <- health.raw %>% 
    mutate(case.status = 0) %>% 
    mutate(date = case.date + lag.week * 7,
           date.month = lubridate::month(date)) %>% 
    filter(date.month == case.month)
  health.data.all <- rbind(health.data.all, control1, control2)
  
}
dim(health.data.all) # 1145267
nrow(health.data.all)/nrow(health.data.long) # 4.40 ===> correct
## check
sum(health.data.all$case.month==health.data.all$date.month) # 1145267 ===> correct!

### 3. link with the exposure data
str(health.data.all)
# 'data.frame':	1145267 obs. of  6 variables:
# $ ID         : int  1 2 3 4 5 6 7 8 9 10 ...
# $ case.date  : Date, format: "1987-01-01" "1987-01-01" ...
# $ case.month : num  1 1 1 1 1 1 1 1 1 1 ...
# $ case.status: num  1 1 1 1 1 1 1 1 1 1 ...
# $ date       : Date, format: "1987-01-01" "1987-01-01" ...
# $ date.month : num  1 1 1 1 1 1 1 1 1 1 ...
str(exposure.data)
# 'data.frame':	5114 obs. of  5 variables:
# $ date: Date, format: "1987-01-01" "1987-01-02" ...
# $ tavg: num  -0.278 0.556 0.556 -1.667 0 ...
# $ dp  : num  31.5 29.9 27.4 28.6 28.9 ...
# $ pm10: num  27 NA 32.8 40 NA ...
# $ o3  : num  4.38 4.93 3.75 4.29 4.75 ...
## (1) create lagged exposures
data.all <- health.data.all
for (lag.day in 0:6) {
  data.all$date.link <- data.all$date - lag.day
  exposure.data.single.lag <- exposure.data %>% 
    relocate(date, .before = 1)
  colnames(exposure.data.single.lag) <- paste0(colnames(exposure.data.single.lag), "_lag", lag.day)
  colnames(exposure.data.single.lag)[1] <- "date.link"
  data.all <- data.all %>% 
    left_join(exposure.data.single.lag, by = "date.link")
  
}
str(data.all)
saveRDS(data.all, file = paste0(save.path, "Data/Data.1.Chicago_cc_format_withlag.rds"))
data.sub <- data.all %>% 
  dplyr::select(-case.date, -case.month, -date.month)
saveRDS(data.sub, file = paste0(save.path, "Data/Data.2.Chicago_cc_format_withlag_clean.rds"))
## create fake county
fake.county <- data.frame(ID=unique(data.sub$ID))
fake.county$fake_county <- sample(c("a","b","c","d","e","f","g"), size = nrow(fake.county), replace = TRUE)
table(fake.county$fake_county)
# a     b     c     d     e     f     g 
# 37384 37214 37059 37225 37140 37096 37196 
data.sub.withcounty <- data.sub %>% 
  left_join(fake.county, by = "ID") %>% 
  relocate(fake_county, .before = 1)
table(data.sub.withcounty$fake_county)
# a      b      c      d      e      f      g 
# 164472 163763 162936 163743 163585 163212 163556 
saveRDS(data.sub.withcounty, file = paste0(save.path, "Data/Data.3.Chicago_cc_format_withlag_clean_with_fake_county.rds"))


### 4. exposure distribution information
# PM10 had 251 NA out of 5114 observations
for (expi in c("tavg", "dp", "pm10", "o3")) {
  single.data <- exposure.data %>% 
    dplyr::select(all_of(expi))
  assign(paste0(expi,"_dist"),quantile(single.data,seq(0,1,0.001),na.rm=TRUE))
}
pct <- seq(0,1,0.001)*100
dist.results <- cbind(pct,do.call(cbind,lapply(paste0(c("tavg", "dp", "pm10", "o3"),"_dist"),get)))
colnames(dist.results) <- c("pct", "tavg", "dp", "pm10", "o3")
saveRDS(dist.results, file = paste0(save.path, "Results/Results.0.Dist_all.rds"))
