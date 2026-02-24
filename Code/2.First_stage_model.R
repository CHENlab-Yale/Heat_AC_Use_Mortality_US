rm(list = ls())
gc()

### 0. load the libraries
library(dplyr)
library(dlnm)
library(splines2)
library(survival)


### 1. load the data
data.path <- save.path <- "~/"
raw.data <- readRDS(paste0(data.path, "Data/Data.3.Chicago_cc_format_withlag_clean_with_fake_county.rds"))
temp.var <- "tavg"
humid.var <- "dp" # dew point

### 2. exposure distribution
exp.dist <- readRDS(paste0(data.path, "Results/Results.0.Dist_all.rds"))
exp.dist <- as.data.frame(exp.dist)


### 3. basic settings
## Notes
# county-specific distributions should be used (moved into the coming for loop)
# but in this example, we only had one distribution in the Chicago data
temp.var.knots <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(10,75,90)),temp.var])
temp.boundary <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(0, 100)),temp.var])
humid.var.knots <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(25,50,75)),humid.var])
humid.boundary <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(0, 100)),humid.var])
formula.main<-"case.status ~ cbtemp + cbhumid + strata(ID)"
max.lag <- 6


### 4. county-specific model
county.seq <- sort(unique(raw.data$fake_county))
coef.all <- NULL
var.all <- list()
time1 <- Sys.time()
for (county.i in county.seq) {
  cat("County: ", county.i, "\n")
  data.single.county <- raw.data %>% 
    filter(fake_county == county.i)
  ## (1) extract the lagged exposure data and set up the splines
  temp.lag <- data.single.county %>% dplyr::select(all_of(paste0(temp.var,"_lag",0:max.lag)))
  humid.lag <- data.single.county %>% dplyr::select(all_of(paste0(humid.var,"_lag",0:max.lag)))
  cbtemp <- dlnm::crossbasis(temp.lag, lag = max.lag, argvar = list(fun = "ns", knots = temp.var.knots,
                                                                    Boundary.knots = temp.boundary),
                             arglag = list(fun = "ns", knots = 3, df = 3))
  cbhumid <- dlnm::crossbasis(humid.lag, lag = max.lag, argvar = list(fun = "ns", knots = humid.var.knots, 
                                                                      Boundary.knots = humid.boundary), 
                              arglag = list(fun = "ns", knots = 3, df = 3))
  ## (2) run the model
  main.model <- clogit(as.formula(formula.main), data=data.single.county, method = "breslow")
  ## (3) output
  results.coefficients <- as.data.frame(t(as.data.frame(main.model$coefficients)))
  results.var <- as.data.frame(main.model$var)
  colnames(results.var) <- colnames(results.coefficients)
  results.var$variable <- colnames(results.coefficients)
  rownames(results.coefficients) <- NULL
  results.coefficients$county <- county.i
  
  ## (4) combine across counties
  coef.all <- rbind(coef.all, results.coefficients)
  var.all[[paste0("county_",county.i)]] <- results.var
}
cat(print(Sys.time() - time1), "\n") # 20 sec
saveRDS(coef.all, file = paste0(save.path, "Results/Results.1st.coef_all.rds"))
saveRDS(var.all, file = paste0(save.path, "Results/Results.1st.var_all.rds"))

