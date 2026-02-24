rm(list = ls())
gc()

library(dplyr)
library(mvmeta)
library(openxlsx)
library(dlnm)
library(splines2)


### 0. basic settings
set.seed(0219)
data.path <- save.path <- "~/"



### 1. read in the association estimates
est.data <- readRDS(paste0(save.path, "Results/Results.1st.coef_all.rds"))
var.data <- readRDS(paste0(save.path, "Results/Results.1st.var_all.rds"))
county.seq <- est.data$county


### 2. clean the data
## (1) point estimate
est.matrix <- est.data
est.matrix <- est.matrix %>% 
  select(county, contains("cbtemp"))
## (2) variance
var.list <- list()
for (county.i in est.matrix$county) {
  single.var <- var.data[[paste0("county_",county.i)]]
  single.var <- single.var %>% 
    filter(grepl("cbtemp",variable)) %>% 
    select(contains("cbtemp"))
    
  var.list[[county.i]] <- single.var
}
est.matrix$county <- NULL
est.matrix <- as.matrix(est.matrix)


### 3. meta-analysis
## (1) create fake AC data
ac.fake <- rnorm(n=nrow(est.matrix), mean = 1300, sd=200)
saveRDS(ac.fake,file = paste0(save.path, "Data/Data.Z.Fake_AC.rds"))
cb.ac <- dlnm::onebasis(ac.fake, fun = "ns", knots = median(ac.fake), Boundary.knots = range(ac.fake))
cb.ac.df <- as.data.frame(cb.ac)
# due to the limited number of cases, only one knot was selected for AC
## (2) meta regression
meta.model <- mvmeta(est.matrix~b1+b2, var.list, data=cb.ac.df,
                     method = "reml",
                     control=list(maxiter=500, showiter=TRUE))
meta.coef <- meta.model$coefficients
meta.vcov <- meta.model$vcov
model.summary <- summary(meta.model)
heter.sum <- model.summary$qstat
heter.df <- do.call(rbind, heter.sum[1:3])
saveRDS(meta.coef, file = paste0(save.path,"Results/Results.2nd.1.est.rds"))
saveRDS(meta.vcov, file = paste0(save.path,"Results/Results.2nd.2.vcov.rds"))
saveRDS(heter.df, file = paste0(save.path,"Results/Results.2nd.3.heterogeneity_across_counties.rds"))
