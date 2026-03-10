rm(list = ls())
gc()

library(dplyr)
library(mvmeta)
library(openxlsx)
library(dlnm)
library(splines2)


### 0. basic settings
data.path <- save.path <- "~/"


### 1. read in the association estimates
meta.coef <- readRDS(paste0(data.path,"Results/Results.2nd.1.est.rds"))
meta.vcov <- readRDS(paste0(data.path,"Results/Results.2nd.2.vcov.rds"))


### 2. ac data and temperature distribution
ac.fake <- readRDS(file = paste0(data.path, "Data/Data.Z.Fake_AC.rds"))
exp.dist <- readRDS(paste0(data.path, "Results/Results.0.Dist_all.rds"))
exp.dist <- as.data.frame(exp.dist)
temp.var <- "tavg"
temp.var.knots <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(10,75,90)),temp.var])
temp.boundary <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(0, 100)),temp.var])
prediction.pct.ofinterest <- 0:100 # rough calculation; could be set as seq(0,100,0.1) to get more precise estimates
temp.pred.seq <- exp.dist[which(round(as.numeric(exp.dist$pct),2) %in% prediction.pct.ofinterest),temp.var, drop=FALSE]
# round(, 2): floating-point precision issue
temp.median <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(50)),temp.var])
temp.mmt.boundary <- as.numeric(exp.dist[which(as.numeric(exp.dist$pct) %in% c(1, 99)),temp.var])



### 3. calculate cbtemp coef and var at different ac levels
## (1) predictor matrix
ac.ofinterest <- 1300 # any AC value of interest
cb.ac.ofinterest <- onebasis(ac.ofinterest, fun = "ns", knots = median(ac.fake), Boundary.knots = range(ac.fake))
cb.ac.ofinterest <- as.data.frame(cb.ac.ofinterest)
predictor.matrix <- cbind(1, cb.ac.ofinterest)
colnames(predictor.matrix)[1] <- "Intercept"
## (2) coefficients for cbtemp
cbtemp.coef <- as.matrix(predictor.matrix) %*% as.matrix(meta.coef)
## (3) variance for cbtemp.coef
num.temp.basis <- ncol(meta.coef)
num.var <- nrow(meta.coef)
cbtemp.coef.var <- matrix(NA, nrow = num.temp.basis, ncol = num.temp.basis)
for(row.i in 1:num.temp.basis){
  for (col.j in 1:num.temp.basis) {
    var.part <- meta.vcov[((row.i-1)*num.var+1):(row.i*num.var),
                           ((col.j-1)*num.var+1):(col.j*num.var)]
    cbtemp.coef.var[row.i,col.j] <- as.matrix(predictor.matrix) %*% as.matrix(var.part) %*% t(as.matrix(predictor.matrix))
  }
}
colnames(cbtemp.coef.var) <- colnames(cbtemp.coef)




### 4. predict and output
## (1) create the basis
cbtemp <- crossbasis(temp.pred.seq[,rep(1,times=7)], lag = 6, 
                     argvar=list(fun = "ns",
                                 knots = temp.var.knots, 
                                 Boundary.knots = temp.boundary),
                     arglag = list(fun = "ns",
                                   knots = 3))
## (2) prediction - median
pred.median <- crosspred(cbtemp, coef = as.vector(cbtemp.coef), vcov = as.matrix(cbtemp.coef.var),
                         model.link = "log", at=unlist(temp.pred.seq), cen=temp.median)
pred.est.median <- as.data.frame(cbind(prediction.pct.ofinterest, pred.median$predvar, pred.median$allRRfit, pred.median$allRRlow, pred.median$allRRhigh))
colnames(pred.est.median) <- c("temp.pct","predvar",c("OR","Low","High"))
#@ (3) prediction - MMT
pred.est.median.sub <- pred.est.median %>% 
  filter(temp.pct>=1 & temp.pct<=99)
MMT <- pred.est.median.sub$predvar[which.min(pred.est.median.sub$OR)]
MMT.pct <- pred.est.median.sub$temp.pct[which.min(pred.est.median.sub$OR)]
pred.mmt <- crosspred(cbtemp, coef = as.vector(cbtemp.coef), vcov = as.matrix(cbtemp.coef.var),
                      model.link = "log", at=unlist(temp.pred.seq), cen=MMT)
pred.est.mmt <- as.data.frame(cbind(prediction.pct.ofinterest,pred.mmt$predvar, pred.mmt$allRRfit, pred.mmt$allRRlow, pred.mmt$allRRhigh))
colnames(pred.est.mmt) <- c("temp.pct","predvar",c("OR","Low","High"))
pred.est.mmt$MMT <- MMT
pred.est.mmt$MMT.pct <- MMT.pct
## (4) plot
p <- ggplot()+
  geom_hline(yintercept = 1, color=grey(0.5))+
  geom_vline(xintercept = temp.mmt.boundary,colour=grey(0.5),linetype=2)+
  geom_smooth(data=pred.est.mmt, aes(x=predvar, y=OR, ymin=Low, ymax=High),
              stat="identity", alpha=0.1, linewidth = 0.4)+
  geom_point(data = pred.est.mmt %>% dplyr::select(MMT) %>% distinct(),
             aes(x = MMT, y=1), color="black", size=1)+
  scale_x_continuous(breaks = seq(-25, 35, 5)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  ylab("Odds ratio at lag 0-6 days")+
  xlab("Daily mean temperature (°C)")
ggsave(p, filename = paste0(save.path, "Results/Figure.Temp_mortality_curve_ac1300.tiff"),
       width = 6.5, height = 5, unit = "in") # 1 sec
