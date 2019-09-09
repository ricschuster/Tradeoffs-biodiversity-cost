library(unmarked)
library(MuMIn)

library(rgdal)
library(rgeos)
library(foreign)

library(foreach)

source("f.AIC_cut.occu.sig.used.15.10.23.R")
source("getDesign.R")

ObsCovarSet <- c("DAY", "TIME", "COUNT_TYPE", "EFFORT_HRS", "EFFORT_DISTANCE_KM", "NUMBER_OBSERVERS")

covar <- read.csv("./input/eBird_occu_covar_site_obs.csv")
covar <- covar[,-c(3:22)]
scov <- read.csv("./input/eBird_landcape_covariates_scaled.csv")

#correlation matrix suggests removing Urban
scov <- scov[,-which(names(scov) %in% c("URB_1km","URB_100m"))]
SiteCovarSet <- names(scov)[-1]
covar <- cbind(covar,scov[,-1])

########################################################
########################################################
##### To go in loop
########################################################
########################################################
files <- list.files("./input",pattern="_eBird")

counter= 33

bird <- read.csv(sprintf("./input/%s",files[counter]))

tt <- strsplit(files[counter],"[.]")
bird_name <- unlist(tt[[1]])[1]

NewTab <- cbind(bird,covar[,-1])
MyOccuData <- formatWide(dfin=NewTab, type="unmarkedFrameOccu")

 #
 #Now fit a few different models for fun.
um.frame <- MyOccuData
fm0 <- occu(~1 ~1, data=um.frame)
#rm(um.frame)

tt <- paste(SiteCovarSet,collapse="+")
form <- as.formula(paste("~. ~. + ", tt))

#maybe not needed
#fm.og <- update(fm0,form)
#fm.d <- f.AICc.occu.sig(fm.og, ObsCovarSet, max.iter=30, detocc = 2, AICcut = 2, p.crit = 0.15, flag_se = TRUE)

fm.d0 <- f.AICc.occu.sig(fm0, ObsCovarSet, max.iter=30, detocc=2, AICcut=2, um.frame=MyOccuData)

dd <- occu.subset(fm.d0$modlst,cutoff=7)

if (length(dd)>1){
  aa <- model.avg(dd)
  tt <- attributes(aa$importance)$names[grep("p\\(", attributes(aa$importance)$names)]
  tt2<-substring(tt, 3)

  #finally, the names of the variables to include for modelling the occu side
  tt2 <- substr(tt2, 1, nchar(tt2)-1)
} else {
  tt2 <- strsplit(sub("~","",as.character(dd[[1]]@formula[2])),"\\s\\+\\s")[[1]]
}

tt2 <- paste0(tt2,sep="+",collapse="")
tt2 <- substr(tt2, 1, nchar(tt2)-1)

##################################################

form.p <- as.formula(paste("~. + ", tt2, "~. "))
#update null model with p variables from above
#equals occu null model
fm.p <- update(fm0, form.p)
flush.console()
##################################################
#start occu part of the models
f.occ <- f.AICc.occu.sig(fm.p, SiteCovarSet, max.iter=20, AICcut=2, ,um.frame=MyOccuData)

save.image(sprintf("%s_01.RData",bird_name))

f.occ.sig <- occu.rem.non.sig(f.occ,p.crit=0.15)
f.occ.sig.sub <- occu.subset(f.occ.sig,cutoff=7)

if (length(f.occ.sig.sub)>1){
  f.occ.avg <- model.avg(f.occ.sig.sub)
  f.occ.avg <- add.coefmat(f.occ.avg)
  designMats <- getDesign(MyOccuData, formula(f.occ.avg))
  fixed <- f.occ.avg$coefmat[,1]
  se <- f.occ.avg$coefmat[,2]
  est <- cbind(fixed,se)
  psi <- as.data.frame(est[grep("psi\\(",row.names(est)),])
  p <- as.data.frame(est[grep("p\\(",row.names(est)),])


} else {
  f.occ.avg <- f.occ.sig.sub
  designMats <- getDesign(MyOccuData, formula(f.occ.avg[[1]]))
  fixed <- coef(f.occ.avg[[1]])
  se <- sqrt(diag(vcov(f.occ.avg[[1]])))
  est <- cbind(fixed,se)
  psi <- as.data.frame(est[grep("psi\\(",row.names(est)),])
  p <- as.data.frame(est[grep("p\\(",row.names(est)),])
}


X <- designMats$X; V <- designMats$V; y <- designMats$y
removed <- designMats$removed.sites
X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
if(is.null(X.offset)) {
    X.offset <- rep(0, nrow(X))
}
if(is.null(V.offset)) {
    V.offset <- rep(0, nrow(V))
}

y <- truncateToBinary(y)
J <- ncol(y)
M <- nrow(y)

if (length(f.occ.sig.sub)>1){
  coefmat.compl <- f.occ.avg$coefmat
  coefmat <- coefmat.compl[,1]
} else {
  coefmat <- coef(f.occ.avg[[1]])
#  ff.se <- sqrt(diag(vcov(f.occ.avg[[1]])))
  ff.se <- vcov(f.occ.avg[[1]])
  ff.se <- diag(ff.se)
  ff.se <- sqrt(ff.se)
  ff.z <- coefmat/ff.se
  ff.p <- 2*pnorm(-abs(ff.z))
  coefmat.compl <- data.frame(Estimate =coefmat,
                              SE = ff.se,
                              z = ff.z,
                              p=ff.p)
  row.names(coefmat.compl) <- names(coefmat)
}
write.csv(coefmat.compl,sprintf("./estimates/%s.csv",bird_name))

rm(f.occ,f.occ.sig,fm.d0)
gc()
save.image(sprintf("%s_02.RData",bird_name))
file.remove(sprintf("%s_01.RData",bird_name))
