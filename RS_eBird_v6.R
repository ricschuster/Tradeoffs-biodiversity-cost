library(unmarked)
library(plyr)
library(here)
library(sp)
library(rgeos)
library(rgdal)

source("Wes.scripts.to.source_rs_edit.R")

BWO_MJ_manual_clip <- read.csv(here("data","BC_WASH_OREGON_March_June_manual_clip.csv"),as.is=T,header=T)
BWO_MJ_study_area_no_atr <- read.csv(here('data',"BC_WASH_OREGON_March_June_Study_area_no_atr.csv"),as.is=T,header=T)
BWO_MJ_study_area <- BWO_MJ_manual_clip[BWO_MJ_study_area_no_atr$Sel == 1,]

PCA <- read.csv(here('data',"ebird_NPLCC_CWNA_Normal_1981_2010MSY_PCA_output_5_seasonal_axis.csv"))
alpha.final <- read.csv(here("data","March_June_NPLCC_red_study_area_Sorted_Counts.csv"))

FilteredLists <- FilterChecklists(BWO_MJ_study_area,MaxHrs = 1.5, MaxDistKm = 5, MinProportion = 0.01)

MinNPerSiteObsrvr <- 1
MaxNPerSiteObsrvr <- 20
FilteredLists2 <- RetainRangeCountsPerLocObsvr(FilteredLists, MinNPerSiteObsrvr, MaxNPerSiteObsrvr)

BirdPlusPCA.all <- merge(x=FilteredLists2[,1:22], y=PCA, by="LOC_ID",all.x=T)
BirdPlusPCA.all <- data.frame(BirdPlusPCA.all,FilteredLists2[,c(1,23:ncol(FilteredLists2))])
BirdPlusPCA.all <- subset(BirdPlusPCA.all, select=-c(LOC_ID.1))

#only keep species that we can map
toMatch <- c(names(BirdPlusPCA.all)[1:43],gsub(" ","_",as.character(alpha.final$Species)),"NLists")
BirdPlusPCA.all.red <- BirdPlusPCA.all[,grep(paste(toMatch,collapse="|"), colnames(BirdPlusPCA.all))]

# rename bird columns to match 4 digit codes
colnames(BirdPlusPCA.all.red)[44:ncol(BirdPlusPCA.all.red)-1] <-  as.character(alpha.final$Alpha.Code[match(colnames(BirdPlusPCA.all.red)[44:ncol(BirdPlusPCA.all.red)-1],
                                                                gsub(" ","_",as.character(alpha.final$Species)))])
#Prepare unmarked files
MaxNPerSiteObsrvr <- 10
SiteIDName <- "ID"
ResponseName <- "BRCR"
SiteCovarSet <- c("NLists", names(PCA)[-1])
ObsCovarSet <- c("OBSERVER_ID","DAY", "TIME", "COUNT_TYPE", "EFFORT_HRS", "EFFORT_DISTANCE_KM", "NUMBER_OBSERVERS")

#Not sure why LOC_ID and OBSERVER_ID combination in Wes's script
# modfied in his script: Wes.scripts.to.source_rs_edit.R)
ID <- as.data.frame(BirdPlusPCA.all.red$LOC_ID, stringsAsFactors=F)
names(ID) <- "ID"
BirdPlusPCA.all.red <- cbind(ID, BirdPlusPCA.all.red)

resp <- colnames(BirdPlusPCA.all.red)[45:ncol(BirdPlusPCA.all.red)-1][order(colnames(BirdPlusPCA.all.red)[45:ncol(BirdPlusPCA.all.red)-1])]

########################################################################
# Resolve close locations
########################################################################
#what distance to lump observations under
cut.dist <- 100

projA <- CRS("+proj=utm +zone=10 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
aa <- SpatialPointsDataFrame(coords=BirdPlusPCA.all.red[,c("X","Y")],data=BirdPlusPCA.all.red, proj4string = CRS("+proj=longlat"))
dd_UTM <- spTransform(aa, CRSobj=projA)
coord.df.orig <- coord.df <- data.frame(LOC_ID=dd_UTM@data$LOC_ID,dd_UTM@coords)

run <- 1

#for(ii in 1:184){
while(run) {

  tt <- dist(coord.df[c("X","Y")])
  attributes(tt)$Labels <- coord.df$LOC_ID
  tt2<-as.matrix(tt)
  tt3 <- tt2[!duplicated(coord.df$LOC_ID),!duplicated(coord.df$LOC_ID)]
  tt3[lower.tri(tt3,diag=T)] <- NA
  tt3[tt3>cut.dist] <- NA

  lng <- vector()
  for (ii in 1:length(tt3[,1])){
    lng[ii]<-sum(!is.na(tt3[ii,]))
  }

  tt5 <- as.data.frame(tt3[lng>0,])
  if (sum(lng)>1){

    #remove NA's
    tt6.list <- Map(Filter, list(Negate(is.na)), split(tt5, seq(nrow(tt5))))
    pnts <- c(rownames(tt6.list[[1]]),colnames(tt6.list[[1]]))

  } else {

    pnts <- c(rownames(tt3)[lng>0],rownames(tt5)[!is.na(tt5)])
    run <- 0
  }


  lala <- coord.df[coord.df$LOC_ID %in% pnts,]
  lala <- lala[!duplicated(lala$LOC_ID),]

  if(length(pnts) >= 2){
    #find midpoint between 2 points or centriod
    X.new <- mean(lala$X)
    Y.new <- mean(lala$Y)
    coord.df[coord.df$LOC_ID %in% pnts,"X"] <- X.new
    coord.df[coord.df$LOC_ID %in% pnts,"Y"] <- Y.new
    coord.df[coord.df$LOC_ID %in% pnts,"LOC_ID"] <- lala$LOC_ID[1]
  } else {
    #here should be an error message
  }
}

#update LOC_ID, X and Y
BirdPlusPCA.all.red$LOC_ID <- BirdPlusPCA.all.red$ID <- as.character(coord.df$LOC_ID)
BirdPlusPCA.all.red$X <- BirdPlusPCA.all.red$LATITUDE <- coord.df$X
BirdPlusPCA.all.red$Y <- BirdPlusPCA.all.red$LONGITUDE <- coord.df$Y

mm <- data.frame(table(BirdPlusPCA.all.red$LOC_ID))
names(mm) <- c("LOC_ID","NLists2")
mm[,1]<-as.character(mm[,1])
BirdPlusPCA.all.red <- join(BirdPlusPCA.all.red,mm,by="LOC_ID")
BirdPlusPCA.all.red$NLists <- BirdPlusPCA.all.red$NLists2

birds <- list()
birds <- lapply(resp, function(x) MakeUnmarkedWideInput(InputData = BirdPlusPCA.all.red,
                                      SiteIDName = SiteIDName,
                                      ResponseName = x,
                                      SiteCovarSet = vector(),
                                      ObsCovarSet = vector(),
                                      MaxNPerSiteObsrvr = MaxNPerSiteObsrvr))

names(birds) <- resp
#used lapply instead
#for (ii in 1:length(resp)){
#    birds[[ii]] <- MakeUnmarkedWideInput(InputData = BirdPlusPCA.all.red,
#                                      SiteIDName = SiteIDName,
#                                      ResponseName = resp[ii],
#                                      SiteCovarSet = vector(),
#                                      ObsCovarSet = vector(),
#                                      MaxNPerSiteObsrvr = MaxNPerSiteObsrvr)
#}


ll <- data.frame(lapply(resp, function(x) rowSums(birds[[x]][,2:ncol(birds[[x]])],na.rm=T)))
names(ll) <- resp
ll[ll>1] <- 1

presence.df <- data.frame(LOC_ID=birds[[1]][,1],ll)
write.csv(presence.df,"presence.df.complete.csv",row.names=F)

#Now run the function.
NewWideTable <- MakeUnmarkedWideInput(InputData = BirdPlusPCA.all.red,
                                      SiteIDName = SiteIDName,
                                      ResponseName = "BRCR",
                                      SiteCovarSet = SiteCovarSet,
                                      ObsCovarSet = ObsCovarSet,
                                      MaxNPerSiteObsrvr = MaxNPerSiteObsrvr)

covar <- NewWideTable[,-(2:11)]
setwd("./occu_data")
write.csv(covar,"eBird_occu_covar_site_obs.csv",row.names=F)

for (ii in 1:length(birds)){
  write.csv(birds[[ii]],sprintf("%s_eBird.csv",names(birds)[ii]),row.names=F)
}









