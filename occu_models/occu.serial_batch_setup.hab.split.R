setwd("D:/Work/NPLCC/occu")
files <- list.files("./input",pattern="_eBird")
birds <- sub(".csv","",files)
input <- readLines("./batch_hab/occu.serial_batch_setup.hab.R")

batches <- 10
cutoff <- ceiling(length(birds)/batches)
bt <- 1
rw <- 0
outlst <- as.list(rep(0,batches))
txt <- "counter=1"

for(ii in 1:length(files)){
  if (ii%%cutoff) {
    rw <- rw + 1
  } else {
    write(outlst[[bt]],sprintf("./batch_hab/aa_batch%s.sh",bt))
    bt <- bt + 1
    rw <- 1
  }
  input_mod <- input
  input_mod[grep(txt, input_mod)]<-paste("counter=",ii)
  write(input_mod,sprintf("./batch_hab/%s.R",birds[ii]))
  outlst[[bt]][rw] <- sprintf("R CMD BATCH --no-save /SciBorg/array0/schuster/batch_hab/%s.R",birds[ii])
}


