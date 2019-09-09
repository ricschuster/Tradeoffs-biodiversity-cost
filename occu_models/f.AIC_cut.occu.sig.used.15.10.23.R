f.AICc.occu.sig <- function(start.model, blocks, max.iter=NULL, detocc = 1, AICcut = 1, p.crit = 0.15, flag_se = FALSE, um.frame="", print.log=TRUE,...){
	# f.AICc.occu.sig: a function for "stepwise" regression using occupancy models of package unmarked 
	# start.model: initial model e.g. occu(~1~1, UMF)
	# detocc: if set to 1 (default) runs the function on the occupancy side; 2 does the detectability
	# foreback: if TRUE covariates are allowed to be dropped; if FALSE only adding covariates allowed
  # Parts based on Forward.lmer by Rense Nieuwenhuis (http://www.rensenieuwenhuis.nl/r-sessions-32/)
  # with some additions by for Nick Isaac 
	# Author: Richard Schuster (mail@richard-schuster.com 
	# 08 October 2012

  modlst <- c(start.model)
  x <- 2
  if (detocc == 1) coeff <- length(start.model@estimates@estimates$state@estimates)
  else coeff <- length(start.model@estimates@estimates$det@estimates)                                                                               

  best <- FALSE
	model.basis <- start.model
  keep <- list(start.model)
  AICmin <- AIClst <- AICc(start.model)
  #cutoff for when to exclude values
  cutoff <- 20	

  #critical z value alhpa = 0.3 (1.036433); 0.15 (1.439531); 0.1 (1.644853); 0.05 (1.959964)
  zc = qnorm(1-(p.crit/2))


	# Maximum number of iterations cannot exceed number of blocks, but this is also the default
	if (is.null(max.iter) | max.iter > length(blocks)) max.iter <- length(blocks)

	# Setting up the outer loop
	for(ii in 1:max.iter){
		models <- list()
    coeff <- coeff + 1
    cnt <- 1
    for (jj in 1:length(keep)) {
      # Iteratively updating the model with addition of one block of variable(s)
  		for(kk in 1:length(blocks)) {

        #check if blocks[kk] is already in the model, if so go to next kk
        if (detocc == 1){
          if(blocks[kk] %in% names(keep[[jj]]@estimates@estimates$state@estimates)) next
        } else {
          if(blocks[kk] %in% names(keep[[jj]]@estimates@estimates$det@estimates)) next
        }


        if (detocc == 1) form <- as.formula(paste("~. ~. + ", blocks[kk]))
        else form <- as.formula(paste("~. + ", blocks[kk], "~. "))

        if (class(dummy <- try(update(keep[[jj]], form, data=um.frame))) == "unmarkedFitOccu") {
          flag <- 0

          #check if model converged
          if(dummy@opt$converge!=0) flag <- 1

          #check if there is any NAN's in the SE's of the occupancy side
          if (flag == 0 && detocc == 1) {
            if(any(diag (vcov(dummy@estimates@estimates$state)) <0) ||
               any(sqrt(diag (vcov(dummy@estimates@estimates$state))) > cutoff)) {
                flag <- 1
            }
          }

          #check if there is any NAN's in the SE's of the detection side
          if (flag == 0 && detocc == 2) {
              if(any(diag (vcov(dummy@estimates@estimates$det))< 0) ||
                 any(sqrt(diag (vcov(dummy@estimates@estimates$det))) > cutoff)) {
                flag <- 1
            }
          }

          #check for repeat models
          if (flag == 0) {
            for (bb in 1:length(AIClst)) {
              if (round(AICc(dummy),digits=6) == round(AIClst[bb], digits=6)) {
                flag <- 1
                break
              }
            }
          }

          #Remove models where z < zc or SE = 0 for a beta > 0
          if (flag == 0 && flag_se == TRUE) {
            if (detocc == 1){
              if(any(abs(dummy@estimates@estimates$state@estimates[-1]/
                         sqrt(diag(vcov(dummy@estimates@estimates$state)))[-1]) <
                     zc)) flag <- 1
            } else {
              if(any(abs(dummy@estimates@estimates$det@estimates[-1]/
                         sqrt(diag(vcov(dummy@estimates@estimates$det)))[-1]) <
                     zc)) flag <- 1
            }
          }
        }
        else {flag <- 1}

        #add dummy model to the model list if it passes all previous tests
        if (flag == 0) {
            models[[cnt]] <- dummy
            modlst[[x]] <- models[[cnt]]
            AIClst <- c(AIClst, AICc(models[[cnt]]))
            x <- x + 1
            cnt <- cnt + 1
        }
      }
    }

    if(length(LL <- unlist(lapply(models, function(x) {AICc(x)}))) == 0) break
    keep <- list()
    k = 1
    cont <- 0
    #check for improvement in AIC, if none stop loop
		for (mm in order(LL, decreasing=FALSE)) {

      if (LL[mm]< AICmin + AICcut) {
        if (detocc == 1) {
          if (length(models[[mm]]@estimates@estimates$state@estimates) == coeff) {
            keep[[k]]<- models[[mm]]
            k <- k + 1
            if (LL[mm]<AICmin) {
              AICmin <- LL[mm]
              cont <- 1
            }
          }
        }
        else {
          if (length(models[[mm]]@estimates@estimates$det@estimates) == coeff) {
            keep[[k]]<- models[[mm]]
            k <- k + 1
            if (LL[mm]<AICmin) {
              AICmin <- LL[mm]
              cont <- 1
            }
          }
        }
      }
      else break
    }
    rm(models)
    gc()
    if (length(keep) == 0) break                                      
  }  
  ## Create Model List
  fitlst <- fitList(fits = modlst)
  modsel <- modSel(fitlst,nullmod=NULL)

	## Return the gathered output
  return(list(model=model.basis, modlst=modlst, fitlst=fitlst, modsel=modsel))
}

##################################################
##################################################
##### Post-processing functions              #####
##################################################
##################################################

occu.subset <- function(model,cutoff=2) {
#subsetting the above function output so it can be used with
#model.avg from package MuMIn
#2015-08-04
  mdlst <- model
  f.occ.aic <- unlist(lapply(mdlst, function(x) AIC(x)))
  min.aic <-min(f.occ.aic)
  f.occ.delta <- f.occ.aic - min.aic
  return(mdlst[f.occ.delta<cutoff])
}

occu.rem.non.sig <- function(model,p.crit = 0.15){
  #critical z value alhpa = 0.3 (1.036433); 0.15 (1.439531); 0.1 (1.644853); 0.05 (1.959964)
  zc = qnorm(1-(p.crit/2))

  model <- model$modlst

  sig <- unlist(lapply(model, function(x) all(abs(x@estimates@estimates$state@estimates[-1]/
                               sqrt(diag(vcov(x@estimates@estimates$state)))[-1]) > zc)))

  return(model[sig])
}

##################################################
##################################################
##### Functions from package MuMIn           #####
##################################################
##################################################
.coefarr.avg <-
function(cfarr, weight, revised.var, full, alpha) {
	weight <- weight / sum(weight)
	nCoef <- dim(cfarr)[3L]
	if(full) {
		nas <- is.na(cfarr[, 1L, ]) & is.na(cfarr[, 2L, ])
		cfarr[, 1L, ][nas] <- cfarr[, 2L, ][nas] <- 0
		#cfarr[, 1L:2L, ][is.na(cfarr[, 1L:2L, ])] <- 0
		if(!all(is.na(cfarr[, 3L, ])))
			cfarr[ ,3L, ][is.na(cfarr[ , 3L, ])] <- Inf
	}

	avgcoef <- array(dim = c(nCoef, 5L),
		dimnames = list(dimnames(cfarr)[[3L]], c("Estimate",
		"Std. Error", "Adjusted SE", "Lower CI", "Upper CI")))
	for(i in seq_len(nCoef))
		avgcoef[i, ] <- par.avg(cfarr[, 1L, i], cfarr[, 2L, i], weight,
			df = cfarr[, 3L, i], alpha = alpha, revised.var = revised.var)

	avgcoef[is.nan(avgcoef)] <- NA
	return(avgcoef)
}


	.makecoefmat <- function(cf) {
		no.ase <- all(is.na(cf[, 3L]))
		z <- abs(cf[, 1L] / cf[, if(no.ase) 2L else 3L])
		pval <- 2 * pnorm(z, lower.tail = FALSE)
		cbind(cf[, if(no.ase) 1L:2L else 1L:3L],
			`z value` = z, `Pr(>|z|)` = zapsmall(pval))
	}

add.coefmat <- function(object){
	is.arm <- ncol(object$msTable) == 6L && (colnames(object$msTable)[6L] == "ARM weight")

	weight <- object$msTable[, if(is.arm) 6L else 5L]

	object$coefmat <- .makecoefmat(.coefarr.avg(object$coefArray, weight,
			attr(object, "revised.var"), TRUE, 0.05))
  return(object)
}

