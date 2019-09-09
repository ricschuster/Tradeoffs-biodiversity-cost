#This R script contains the work needed for the third part of the outlined
# manuscript on occupancy modelling with eBird data: the actual implementation
# and fitting of occupancy models to data from eBird, with specific attention
# being paid to whether and how estimates vary with the number of checklists
# from a location.  This emphasis is placed to see whether anything specific
# needs to be done in order to take advantage of the locations with one-off
# checklists.
#As part of this work, I want to do cross-validation of model fitting, which
# will mean (I think) creating a function to do the cross-validation
# model fit calculations...and now that I think about it also a function to
# select the subsample to be withheld fo testing.
#This script was started on 10 january 2015 by Wesley Hochachka, with bits
# borrow from other files used for more preliminary looks at these data.

#######################################################################################
########### Multi-Step Filtering of ERD Data Further After Inputting into R ###########
#######################################################################################

##################################################################
# Filter Out Unwanted Data by Column Values & Low Row Prevalence #
##################################################################
#Define a function that filters the input data, retaining only those rows (checklists)
# and columns (species) that the analyst deems suitable for further use.  This
# filtering is done based on the following criteria:
#  (1) Count protocol: hardwired removal of rows not from stationary or travelling counts
#  (2) Count duration: removal of counts longer in time than user-specified maximum
#  (3) Count distance: removal of counts longer than user-specified distance
#  (4) Species prevalence: removal of columns for any species represented in a
#       smaller proportion of checklists than specified by the user.
# As inputs, this function requires the name of the input data.frame, the
# threshold maximum count duration in hours, the maximum threshold count
# distance travelled in kilometers, and the minimum proportion (i.e. zero to one
# value) of checklists on which a species is reported.
#Additionally, make a slight change to content of the field describing count
# protocols, by turning the code values into human-intelligible labels.
FilterChecklists <- function(InputData, MaxHrs, MaxDistKm, MinProportion){
  #The first task for this function is to remove data from protocols other than
  # stationary ("P21") and travelling ("P22"), and counts from longer times or
  # distances than specified in the input variables.  This filtering is done by
  # appending a variable to be treated as boolean to the input data.frame, and
  # then running each row through a set of 3 tests that will flip the boolean
  # from "keep" to "remove" states (i.e. from the default value of 1 to zero).
  Keep <- rep(1, nrow(InputData))
  InputData <- cbind(InputData, Keep)
  #
  InputData[InputData$COUNT_TYPE != "P21" & InputData$COUNT_TYPE != "P22",]$Keep <- 0
  #
  InputData[InputData$EFFORT_HRS > MaxHrs,]$Keep <- 0
  #
  InputData[InputData$EFFORT_DISTANCE_KM > MaxDistKm,]$Keep <- 0
  #
  Filtered_1 <- InputData[InputData$Keep == 1,]
  #
  #Now change the labels for the count protocols to make sense to mere mortals.
  Filtered_1[Filtered_1$COUNT_TYPE == "P21",]$COUNT_TYPE <- "stationary"
  Filtered_1[Filtered_1$COUNT_TYPE == "P22",]$COUNT_TYPE <- "travelling"
  #
  #Now that we have trimmed down the number of rows in the data, the next step is
  # to trim out columns that are not needed, and specifically columns of data from
  # species that have such a low rate of reporting as to make them uninteresting
  # for use.  We do this by moving through the columns of bird records and one by
  # one count the number of checklists with each species and see if the proportion of
  # checklists with the focal species is at or above the minimum threshold for
  # retention.  The data from each species that exceeds this threshold are appended
  # to the base data file of count-related information.  In the process, counts of
  # birds are converted to presence/absence information.  Because many species have
  # checklists in which a species is just reported as being present but without any
  # counts (a value of "X") many of the columns of ostensibly numeric information
  # are being treated by R as character data, so part of this processing changes the
  # data back to numeric.  NOTE: for abundance modelling, this part of the function
  # would need to be modified to delete or otherwise appropriately deal with records
  # in which the species was reported as merely present.
  #NOTE: This part of the function is somewhat fragile, because it is hard-coded to
  # assume that the first 19 columns of the input data.frame contain all of the
  # basic information describing the count process (e.g., location, date, time, effort),
  # and columns 20 to the end as data on the counts of individual bird species (one
  # per column).  This is the format of the 2013 version of the ERD data.  However,
  # future modifications might cause changes in the number of initial columns of
  # count-descriptor information, which would cause the function to behave
  # inappropriately.
  ProportionDenominator <- nrow(Filtered_1)
  Filtered_2 <- Filtered_1[,1:19]
  #
  for (i in 20:(ncol(Filtered_1) - 1)){
    FocalColumn <- as.character(Filtered_1[,i])  #insure character data to start
    FocalColumn[FocalColumn != "0"] <- "1"       #give value of 1 to all non-zero counts
    FocalColumn <- as.numeric(FocalColumn)       #convert data to numeric
    if (sum(FocalColumn)/ProportionDenominator >= MinProportion){     #test if >threshold prevalence
      Filtered_2 <- cbind(Filtered_2, FocalColumn)                    #append new column to data.frame
      names(Filtered_2)[ncol(Filtered_2)] <- names(Filtered_1)[i]     #name new column with species name
    }
  }
  #Send completed data.frame back to main program.
  return(Filtered_2)
}

########################################################
# Filter Out Unwanted Counts Per Site and Rare Species #
########################################################
#Now we have to deal with the instances in which an observer has submitted a larger
# than MaxNPerSiteObsrvr number of checklists.  The approach that I am taking is to
# select a random subset of these lists to use, discarding the others.  We start by
# splitting the data into two groups: those person-site combinations with up to the
# the maximum number of checklists, and those with more than the maximum.  Then, for
# this latter group of person-site combinations we assign a random number to
# each of the records, and finally delete all ranks above the maximum number that
# we want to keep.  Unfortunately, R's "rank" function is global in scope, not
# allowing us to easily assign ranks to each separate observer-location combination.
# As best as I can tell we can't easily take advantage of an SQL-based solution
# because RSQLite also does not have a "rank" function, so we need to use
# "aggregate" to do the job for us.  Using "aggregate" is somewhat painful,
# because "aggregate" returns a separate list of ranks for each person-location-year,
# and I need to convert these to a vector using the "unlist" function before
# adding the new column and removing the rows with rank values over the maximum
# allowed.  Additional funkiness is that aggregate seems to be ordering its output
# in alphabetical order of the factor numbers of a factor version of the grouping
# variable, PersonLocYr, and not of the character-value alphabetical order.  So,
# I need to insure that the data in DataSubset_ExcessLists are ordered by this
# same way.  Only then can I concatinate the two subsets of rows back together
# into a single data.frame, and then recalculate the NTimesPerLoc values for
# this subset of the data.  From this we can finally build the data objects
# for occupancy modelling.
#While we are at it, we will also remove the data from location-observer
# combination with fewer than the minimum number of lists in the set of data.
#All of this is wrapped into function so that it can be applied to multiple,
# independent subsets of data.
RetainRangeCountsPerLocObsvr <- function(InputData, MinNPerSiteObsrvr, MaxNPerSiteObsrvr){
  #
  #First, we need to add a count of the number of checklists in the data set for each
  # combination of location (LOC_ID) and observer (OBSERVER_ID).  The four steps are:
  # (1) use the "aggregate" function to count the times that each LOC_ID - OBSERVER_ID
  # combination was seen, (2) name the columns in the output for merging, (3) merge the
  # new information to append the counts to the ends of the original data lines, and
  # (4) tidy up.
  #  NPerLoc <- aggregate(InputData$LOC_ID, by=list(InputData$LOC_ID, InputData$OBSERVER_ID), FUN="length")
  NPerLoc <- aggregate(InputData$LOC_ID, by=list(InputData$LOC_ID), FUN="length") #RS edit
  names(NPerLoc) <- c("LOC_ID", "NLists")
  InputData <- merge(x=InputData, y=NPerLoc, by="LOC_ID")
  rm(NPerLoc)
  #
  #Now, divide the data into two subsets based on whether we need to remove excess lists.
  DataSubset_EnoughLists <- InputData[InputData$NLists <= MaxNPerSiteObsrvr,]
  DataSubset_ExcessLists <- InputData[InputData$NLists > MaxNPerSiteObsrvr,]
  #
  #For the subset that needs trimming, add a column of random numbers to be converted
  # into ranks, and then sort the original data by a factor-coded version of the location
  # ID and observer ID to insure that we can concatinate the column of rank information back
  # onto the original data.  The data should already be sorted appropriately, so this
  # sorting is a bit of paranoid self defense.
  DataSubset_ExcessLists <- cbind(DataSubset_ExcessLists, runif(n=nrow(DataSubset_ExcessLists)))
  names(DataSubset_ExcessLists) <- c(names(InputData), "RandNum")
  DataSubset_ExcessLists <- DataSubset_ExcessLists[order(as.factor(DataSubset_ExcessLists$LOC_ID), as.factor(DataSubset_ExcessLists$OBSERVER_ID)),]
  #
  #Create the ranks based on location and observer groups, and then extract this information from
  # the list object that is returned by "aggregate" in the form of a vector using "unlist".  Note
  # that even though the data have been sorted first by LOC_ID and then by OBSERVER_ID, these two
  # sorting variables are given in reverse order in the "by" option of aggregate.  We need to do
  # this, because the list produced as output to "aggregate" is sorted in reverse order to the
  # order in which the variables are stated in the "by" list.
  RanksList <- aggregate(x=DataSubset_ExcessLists$RandNum,
                         by=list(DataSubset_ExcessLists$OBSERVER_ID, DataSubset_ExcessLists$LOC_ID),
                         FUN="rank", simplify=T)
  Ranks <- unlist(RanksList$x)
  #
  #Join the column of ranks back onto the original data, and remove those records with
  # ranks greater than the maximum desired number of records per location-observer.
  # Remove the two columns (random number and rank) that are no longer needed, and
  # add back the records onto the subset of the data that was not in need of removing
  # excess records.
  DataSubset_ExcessLists <- cbind(DataSubset_ExcessLists, Ranks)
  DataSubset_ExcessLists <- DataSubset_ExcessLists[DataSubset_ExcessLists$Ranks <= MaxNPerSiteObsrvr,]
  DataSubset_ExcessLists <- DataSubset_ExcessLists[, 1:(ncol(DataSubset_ExcessLists) - 2)]
  DataSubset_new <- rbind(DataSubset_EnoughLists, DataSubset_ExcessLists)
  #
  #Tidy up.
  rm(RanksList, Ranks, DataSubset_EnoughLists, DataSubset_ExcessLists)
  #
  #Trim off the data from location-observer combinations with fewer than specified minimum
  # number of lists.
  DataSubset_new <- DataSubset_new[DataSubset_new$NLists >= MinNPerSiteObsrvr,]
  #
  #The counts NLists are now wrong for any of the site-observer combinations
  # for which only a random subset of size MaxNPerSiteObsrvr was retained.
  # So, we need to recount the numbers of checklists per site-observer
  # combination, before returning the updated data set for further use.
  NPerLoc <- aggregate(DataSubset_new$LOC_ID,
                       #                       by=list(DataSubset_new$LOC_ID, DataSubset_new$OBSERVER_ID),
                       by=list(DataSubset_new$LOC_ID), #RS edit
                       FUN="length")
  names(NPerLoc) <- c("LOC_ID", "NLists_new")
  DataSubset_new <- merge(x=DataSubset_new, y=NPerLoc, by="LOC_ID")
  DataSubset_new <- DataSubset_new[,-(which(names(DataSubset_new) == "NLists"))]
  names(DataSubset_new)[which(names(DataSubset_new) == "NLists_new")] <- "NLists"
  rm(NPerLoc)
  #Ship out the final data.frame for use outside of the function.
  return(DataSubset_new)
}
##################################################
# Retain Columns of Counts for Subset of Species #
##################################################
#
#Create a function to simplify the extraction of columns from the larger
# data.table so that only data from pre-specified list of species are retained,
# along with all of the count-description information.  As part of this, we
# will create a new, single-column identifier that denotes a unique combination
# of location and observer.  We make this new column as a simple concatenation
# of the LOC_ID and OBSERVER_ID values.  This is another function that is somewhat
# fragile, in the sense that we are hardwiring in the assumption that the first
# seventeen columns of the input data contain all of the count-description data;
# while this is true for the ERD 2013 data this could change in the future.
FilterSelectSpecies <- function(InputData, Species){
  #Create reduced-column data.frame with the count-description data and the count
  # data from the first species in the list of desired species.  This will serve as
  # the base on which to append additional columns (if needed) of count data from
  # other species in that list.
  #
  #Identify the column in which the first focal species' count data are located.
  SppColIndex <- which(names(InputData) == Species[1])
  #
  #Now bind together three "groups" of columns: (1) the basic count-description data,
  # (2) the tally of the number of records for each location-observer combination (the
  # last column of the input data.frame), and (3) the count data from the first (and
  # possibly only) species in the list of species whose count data are to be extracted.
  # Once the base of the output data.frame is created, label with columns with their
  # appropriate names and tidy up.
  Output <- data.frame(InputData[,1:17], InputData[,ncol(InputData)], InputData[,SppColIndex])
  names(Output) <- c(names(InputData[,1:17]),
                     names(InputData[ncol(InputData)]),
                     names(InputData[SppColIndex]))
  rm(SppColIndex)
  #
  #Now loop through the rest of the species in the list of desired species (testing
  # if there is more than one species!), and append additional columns of species'
  # data one by one.
  if (length(Species) > 1){
    for (i in 2:length(Species)){
      SppColIndex <- which(names(InputData) == Species[i])
      Output <- cbind(Output, InputData[,SppColIndex])
      names(Output) <- c(names(Output[,1:(length(Output)-1)]), names(InputData[SppColIndex]))
    }
    rm(SppColIndex, i)
  }
  #As a final step add to the front of the output data.frame a 1-column
  # unique ID value for each combination of observer and location, and
  # add this to the data to be output from this function.
  #   ID <- as.data.frame(paste(Output$LOC_ID, Output$OBSERVER_ID, sep="_"), stringsAsFactors=F)
  ID <- as.data.frame(Output$LOC_ID, stringsAsFactors=F) #RS edit
  names(ID) <- "ID"
  Output <- cbind(ID, Output)
  #
  #Return this data.frame.
  return(Output)
}


#######################################################################################
################# Creating Data.Frame Suitable for Use With unmarked ##################
#######################################################################################

#Once the primordial data.frames have been created, we need to take the one-observation-
# per-row data, and turn them into a form in which the repeated observations from a
# single location form a series.  That is the work for this next section of the script.
# Our presumed target will be a data.frame object that can be fed into the unmarked
# package.  While it is possible for unmarked to use one-observation-per-row data,
# using unmarked's "formatLong" function, this is only possible of there are
# no site-level covariates (i.e. predictors specific to a site, e.g., elevation, that
# do not change from observation to observation within that site).  Because this
# is not always the case, the script below will create a data.frame that is more
# universally useable, to be processed by the "formatWide" function of unmarked
# for analysis.
#
#load the unmarked library, and set the maximum number of observation
# periods for data from each site.

#Create a function to produce an output file in the format of an unmarked "wide"
# input data.frame.
MakeUnmarkedWideInput <- function(InputData, SiteIDName, ResponseName, SiteCovarSet, ObsCovarSet, MaxNPerSiteObsrvr){
  #Create single-item data.frame from which to build units of repeated counts
  # and their covariates.
  BuildingBlock <- data.frame(x = NA)
  #
  #From BulidingBlock create a template for the site identifier, and another for
  # the vectors of response variable and observation-specific covariations that
  # have the length of the number of maximum number of possible observation
  # periods.
  SingletonTemplate <- BuildingBlock
  #
  GroupTemplate <- BuildingBlock
  for (i in 2:MaxNPerSiteObsrvr){
    GroupTemplate <- cbind(GroupTemplate, BuildingBlock)
  }
  rm(BuildingBlock)
  #
  #Copy these templates to form the base on which to build a complete
  # record of the response plus observation-level covariates, filled
  # with NA values.  For use with unmarked, the format for data should
  # be the following:
  #  (i)   The site-identifier variable (optional) named "site"
  #  (ii)  The sets of response-variable values for all observations at
  #         that site.
  #  (iii) Columns for site-level covariates (i.e. covariates whose values
  #         do not change among observation periods)
  #  (iv)  Sets of columns for the observation-level covariates (i.e.
  #         covariates that can change in value among the observation
  #         periods within each site).  For each covariate there needs to
  #         be a set equal in length to the maximum number of observation
  #         periods.  The names of covariates within each set have the same
  #         prefix, a dot, and then na integer number from 1 to the maximum
  #         number (e.g., Time.1, Time.2, Time.3 ... Time.n).
  # Our first step in building this structure will be to assemble together
  # parts (i) and (ii).  Because the site label/column is optional for data
  # intended for unmarked, we need to test if there is an identifier variable.
  # We do this by testing if the R object "SiteIDName" exists.
  if (exists("SiteIDName")){
    BlankRecord <- cbind(SingletonTemplate, GroupTemplate)
  } else if (!exists("SiteIDName"))
    BlankRecord <- GroupTemplate
  #
  #Assign names of individual items in the GroupTemplate, initially for the
  # response variable, following the naming convention that unmarked
  # will need: variable name, followed by a dot, followed by the integer
  # number of that observation in the sequence.  As before, we test to
  # make sure that we actually have a variable to be used as a site
  # identifier, and the counter used in the for-loop for creating the
  # column name labels deals with the possible absence of the site
  # identifier by using the "1+j" construct.
  if (exists("SiteIDName"))
    names(BlankRecord)[1] <- "site"
  ColNameBase <- "y"
  if (exists("SiteIDName")){
    j <- 1
  } else j <- 0
  for (i in (1+j):(MaxNPerSiteObsrvr+j)){
    names(BlankRecord)[i] <- paste(ColNameBase, (i-1), sep=".")
  }
  rm(i, j, ColNameBase)
  #
  #Now we add component (iii), the site-level covariates and name them.
  # The list of site-level covariates needs to have been specified already
  # in a vector of names SiteCovarSet.  Note that we are cutting a corner
  # a wee bit, but recycling and renaming the single variable in
  # SingletonTemplate for each of the site-level covariates.
  if(length(SiteCovarSet)){
    for (i in 1:length(SiteCovarSet)){
      ColName <- SiteCovarSet[i]
      names(SingletonTemplate) <- ColName
      BlankRecord <- cbind(BlankRecord, SingletonTemplate)
    }
    rm(i, ColName)
  }
  #
  #Now we add the parts of component (iv) to the structure. We add each
  # observation-level covariate in turn as its complete set and name the
  # observation-specific covariate variable.  The list of names of
  # observation-level covariates needs to have already been specified in
  # the vector ObsCovarSet.  Note that we are cutting corners a bit by
  # recycling and renaming the components of GroupTemplate for each of
  # the observation-level covariates in turn.
  if(length(ObsCovarSet)){
    for (i in 1:length(ObsCovarSet)){
      ColNameBase <- ObsCovarSet[i]
      for (j in 1:MaxNPerSiteObsrvr){
        names(GroupTemplate)[j] <- paste(ColNameBase, j, sep=".")
      }
      BlankRecord <- cbind(BlankRecord, GroupTemplate)
    }
    rm(i, j, ColNameBase, SingletonTemplate, GroupTemplate)
  }
  #
  #
  #Prepare the input data for formatting for unmarked.  First sort the data
  # in the order required for reformatting.
  InputData <- InputData[order(InputData$ID,
                               InputData$DAY,
                               InputData$TIME),]
  #
  #
  #Create a copy of the template to be filled and appended to the output table.
  NewWideRow <- BlankRecord
  #
  #Now we can assemble our data.frame containing one row of data for the
  # complete set of observations from each of the site-observer combinations.
  # We do this by using an outer "while" loop that moves through the input data
  # for each site-observer set of lines in turn.  For each new site-observer
  # combination the first data to population this new record are the site ID
  # variable and then information from the site-level covariates.  We then
  # step through the records of this site with a middle "for" loop that walks
  # through the lines that contain data from each observation periods for that
  # specific location-observer combination, making  use of NLists (the variable
  # that tells the number of observation periods) in the middle loop's counter
  # of the records for each site-observation combination.  Within each of the
  # input lines we need a final, inner-most "for" loop that steps through each
  # observation period's response (detected/not detected) information and the
  # observation-level covariates in turn within each of the lines specified
  # in the middle loop.
  #Note that this process only works when the input data are sorted by the
  # ID column (the concatenated information from LOC_ID and OBSERVER_ID).
  # While this should already be the case, out of an overabundance of paranoia
  # we will be sorting the data again, and in doing this sorting not only by
  # the ID variable but also by date and time of day within ID (mostly just
  # because it feels more orderly to have the final record presenting
  # observations in chronological order).
  #
  #We use "i" as the counter in the outer loop, working down through the first
  # lines of each group of input records, starting with the very first line of the
  # input file. We only stop when the next possible line number is beyond the
  # end of the file.  This starting and stopping is set and controlled by the
  # next two lines.
  i <- 1
  while (i <= nrow(InputData)){        #outer loop start: stepping through unique site/observer combinations
    if (exists("SiteIDName"))         #if there is site ID information, add to record being built
      NewWideRow[1] <- InputData[i, which(names(InputData) == SiteIDName)]
    if(length(SiteCovarSet)){
      for (s in 1:length(SiteCovarSet)){            #population record with site-level covariate information
        NewWideRow[1+MaxNPerSiteObsrvr+s] <- InputData[i, which(names(InputData) == SiteCovarSet[s])]
      }
      rm(s)
    }
    for (j in 2:(ifelse(InputData[i,]$NLists > MaxNPerSiteObsrvr, MaxNPerSiteObsrvr + 1,InputData[i,]$NLists + 1))){   #middle loop start: each observation in site-observer combination
      NewWideRow[j] <- InputData[i+j-2, which(names(InputData) == ResponseName)]   #response variable data
      if(length(ObsCovarSet)){
        for (k in 1:length(ObsCovarSet)){       #inner loop: each obs-level covariate within site-observer
          NewWideRow[j+length(SiteCovarSet)+(MaxNPerSiteObsrvr*k)] <- InputData[i+j-2, which(names(InputData) == ObsCovarSet[k])]
        }
      }
    }
    #      if (i == 1 & (j-1) == InputData[i,]$NLists){         #create new table if this is the first output record
    if (i == 1 & (j-1) == ifelse(InputData[i,]$NLists > MaxNPerSiteObsrvr, MaxNPerSiteObsrvr,InputData[i,]$NLists)){
      NewWideTable <- NewWideRow
      NewWideRow <- BlankRecord
      #      } else if (i > 1 & (j-1) == InputData[i,]$NLists){     #append latest record to output table, if it exists
    } else if (i > 1 & (j-1) == ifelse(InputData[i,]$NLists > MaxNPerSiteObsrvr, MaxNPerSiteObsrvr,InputData[i,]$NLists)){     #append latest record to output table, if it exists
      NewWideTable <- rbind(NewWideTable, NewWideRow)
      NewWideRow <- BlankRecord
    }
    i <- i + InputData[i,]$NLists     #increment the counter to the start of the next location-observer block
  }
  rm(i, j, k, NewWideRow)  #tidy up
  #
  #Any predictor variables that are represented by character data need to be converted into factors
  # for analysis.  While it might seem logical (to me at least) to convert the character variables
  # prior to assembling the unmarked wide-format data.frame, for some reason that I do not understand
  # R seems to be back-converting factor data to character in some cases (I think it has to do with
  # a column's data type being set in the first row, the presence of missing values not being
  # factors, and causing newly appended rows to have their data coerced back to character from
  # factors).  Anyway, as a result, we need to do a post-assembly coersion of all character data to
  # R's factor data type.  Here we do this by identifying the first column of predictor data (the one
  # immediately after the last "y.*" column), and then from there to the last column testing if each
  # column has the data type "character" and if so converting it to the factor data type.
  #Identify the column number for the first column of predictors
  PredictCol1 <- (which(names(NewWideTable) == paste("y.", MaxNPerSiteObsrvr, sep="")) + 1)
  #Loop through the columns in turn, and change character data types to factors.
  if(PredictCol1 <- length(NewWideTable)){
    for (i in PredictCol1:length(NewWideTable)){
      if (is.character(NewWideTable[,i]))
        NewWideTable[,i] <- as.factor(NewWideTable[,i])
    }
    rm(i)
  }
  #
  #Finally, we can return the table.
  return(NewWideTable)
}
