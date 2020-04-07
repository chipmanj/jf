# findClosestRecord.R
# J Chipman
#
# minFun: returns NA if all observations are NA; o/w minimum excluding NAs
# maxFun: returns NA if all observations are NA; o/w maximum excluding NAs
#
# findClosestRecord: Finds the closest observation to a set point in time
#         within a window of time.  Function also has capability of finding
#         the furthest record within a window of time.
#
# featureDat:        A long dataset of the feature of interest.
# index:             A dataset of all observations in cohort.  Each observation observed once.
# mergeby:           Key to merge featureDat and index.
# featureCol:        Column within featureDat from which the feature may be identified.
#                    featureCol may include many features though only one is of interest.
# feature:           The feature of interest within featureCol.
# featureEvaluation: The value of the feature.  Example the feature may be BodyHeight
#                    and the featureEvaluation may be 72 inches.
# featureDate:       When each feature was observed.
# featureName:       The name which will be saved to the dataset.  If it's a baseline feature,
#                    the feature name will be: Bl_[featureName].
# featureType:       Type of feature which may be helpful when customizing the code.
#                    (I may delete this column for more generalizability).
# t0Col:             Column identifying time zero.
# t0PlusDays:        Default = 0. Some applications will have the true time zero be some days after a certain event.
#                    Example 30 days after surgery.
# daysRangeLo:       Low end of window days.  When Null, observations can have occurred anytime previously without a limit.
# daysRangeHi:       Hi end of window days.  When Null, observations can have occurred anytime after without a limit.
# furthesFromT0:     Default = FALSE and finds closest record to t0Col + t0PlusDays.  When TRUE, finds furthest record.
# daysReduction:     Default = NA when there are not multiple, differing evalutions on the same day.
#                    "min", "max", or "median" finds minimum, maximum, or median of same day evaluations.

minFun <- function(x){
  if(sum(is.na(x))==length(x)){
    return(NA)
  } else {
    return(min(x,na.rm=TRUE))
  }
}

maxFun <- function(x){
  if(sum(is.na(x))==length(x)){
    return(NA)
  } else {
    return(max(x,na.rm=TRUE))
  }
}

findBlFollup <- function(featureDat, index, mergeby,
                              featureCol=NA, feature=NA, featureEvaluation=NA, featureDate,
                              featureName, featureType,
                              t0Col, t0PlusDays = 0,
                              daysRangeLo, daysRangeHi, furthestFromT0=FALSE,
                              dayReduction = NA){


 #### Pre-Processing ####
 ####--------------------

  # Keep only mergeby field as common field between featureDat and index datasets
  if(sum(is.element(colnames(featureDat), colnames(index))) > 1){
    commonFields <- which(colnames(featureDat) %in% colnames(index))
    keepMergeBy  <- which(colnames(featureDat) == mergeby)
    dropFields   <- commonFields[commonFields != keepMergeBy]
    featureDat   <- featureDat[,-dropFields]
  }

  # If type = gap90check, create a column for presence of any record
  if(featureType=="gap90check"){
    featureDat[,"ptEngage"] <- 1
    featureEvaluation       <- "ptEngage"
  }

  # Some features (such as code hits) are evaluated as present or not.
  # They do not have a further evaulation field such as labs.
  # When no additional featureEvaluation field provided, set equal to featureCol
  if(is.na(featureEvaluation)){
    featureEvaluation <- featureCol
  }

  # Subset feature data to feature of interest
  if(!is.na(featureCol)){
    featureDat <- featureDat[ featureDat[,featureCol] %in% feature &
                                !is.na(featureDat[,featureCol]), ]
  }

  # Subset index data to those with observed t0Col column
  indexT0    <- index[ !is.na(index[,t0Col]), ]


  #### Reduce a participant's observations to a single record closest to Baseline ####
  #### -------------------------------------------------------------------------------

  # Merge in t0 date
  indexRcdMerge <- merge(indexT0[,c(mergeby,t0Col)], featureDat, by=mergeby, all.x = TRUE)


  # Generic column called closestRecord from features date column
  indexRcdMerge[,"closestRecord"]    <- as.Date(indexRcdMerge[,featureDate])


  # Find difference (and abslute difference) between record and t0 + t0PlusDays
  indexRcdMerge[,"DiffBl_t0PlusDays"] <- indexRcdMerge[,"closestRecord"] -
                                          (indexRcdMerge[,t0Col] + t0PlusDays)

  # Find closest Rcd to t0 + t0PlusDays
  # Note: allows records which do not have an observed date.  Example, a cancer
  #       diagnosis but with unknown date is an exclusion criteria as it may have
  #       been prior to time zero.
  # Note2: With flatiron data, all labs / vitals have dates.  Some code hits do not have dates
  dropTooEarlyBl  <- which(indexRcdMerge[,"DiffBl_t0PlusDays"] < daysRangeLo)
  dropTooLateBl   <- which(indexRcdMerge[,"DiffBl_t0PlusDays"] > daysRangeHi)
  indexRcdMergeBl <- indexRcdMerge[-c(dropTooEarlyBl,dropTooLateBl),]


  # Also calculate absolute difference from t0Col + t0From Baseline
  # (The window for finding closest record may be before and after t0Col + t0PlusDays)
  indexRcdMergeBl[,"AbsDiffBl_t0PlusDays"] <- abs(indexRcdMergeBl[,"DiffBl_t0PlusDays"])


  # If multiple feature evaluations per day reduce to one record per date
  # General strategy
  # vitals: median of the day
  # labs: worst (least healthy) measure of day
  #  - max for all but gfr and albumin (use min for these two)
  #  - urine protein is categorized but in order of the numeric provided
  if(!is.na(dayReduction)){

    # New key merging primary key with date of obersved evaluation
    indexRcdMergeBl$pidDate <- paste0(indexRcdMergeBl[,mergeby], indexRcdMergeBl$closestRecord)

    # Single day reduction
    if(dayReduction=="median"){

      valReduced <- c(tapply(X     = indexRcdMergeBl[,featureEvaluation],
                             INDEX = indexRcdMergeBl[,"pidDate"],FUN = median))

    } else if (dayReduction=="max"){

      valReduced <- c(tapply(X     = indexRcdMergeBl[,featureEvaluation],
                             INDEX = indexRcdMergeBl[,"pidDate"],FUN = maxFun))

    } else if(dayReduction=="min"){

      valReduced <- c(tapply(X     = indexRcdMergeBl[,featureEvaluation],
                             INDEX = indexRcdMergeBl[,"pidDate"],FUN = minFun))

    } else {

      print(paste0("lab/vital [",featureName,"] not reduced by date and should be."))

    }


    # Merge in reduced, single day evaluations
    singleDaysEvals <- data.frame(pidDate=names(valReduced), valReduced)
    indexRcdMergeBl <- merge(indexRcdMergeBl, singleDaysEvals, by="pidDate")

    # New name of featureCol after reducing to single observation per day
    featureEvaluation <- "valReduced"

  }


  # Sort by:
  # 1) observed feature evaluation
  # 2) absolute difference from baseline
  # 3) difference from baseline (obs -2 days before baseline will be preffered over 2 days from baseline)
  # decreasing allows for finding the days closest to baseline or furthest from baseline
  indexRcdMergeBl <- indexRcdMergeBl[order(indexRcdMergeBl[,mergeby],
                                           is.na(indexRcdMergeBl[,featureEvaluation]),
                                           indexRcdMergeBl[,"AbsDiffBl_t0PlusDays"],
                                           indexRcdMergeBl[,"DiffBl_t0PlusDays"],
                                           decreasing = furthestFromT0),]

  # new object closestRcdBl of each participant's closest record
  closestRcdBl0  <- indexRcdMergeBl[!duplicated(indexRcdMergeBl[,mergeby]),]

  # Merge with index data to get one record per participant
  closestRcdBl <- merge(index[,c(mergeby,t0Col)],
                        closestRcdBl0[,c("patientid",featureEvaluation,featureDate,"closestRecord")],
                        by    = mergeby,
                        all.x = TRUE)



  # Rename variables and create indicators for comorbidities, rx fills, and outcome
  # Outpatient encounters named above as numOutpatientDischargeDates
  blPrefix <- "Bl_"

  if(featureType %in% c("dxHistory","rxHistory")) {

    blName <- paste0(blPrefix,featureName)

    # Diagnosis / Utilization but no date
    whichUNK <- which(!is.na(closestRcdBl[,featureDate]) &  is.na(closestRcdBl[,"closestRecord"]))
    # Diagnosis / Utilization in window
    whichYES <- which(!is.na(closestRcdBl[,featureDate]) & !is.na(closestRcdBl[,"closestRecord"]))

    closestRcdBl[whichUNK,blName]  <- -1
    closestRcdBl[whichYES,blName]  <-  1
    if(length(c(whichUNK,whichYES)) > 0){
      closestRcdBl[-c(whichUNK,whichYES),blName]  <-  0
    } else closestRcdBl[,blName]  <-  0

  } else if(featureType %in% c("vitalHistory","labHistory","ecogHistory")){

    blName                 <- paste0(blPrefix,featureName)
    closestRcdBl[,blName]  <- closestRcdBl[,featureEvaluation]

  } else if(featureType == "gap90check"){

    blName                 <- paste0("g90_",featureName)
    closestRcdBl[,blName]  <- closestRcdBl[,featureEvaluation]

  } else if(featureType == "pdl1"){

    blName                 <- paste0(blPrefix,featureName)
    closestRcdBl[,blName]  <- closestRcdBl[,featureEvaluation]
    "Tested and unknown"   -> closestRcdBl[closestRcdBl[,blName] %in%
                                             c("No interpretation given in report",
                                               "Results pending",
                                               "Unknown"),blName]
    "Not tested or missing" -> closestRcdBl[is.na(closestRcdBl[,blName]),blName]


  } else if(featureType == "follup"){

    blName                 <- paste0("last_",featureName)
    closestRcdBl[,blName]  <- closestRcdBl[,featureEvaluation]

  }


  # Rename date column
  closestRcdBl[,paste0(blName,"_Date")] <- closestRcdBl[,"closestRecord"]


  # Reorder by original rownames
  # Keep only feature value and feature date
  out <- closestRcdBl[order(as.numeric(rownames(closestRcdBl))),
                      c("patientid",blName,paste0(blName,"_Date"))]


  return(out)

}

# Example
#
# diagnoses <- list(list(name="Heart",      code=codeHeart),
#                   list(name="Kidney",     code=codeKidney),
#                   list(name="Neuropathy", code=codeNeuropathy),
#                   list(name="Hearing",    code=codeHearing),
#                   list(name="Carcinoma",  code=codeCancers))
#
# print("using diagnosisdate for previous hx of carcinoma")
# for(d in 1:length(diagnoses)){
#
#   # Look for any carcinoma diagnosis prior to *diagnosisdate* and without time limit
#   if(diagnoses[[d]]$name == "Carcinoma"){
#     dx_t0Col  <- "diagnosisdate"
#     daysLo <- NULL
#   } else {
#     dx_t0Col  <- t0Col
#     daysLo <- daysRangeLo
#   }
#
#   blToMerge <- findBlFollup(featureDat  = dxData, index = ptData, mergeby = "patientid",
#                             featureCol  = "diagnosiscodeAny", feature     = diagnoses[[d]]$code,
#                             featureDate = "diagnosisdateAny", featureName = diagnoses[[d]]$name,
#                             featureType = "dxHistory",
#                             t0Col       = dx_t0Col,
#                             daysRangeLo = daysLo,
#                             daysRangeHi = daysRangeHi)
#
#   ptData <- merge(ptData, blToMerge, by="patientid")
# }

