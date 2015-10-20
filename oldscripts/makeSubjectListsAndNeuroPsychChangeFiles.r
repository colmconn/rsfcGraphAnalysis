rm(list=ls())
graphics.off()

source("scoreMasc.r")

####################################################################################################
### START OF FUNCTIONS
####################################################################################################

fixDates <- function (inData) {
    ## this complicated looking regexp stuff cleans up the years in DOB
    ## and MRI with 4 digits to be just 2 digits
    ## month day year
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$MRI)
    
    ## now convert to year/month/day
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$MRI)
    
    inData$DOB=as.Date(inData$DOB, "%y/%m/%d")
    inData$MRI=as.Date(inData$MRI, "%y/%m/%d")

    return(inData)
}

computeAge <- function(inData) {

    age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
    age.in.weeks=as.numeric(age.in.weeks)

    inData$age.in.years=age.in.weeks/52

    return(inData)
}

computeMascScore <- function (inData) {
    inData.dim=dim(inData)
    for (r in seq(1, inData.dim[1]) ) {
        ## cat("##################################################\n")
        subjectNumber=inData[r, "subject"]
        gender=inData[r, "Gender"]
        age=round(inData[r, "age.in.years"], 0)
        old.masc.tscore=inData[r, "MASC.tscore"]

        new.masc.tscore=scoreMasc(gender, age, inData[r, "MASC.total"])
        if (is.na(new.masc.tscore) ) {
            warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, inData[r, "MASC.total"]))
        }
        
        inData[r, "MASC.tscore"]=new.masc.tscore

        ## cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f new MASC tscore=%0.0f\n",
        ##             r, subjectNumber, gender, age, inData[r, "MASC.total"], old.masc.tscore, new.masc.tscore))

    }
    return (inData)
}

####################################################################################################
### END OF FUNCTIONS
####################################################################################################

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't determine the number of CPUs in this architecture defaulting to", ncpus, "\n"))
    stop(paste("Sorry can't set data directories for this computer\n"))
}

data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/admin"))
config.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/config"))

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "", " NO SUBJECT "))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

############################################################
cdrsrFilename=file.path(admin.data.dir, "CDRSR.csv")
cat("*** Reading", cdrsrFilename, "\n")
cdrsr=read.csv(cdrsrFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "", " NO SUBJECT ", "Test_Subject"))
cat(sprintf("*** Read CDRSR data for %s unique subjects\n",  length(unique(cdrsr$SubjID))))

## dump non compete cases
cdrsr=cdrsr[complete.cases(cdrsr), ]
## dump all but time points A and C
cdrsr=cdrsr[cdrsr$timepoint %in% c("A", "C"), ]

cdrsr$id=with(cdrsr, paste(SubjID, timepoint, sep="_"))

## now find only those subjects with time points A and C
bothTimePoints=vector(mode="character", length=length(cdrsr$SubjectID))
nn=1
for (subject in cdrsr$SubjID) {
##    cat(subject, "\n")
    if ( paste(subject, "A", sep="_") %in% cdrsr$id &&
         paste(subject, "C", sep="_") %in% cdrsr$id) {
        bothTimePoints[nn]=subject
        nn=nn+1
    }
}
bothTimePoints=unique(bothTimePoints)
cat(bothTimePoints, "\n")
cdrsr=cdrsr[which (cdrsr$SubjID %in% bothTimePoints) , ] 
cdrsr$id=NULL
cdrsr$CDRSR.raw=NULL

cdrsr.change=reshape(cdrsr,
    timevar="timepoint",
    idvar=c("SubjID"),
    direction="wide")
cdrsr.change$CDRSR.diff=with(cdrsr.change, CDRSR.tscore.C - CDRSR.tscore.A)
cdrsr.change$CDRSR.percentDiff= with(cdrsr.change, CDRSR.diff / CDRSR.tscore.A)*100


############################################################

mascAndBdiAndCgasAndRrsFilename=file.path(admin.data.dir, "MascAndBdiAndCgasAndRrs.csv")
cat("*** Reading", mascAndBdiAndCgasAndRrsFilename, "\n")
mascAndBdiAndCgasAndRrs=read.csv(mascAndBdiAndCgasAndRrsFilename, header=T,
    na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "", " NO SUBJECT ", "Test_Subject"))
cat(sprintf("*** Read mascAndBdiAndCgasAndRrs data for %s unique subjects\n",  length(unique(mascAndBdiAndCgasAndRrs$SubjID))))


radsAndCdiFilename=file.path(admin.data.dir, "RadsAndCdi.csv")
cat("*** Reading", radsAndCdiFilename, "\n")
radsAndCdi=read.csv(radsAndCdiFilename, header=T,
    na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "", " NO SUBJECT ", "Test_Subject"))
cat(sprintf("*** Read radsAndCdi data for %s unique subjects\n",  length(unique(radsAndCdi$SubjID))))


## dump non compete cases
## mascAndBdiAndCgasAndRrs=mascAndBdiAndCgasAndRrs[complete.cases(mascAndBdiAndCgasAndRrs), ]

## dump all but time points A and C
mascAndBdiAndCgasAndRrs=mascAndBdiAndCgasAndRrs[mascAndBdiAndCgasAndRrs$timepoint %in% c("A", "C"), ]

radsAndCdi=radsAndCdi[radsAndCdi$timepoint %in% c("A", "C"), ]

##mascAndBdiAndCgasAndRrs$id=with(mascAndBdiAndCgasAndRrs, paste(SubjID, timepoint, sep="_"))

## now find only those subjects with time points A and C
## bothTimePoints=vector(mode="character", length=length(mascAndBdiAndCgasAndRrs$ID))
## nn=1
## for (subject in mascAndBdiAndCgasAndRrs$SubjNum) {
## ##    cat(subject, "\n")
##     if ( paste(subject, "A", sep="_") %in% mascAndBdiAndCgasAndRrs$ID &&
##          paste(subject, "C", sep="_") %in% mascAndBdiAndCgasAndRrs$ID) {
##         bothTimePoints[nn]=subject
##         nn=nn+1
##     }
## }

mascAndBdiAndCgasAndRrs=fixDates(mascAndBdiAndCgasAndRrs)
mascAndBdiAndCgasAndRrs=computeAge(mascAndBdiAndCgasAndRrs)

## add gender as we need it to compute MASC
mascAndBdiAndCgasAndRrs$Gender=demographics[match(mascAndBdiAndCgasAndRrs$SubjNum, demographics$ID), c("Gender")]
radsAndCdi$Gender=demographics[match(radsAndCdi$SubjNum, demographics$ID), c("Gender")]

## keep only male and female
mascAndBdiAndCgasAndRrs=mascAndBdiAndCgasAndRrs[which(mascAndBdiAndCgasAndRrs$Gender %in% c("F", "M")), ]
radsAndCdi=radsAndCdi[which(radsAndCdi$Gender %in% c("F", "M")), ]

mascAndBdiAndCgasAndRrs=computeMascScore(mascAndBdiAndCgasAndRrs)

## bothTimePoints=unique(bothTimePoints)
## cat(bothTimePoints, "\n")

## mascAndBdiAndCgasAndRrs=mascAndBdiAndCgasAndRrs[which (mascAndBdiAndCgasAndRrs$SubjNum %in% bothTimePoints) , ]
## eliminate subjects with NO BDI score
## mascAndBdiAndCgasAndRrs=mascAndBdiAndCgasAndRrs[ mascAndBdiAndCgasAndRrs$BDI >= 0 , ]

## now filter out unwanted columns
##mascAndBdiAndCgasAndRrs = mascAndBdiAndCgasAndRrs[, c("ID", "SubjNum", "timepoint", "BDI", "MASC.tscore")]
## stop()
## mascAndBdiAndCgasAndRrs.change=reshape(mascAndBdiAndCgasAndRrs,
##     timevar="timepoint",
##     idvar=c("ID", "SubjNum"),
##     direction="wide")
## mascAndBdiAndCgasAndRrs.change$mascAndBdiAndCgasAndRrs.diff=with(mascAndBdiAndCgasAndRrs.change, mascAndBdiAndCgasAndRrs.tscore.C - mascAndBdiAndCgasAndRrs.tscore.A)
## mascAndBdiAndCgasAndRrs.change$mascAndBdiAndCgasAndRrs.percentDiff= with(mascAndBdiAndCgasAndRrs.change, mascAndBdiAndCgasAndRrs.diff / mascAndBdiAndCgasAndRrs.tscore.A)*100


## This variable controls whether the list of subjects, irrespective
## of whether they have usable MRI data or not, should be generated
canonicalTimeAList=FALSE

## use this switch set to TRUE if you want to make subject lists for
## analysing MRI data at time A and C
canonicalTimeAandCList=FALSE

createCdrsrSubjectLists=FALSE
############################################################

##stop()

## now create a data frame with subject IDs and a column ("use") which
## indicates whether the subject can be included in the analysis or not.
##
## of the rsfc directory exist and rsfcPreprocessed/00_DO_NOT_ANALYSE
## does not exist use will be TRUE. It will be false otherwise. The
## applies only in the else branch fo the if condition below.

## subjects=$( cd $DATA ;  ls -d *_[A] )
if (canonicalTimeAandCList) {
    initial.directory.list=dir(data.dir, pattern="[0-9]{3}_[AC]$", no.. = TRUE)
} else {
    initial.directory.list=dir(data.dir, pattern="[0-9]{3}_A$", no.. = TRUE)
}

excessiveMotionThresholdPercentage=30
cat("****************************************************************************************************\n")
cat("*** Using", excessiveMotionThresholdPercentage, "% as motion cutoff threshold\n")
cat("****************************************************************************************************\n")

use.initial.directories=
    as.data.frame(
        cbind(initial.directory.list,
              gsub("[0-9]{3}_", "", initial.directory.list),
              ## rep("A", length(initial.directory.list)),
              if ( canonicalTimeAList ) {
                  sapply(initial.directory.list, function (x) {
                      ff=file.path(data.dir, x, "rsfcPreprocessed", sprintf("00_DO_NOT_ANALYSE_%s_%dpercent.txt", x, excessiveMotionThresholdPercentage))
                      if ( file.exists(ff) ) {
                          return(TRUE)
                      } else {
                          return(TRUE)
                      }
                  })
              } else {
                  ## the code below actually sets the tru or false
                  ## correctly based on the exisitance of not of the
                  ## DO_NOT_ANALYSE files
                  sapply(initial.directory.list, function (x) {
                      ff=file.path(data.dir, x, "rsfcPreprocessed", sprintf("00_DO_NOT_ANALYSE_%s_%dpercent.txt", x, excessiveMotionThresholdPercentage))
                      dd=file.path(data.dir, x, "rsfc")                      
                      if ( file.exists(dd) && ! file.exists(ff) ) {
                          return(TRUE)
                      } else {
                          return(FALSE)
                      }
                  })
               }
              ))

colnames(use.initial.directories)=c("subject", "timepoint", "use")

use.initial.directories$subject=gsub("_[AC]", "", initial.directory.list, fixed=FALSE)
if ( any(grep("300", use.initial.directories$subject)) ) {
    use.initial.directories$subject[which (use.initial.directories$subject=="300")]="169/300"
}

use.initial.directories=cbind(use.initial.directories, demographics[match(use.initial.directories$subject, demographics$ID), c("Grp", "Gender")])
## 425 is missing group and gender so add it manually
use.initial.directories[use.initial.directories$subject=="425", c("Grp", "Gender")]=c("NCL", "F")

cat("*** Number of subjects before removing those on meds ", dim(use.initial.directories)[1], "\n")

## Now remove subjects that are NOT medication-naive (Prozac, Zoloft, Celexa, Citalopram, Klondpin, Seroquil, Cymbalta)
use.initial.directories <- use.initial.directories[! use.initial.directories$subject %in% c("130", "132", "318", "319", "320", "322", "323", "324", "325", "329", "345", "376", "384"), ]

## Now remove subjects that were on ADHD medication (Focolin, Stratera)
use.initial.directories <- use.initial.directories[! use.initial.directories$subject %in% c("333", "349"), ]

## Now remove ctrls with CDRSR > 54
use.initial.directories <- use.initial.directories[! use.initial.directories$subject %in% c("119", "346", "370"), ]

##now ensure that we have only NCL, MDD and M, and F subjects
use.initial.directories=subset(use.initial.directories, Grp %in% c("NCL", "MDD") & Gender %in% c("M", "F"))

use.initial.directories$id=with(use.initial.directories, paste(subject, timepoint, sep="_"))
rownames(use.initial.directories)=NULL
use.initial.directories=droplevels(use.initial.directories)

cat("*** Number of subjects after removing those on meds ", dim(use.initial.directories)[1], "\n")

cat("*** Number of subjects going into the analysis by group and gender\n")
print(with(use.initial.directories, addmargins(table(Grp, Gender))))

if ( canonicalTimeAList ) {

    ncl.list=subset(use.initial.directories, Grp=="NCL" & use==TRUE & timepoint=="A", select=id)[, 1]
    mdd.list=subset(use.initial.directories, Grp=="MDD" & use==TRUE & timepoint=="A", select=id)[, 1]

    
    ## Write a canonical list of subjects woul could be included in
    ## the analysis. This should be taken as a starting point that
    ## will be whittled down depending on whether the subjects are too
    ## contaminated with motion or are missing behavioral data at one
    ## or more timepoints (e.g., A and C)
    ncl.list.filename=file.path(config.data.dir, "ncl.subjectList.txt")
    mdd.list.filename=file.path(config.data.dir, "mdd.subjectList.txt")

    cat("*** Writing NCL list to " , ncl.list.filename, "\n")
    write.table(ncl.list, ncl.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    cat("*** Writing MDD list to " , mdd.list.filename, "\n")
    write.table(mdd.list, mdd.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
}

if (canonicalTimeAandCList) {

    bothTimePoints=vector(mode="character", length=length(use.initial.directories$ID))
    nn=1
    for (subject in use.initial.directories$subject) {
        ## cat(subject, "\n")
        if (paste(subject, "A", sep="_") %in% use.initial.directories$id &&
            paste(subject, "C", sep="_") %in% use.initial.directories$id) {
            bothTimePoints[nn]=subject
            nn=nn+1
        }
    } ## end of for (subject in ss$SubjNum) {

    cat("*** List of subjects with both time A and C timepoints:\n")
    cat(paste(bothTimePoints, collapse=" "), "\n")

    ss=subset( use.initial.directories, subject %in% bothTimePoints)

    ## further filter by only those subjects that have usable MR data
    ## at time A and C

    bothTimePointsUsable=vector(mode="character", length=dim(ss)[1]/2)
    nn=1
    for (subj in unique(ss$subject) ) {
        if (ss[ss$subject == subj & ss$timepoint=="A", "use"] == TRUE && ss[ss$subject == subj & ss$timepoint=="C", "use"] == TRUE) {
            ## print(ss[ss$subject==subj, ])
            bothTimePointsUsable[nn]=subj
            nn=nn+1
        }
    }

    cat("*** List of subjects with usable timepoint A and C MRI data\n")
    cat(paste(bothTimePointsUsable, collapse=" "), "\n")

    ss=subset(ss, subject %in% bothTimePointsUsable)
    
    ncl.list=subset(ss, Grp=="NCL", select=id)[, 1]
    mdd.list=subset(ss, Grp=="MDD", select=id)[, 1]

    cat ("*** NCL List\n")
    print(ncl.list)
    cat ("*** MDD List\n")    
    print(mdd.list)

    cat("*** Time x Grp x usability table\n")
    print(addmargins(with(ss, table(timepoint, Grp, use))))

    ncl.list.filename=file.path(config.data.dir, "ncl.subjectList.AandC.txt")
    mdd.list.filename=file.path(config.data.dir, "mdd.subjectList.AandC.txt")
    
    cat("*** Writing NCL list to " , ncl.list.filename, "\n")
    write.table(ncl.list, ncl.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    cat("*** Writing MDD list to " , mdd.list.filename, "\n")
    write.table(mdd.list, mdd.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
}


####################################################################################################

## Now generate the lists of subjects for CDRSR analysis. This is done
## seperately because the CDRSR scores are loaded in with one subject
## per row and the CDRSR timepoints A and C in separate columns

bothTimePointSubjects= use.initial.directories[ which (use.initial.directories$subject %in% cdrsr.change$SubjID), ]
cat("*** Number of subjects who have CDRSR at time A and C going into the analysis by group and gender before eliminating those contaminated by motion\n")
print(with(bothTimePointSubjects, addmargins(table(Grp, Gender))))
cat("*** Number of subjects who have CDRSR at time A and C going into the analysis by group and MRI usability\n")
print(with(bothTimePointSubjects, addmargins(table(use, Grp))))

ncl.list=subset(bothTimePointSubjects, Grp=="NCL" & use==TRUE, select=id)[, 1]
mdd.list=subset(bothTimePointSubjects, Grp=="MDD" & use==TRUE, select=id)[, 1]

cdrsr.change.score.filename=file.path(admin.data.dir, "cdrsr.change.csv")
cat("*** Saving CDRSR change scores to", cdrsr.change.score.filename, "\n")
write.table(cdrsr.change, cdrsr.change.score.filename, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=",")
    
if ( createCdrsrSubjectLists ) { 
    ncl.list.filename=file.path(config.data.dir, "ncl.subjectList.with.Cdrsr.AandC.txt")
    mdd.list.filename=file.path(config.data.dir, "mdd.subjectList.with.Cdrsr.AandC.txt")

    cat("*** Writing NCL list to " , ncl.list.filename, "\n")
    write.table(ncl.list, ncl.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    cat("*** Writing MDD list to " , mdd.list.filename, "\n")
    write.table(mdd.list, mdd.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
}


####################################################################################################

makeSubjectList <- function (df, columns) {

    for (col in columns ) {

        cat("####################################################################################################\n")
        cat("### Creating subject lists for", col, "\n")
        
        ss = df[, c("ID", "SubjNum", "timepoint", col)]

        bothTimePoints=vector(mode="character", length=length(ss$ID))
        for (subject in ss$SubjNum) {
            ## cat(subject, "\n")
            if (paste(subject, "A", sep="_") %in% ss$ID &&
                paste(subject, "C", sep="_") %in% ss$ID) {
                bothTimePoints[nn]=subject
                nn=nn+1
            }
        }

        ## keep only thouse subjects who have the behavioral measure at
        ## time points A and C
        cat("*** Removing", sum(with(ss, ! (SubjNum %in% bothTimePoints ))), "subjects who dont have", col, "data for both time points\n")

        ss=subset(ss, SubjNum %in% bothTimePoints)

        nLessThanZeroCases=sum(ss[, col] < 0, na.rm=TRUE)
        cat("*** Removing", nLessThanZeroCases, "subjects with", col, "data less than 0\n")
        ss=subset(ss, col >= 0)

        ## eliminate the ID column before we try to run reshape on it. The
        ## ID column is unnecessary
        ss$ID = NULL
        ss.change=reshape(ss,
            timevar="timepoint",
            idvar=c("SubjNum"),
            direction="wide")
        ss.change$ss.diff       = ss.change[, paste(col, "C", sep=".")] - ss.change[, paste(col, "A", sep=".")] ## with(ss.change, ss.tscore.C - ss.tscore.A)
        ss.change$ss.percentDiff= (ss.change$ss.diff     / ss.change[, paste(col, "A", sep=".")]) * 100

        cat("*** Table of case completeness for", col, "\n")
        print(table(complete.cases(ss.change)))
        cat("*** Removing", sum(!complete.cases(ss.change)), "with missing", col, "data\n")
        ss.change=ss.change[complete.cases(ss.change), ]

        colnames(ss.change) = c(colnames(ss.change)[1:(length(colnames(ss.change))-2)], paste(col, c("diff", "percentDiff"), sep="."))

        ## this file contains measures from both time point A and C
        ## irrespective of whether there is usable MRI data or not
        csv.filename=file.path(admin.data.dir, paste(col, "change", "csv", sep="."))
        cat("*** Saving", col, "change scores to", csv.filename, "\n")
        write.table(ss.change, csv.filename, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=",")

        ss.change$Grp = demographics[match (ss.change$SubjNum, demographics$ID), "Grp"]

        ## now compute the intersection of the subjects with both time
        ## point col (masc or BDI) data and those subjects who have usable
        ## MRI data (i.e., that data that was not contaminated with too
        ## much motion)

        bothTimePointSubjects= use.initial.directories[ which (use.initial.directories$subject %in% ss.change$SubjNum), ]

        cat("*** Number of subjects who have", col, "at time A and C going into the analysis by group and MRI usability\n")
        print(with(bothTimePointSubjects, addmargins(table(use, Grp))))
        
        ## these lists contain the lists of NCL and MDD, respectively, who
        ## have usable MRI data at time point A and time point A and C for
        ## the variable of interest
        ncl.list=subset(bothTimePointSubjects, Grp=="NCL" & use==TRUE, select=id)[, 1]
        mdd.list=subset(bothTimePointSubjects, Grp=="MDD" & use==TRUE, select=id)[, 1]

        ncl.list.filename=file.path(config.data.dir, paste("ncl.subjectList.with", col, "AandC.txt", sep="."))
        mdd.list.filename=file.path(config.data.dir, paste("mdd.subjectList.with", col, "AandC.txt", sep="."))
        
        cat("*** Writing NCL list to" , ncl.list.filename, "\n")
        write.table(ncl.list, ncl.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
        
        cat("*** Writing MDD list to" , mdd.list.filename, "\n")
        write.table(mdd.list, mdd.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
        
    }
}

####################################################################################################
### RADS and CDI
####################################################################################################

##makeSubjectList(mascAndBdiAndCgasAndRrs, c("MASC.tscore", "BDI", "CGAS"))
makeSubjectList(mascAndBdiAndCgasAndRrs, c("RRS"))
##makeSubjectList(radsAndCdi, c("RADS.Total.Tscore", "CDI"))
