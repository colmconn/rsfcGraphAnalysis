#!/usr/bin/Rscript

rm(list=ls())

source("scoreMasc.r")

## Reads the seed file and does the (crude) equivalent of BASH variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character(), quiet=TRUE)
    table=gsub("$DATA", seeds.data.dir, table, fixed=TRUE)

    return (table)
}

## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed
getSeedName <- function(inSeedPath){
    name=basename(inSeedPath)
    if (grepl("\\.nii", name)) {
        return(gsub("\\.nii.*", "", name))
    } else if (grepl("\\+tlrc", name)) {
        return(gsub("\\+tlrc.*", "", name))
    } else {
        return (name)
    }
}

readSubjectOrderTable <- function (inFilename) {
    cat("*** Reading", subjectOrderFilename, "\n")
    subjectOrder=read.csv(subjectOrderFilename, header=T)

    return(subjectOrder)
}

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
}

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}

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


##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################


if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

## ncpus=1

scripts.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts")
data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=normalizePath(file.path(data.dir, "admin"))
config.data.dir=normalizePath(file.path(data.dir, "config"))
group.data.dir=file.path(data.dir, "Group.data")
data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")

seeds.data.dir=file.path(data.dir, "seeds")

demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)

groups=c("mddAndCtrl") ##,"ctrlOnly", "mddOnly")

rvVariable="WASI.Full.4"
## rvVariable="MASC.tscore"
seeds=readSeedsFile(file.path(config.data.dir, "juelich_amygdala_seeds_weights.txt"))

for (grouping in groups) {
    for (seed in seeds) {
        seedName=getSeedName(seed)
        
        cat("####################################################################################################\n")
        cat(sprintf("*** Creating covariate file for the %s seed of the %s grouping\n", seedName, grouping))
        
        ## this file stores the order of the subjects in each of the following BRIK files
        subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", grouping, seedName, "csv", sep="."))
        subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))
        
        if ( rvVariable == "MASC.tscore" ) {
            mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "MASC.total")])
        } else {
            mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", rvVariable)])
        }
        ## ensure that subject is s factor
        mgd$subject=as.factor(mgd$subject)
        mgd=droplevels(mgd)
        
        if ( rvVariable == "MASC.tscore" ) {
            mgd=fixDates(mgd)
            mgd=computeAge(mgd)
            
            mgd=computeMascScore(mgd)
        }
        
        ## print(mgd)

        ## when the covariate.table contains both groups, it (and mgd
        ## above) will be in the order of controls, then MDDs since
        ## the control group is named NCL and it sorts before MDD
        ## alphabetically. The order of subjects will be dictated by
        ## the order in which they are present in the subjectOrder
        ## file
        covariate.table=mgd[, c("subject", rvVariable)]
        covariate.table$subject=paste(as.character(covariate.table$subject), "A", sep="_")
        ## print(covariate.table)

        covariate.filename=file.path(group.data.dir, paste("covariates", grouping, seedName, "tab", sep="."))
        cat("*** Writing covariate file to", covariate.filename, "\n")
        write.table(covariate.table, covariate.filename, quote=FALSE, col.names=TRUE, row.names=FALSE)
    }
}

    
