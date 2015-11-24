rm(list=ls())
graphics.off()

## for the some package
library(car)

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

wasiFilename=file.path(admin.data.dir, "WASI.csv")
wasi=readCsvFile(wasiFilename, "SubID")


e.norm.filename=file.path(group.data.dir, "restingstate.mddAndCtrl.motion.enorm.csv")
enorm=fixSubjectOrderTable(readCsvFile(e.norm.filename, inSubjectColumnName="subject"))

## seeds=readSeedsFile(file.path(config.data.dir, "juelich_amygdala_seeds_weights.txt"))
seeds=readSeedsFile(file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt"))

grouping="mddAndCtrl"

for (seed in seeds) {
    seedName=getSeedName(seed)
    
    cat("####################################################################################################\n")
    cat(sprintf("*** Creating covariate file for the %s seed of the %s grouping\n", seedName, grouping))
    
    ## this file stores the order of the subjects in each of the following BRIK files
    subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", grouping, seedName, "csv", sep="."))
    subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))

    bucketListFilename=file.path(group.data.dir, paste("bucketList", grouping, seedName, "txt", sep="."))
    cat("*** Reading", bucketListFilename, "\n")
    bucketList=read.table(bucketListFilename, header=FALSE, as.is=TRUE)

    mgd=cbind(subjectOrder, bucketList,
        demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI")],
        wasi        [match(subjectOrder$subject, wasi$SubID),      c("Full")],
        enorm       [match(subjectOrder$subject, enorm$subject),   c("enorm")])
     ## wasi        [match(subjectOrder$subject, wasi$SubID),      c("Verbal", "Performance", "Full")])    

    if (any(mgd$subject=="378"))
        mgd[mgd$subject=="378", "Gender"]="F"

    mgd$subject=paste(mgd$subject, "A", sep="_")
    mgd=droplevels(mgd)        
    colnames(mgd)=c("Subj", "InputFile", "Group", "Gender", "DOB", "MRI", "Full", "enorm")
    
    mgd=fixDates(mgd)
    mgd=computeAge(mgd)
    colnames(mgd)=c("Subj", "InputFile", "Group", "Gender", "DOB", "MRI", "Full", "enorm", "age")    

    mgd[, c("age", "Full", "enorm")] = scale(mgd[, c("age", "Full", "enorm")], center=TRUE, scale=FALSE)

    ## don't scale the Full IQ as this differs between groups
    ## mgd[, c("age")] = scale(mgd[, c("age")], center=TRUE, scale=FALSE)

    ## tapply(mgd$Full, list(mgd$Group), mean) # get the means of the
    ## groups then repeat each as manu times as there members in each
    ## group and provide the resultant vector of repeated means to
    ## scale to subtract from each original value. The results is Full
    ## IQ centered within group
    ## mgd[, c("Full")] = mgd[, c("Full")] - rep(tapply(mgd$Full, list(mgd$Group), mean), table(mgd$Group))

    print(some(mgd))
    
    mgd$InputFile=sub("z-score", "z-score.masked", mgd$InputFile, fixed=TRUE)

    ## reorder the columns to suite 3dMVM
    ss=mgd[, c("Subj", "Group", "Gender", "Full", "InputFile")]
    dataTableFilename=file.path(group.data.dir, paste("dataTable", grouping, "group", "and", "gender", seedName, "txt", sep="."))
    cat("*** Writing data table to", dataTableFilename, "\n")
    write.table(format(ss, digits=4), file=dataTableFilename, quote=FALSE, col.names=TRUE, row.names=FALSE, eol=" \\\n")

    ## reorder the columns to suite 3dMVM
    ss=mgd[, c("Subj", "Group", "age", "Full", "InputFile")]
    dataTableFilename=file.path(group.data.dir, paste("dataTable", grouping, "group", "and", "age", seedName, "txt", sep="."))
    cat("*** Writing data table to", dataTableFilename, "\n")
    write.table(format(ss, digits=4), file=dataTableFilename, quote=FALSE, col.names=TRUE, row.names=FALSE, eol=" \\\n")

}
