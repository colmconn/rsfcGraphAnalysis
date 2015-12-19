rm(list=ls())
graphics.off()

####################################################################################################
### Functions
####################################################################################################

readStatsTable <- function (inFilename) {

    cat("*** Reading" , inFilename, "\n")
    statsTable=read.table(inFilename, header=T, sep="")
    ## dump the first column as it's only the file name
    statsTable=statsTable[, -1]
    return(statsTable)
}

readSubjectOrderTable <- function (inFilename) {
    cat("*** Reading", inFilename, "\n")
    subjectOrder=read.csv(inFilename, header=T)

    return(subjectOrder)
}

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
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

####################################################################################################
###
####################################################################################################

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")
admin.data.dir=file.path(data.dir, "admin")

task="restingstate"
usedFwhm="4.2"
groups="mddOnly"
seedName="R_whole_amygdala.3mm"
rvVariable="CDRS.t.score.scaled.diff"

infix=sprintf("regression.fwhm%s.%s.%s.%s.and.%s", usedFwhm, task, groups, seedName, rvVariable)

group.data.dir=file.path(data.dir, paste("Group.data", rvVariable, sep="."))
group.results.dir=file.path(data.dir, paste("Group.results", rvVariable, sep="."))

roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))

## this file stores the order of the subjects in each of the following BRIK files
subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", groups, seedName, "csv", sep="."))
subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))

cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))


roistats=readStatsTable(roistats.filename)
roistats$Sub.brick=NULL

## this file stores the order of the subjects in each of the following BRIK files
subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", groups, seedName, "csv", sep="."))
subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))

cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

mgd=cbind(subjectOrder, roistats)


demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ""))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))


selectedColumns=c("Grp", "Gender", "DOB", "MRI",  "CDRS.tscore", "CoRum", "PSWQ", "RSQ.fixed")

mgd=cbind(mgd, demographics[match(mgd$subject, demographics$ID), selectedColumns])

mgd=computeAge(fixDates(mgd))


parahippocampal.rsq.cor=with(mgd, cor.test(Mean_3, RSQ.fixed, method="spearman", na.action="na.omit"))
cat("*** Spearman's correlation between L Parahippocampal gyrus and RSQ\n")
print(parahippocampal.rsq.cor)
cat("*** There are", sum(is.na(mgd$RSQ.fixed)), "NAs in the RSQ column\n")


parahippocampal.pswq.cor=with(mgd, cor.test(Mean_3, PSWQ, method="spearman", na.action="na.omit"))
cat("*** Spearman's correlation between L Parahippocampal gyrus and PSWQ\n")
print(parahippocampal.pswq.cor)
cat("*** There are", sum(is.na(mgd$PSWQ)), "NAs in the PSWQ column\n")

cdrsr.rsq.cor=with(mgd, cor.test(CDRS.tscore,  RSQ.fixed, method="spearman", na.action="na.omit"))
cat("*** Spearman's correlation between CDRS-R and RSQ\n")
print(cdrsr.rsq.cor)
