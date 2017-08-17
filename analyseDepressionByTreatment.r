rm(list=ls())
##graphics.off()

library(plyr)

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}


if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=normalizePath(file.path(data.dir, "admin"))
config.data.dir=normalizePath(file.path(data.dir, "config"))
seeds.data.dir=normalizePath(file.path(data.dir, "seeds"))


## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
## demographics=readCsvFile(demographicsFilename)
## demographics=rename(demographics, c("Grp"="Group"))

## medsAndTreatment.filename=file.path(admin.data.dir, "medsAndTreatment.csv")
medsAndTreatment.filename=file.path(admin.data.dir, "new_medsAndTreatment.csv")
medsAndTreatment=readCsvFile(medsAndTreatment.filename, "SubjID")

cdrsr.filename=file.path(admin.data.dir, sprintf("new.CDRS.t.score.timepoint.a.and.c.score.csv"))
cdrsr.scores=readCsvFile(cdrsr.filename, "ID")

if ( any(cdrsr.scores$ID=="300") ) {
    cat ("*** Found subject 169/300 as number 300. Fixing.\n")
    cdrsr.scores$ID[cdrsr.scores$ID=="300"]= "169/300"
    cdrsr.scores$ID=as.factor(cdrsr.scores$ID)
}
cdrsr.scores$ID=as.factor(cdrsr.scores$ID)
cdrsr.scores=rename(cdrsr.scores, c(Grp="Group"))

ss.cdrsr.scores=subset(cdrsr.scores, Group=="MDD" & timepoint=="C")

ss.medsAndTreatment=subset(medsAndTreatment, Timepoint=="C", select=c("SubjectNumber", "SubjID", "Timepoint", "Group", "Rx.Summary", "Tx.Summary"))
## print(ss.medsAndTreatment)

## now try to add the meds and treatment info for each subject

mgd=cbind(ss.cdrsr.scores,
    ss.medsAndTreatment[match(ss.cdrsr.scores$ID, ss.medsAndTreatment$SubjID), c("Rx.Summary", "Tx.Summary")])
rownames(mgd)=NULL

mgd$Tx=as.factor(ifelse(mgd$Tx.Summary=="No Treatment", "No Treatment", "Treatment"))
print(summary(mgd$Tx))
print(summary(mgd$Tx.Summary))
print(summary(mgd$Rx.Summary))
print(summary(aov(CDRS.t.score ~ Tx,              data=mgd)))
print(summary(aov(CDRS.t.score ~ Tx + Rx.Summary, data=mgd)))
print(summary(aov(CDRS.t.score ~ Tx * Rx.Summary, data=mgd)))
