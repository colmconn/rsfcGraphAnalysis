rm(list=ls())
graphics.off()

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

admin.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/admin"))
config.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/config"))
group.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.data"))
group.results.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.results"))
scriptsDir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))

filterOutStimulusColumn <- function (inDataFrame) {
    if (dim(inDataFrame)[2] > 1) {
        cat("*** Filtering the data frame because it has more then 1 column\n")
        df=data.frame("subject"=unique(inDataFrame[, 1]))
    }
    return (df)
}

mdd.subjectOrder.filename=file.path(group.data.dir, "subjectOrder.mddOnly.L_BLA.weight.3mm.csv")

ksads.filename=file.path(admin.data.dir, "KSADSDiagnoses_tiffany.csv")

cat("Reading list of MDDs subject in the bucket file: ", mdd.subjectOrder.filename, "\n")
mdd.subjectList=read.csv(mdd.subjectOrder.filename, header=TRUE)
colnames(mdd.subjectList)=c("subject")
cat(sprintf("Read %d subjects from the MDD subject bucket list\n", dim(mdd.subjectList)[1]))

cat("Reading KSADS data from: ", ksads.filename, "\n")
ksads=read.csv(ksads.filename, header=TRUE, skip=1)
cat(sprintf("Read %d subjects from the KSADS file\n", dim(ksads)[1]))

mgd=cbind(mdd.subjectList$subject, ksads[match(mdd.subjectList$subject, ksads$Subject.Number), -c(1:8)])
colnames(mgd)=c("subject", colnames(mgd)[-1])
missing=mgd[!complete.cases(mgd), ]
mgd=mgd[complete.cases(mgd), ]

cat("The following subjects are missing KSADS data:\n")
cat(as.vector(missing$subject), "\n")

for ( cName in colnames(mgd)[-1] ) {
    if (is.factor(mgd[ , cName])) {
        cat("*** Converting", cName, "to vector\n")
        mgd[, cName] = as.integer(mgd[, cName])
    }
    if (cName == "COMO.") {
        cat(sprintf("%s: %d (%s)\n", cName,              sum(mgd[, cName], na.rm=TRUE), paste(as.character(mgd[mgd[, cName] == 1, "subject"]), collapse=" ")))
        cat(sprintf("%s: %d (%s)\n", paste ("NO", cName), sum(1 - mgd[, cName], na.rm=TRUE), paste(as.character(mgd[mgd[, cName] == 1, "subject"]), collapse=" ")))
    } else {
        cat(sprintf("%s: %d (%s)\n", cName, sum(mgd[, cName], na.rm=TRUE), paste(as.character(mgd[mgd[, cName] == 1, "subject"]), collapse=" ")))
    }
}
