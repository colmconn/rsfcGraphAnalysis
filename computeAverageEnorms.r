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
data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
group.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.data"))
group.results.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.results"))
scriptsDir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))

ctrl.subjectOrder.filename=file.path(config.data.dir, "flex.rsfc.ctrlOnly.csv")
mdd.subjectOrder.filename=file.path(config.data.dir,  "flex.rsfc.mddOnly.csv")

cat("Reading list of CONTROLs subject in the bucket file: ", ctrl.subjectOrder.filename, "\n")
ctrl.subjectList=read.csv(ctrl.subjectOrder.filename, header=FALSE)
colnames(ctrl.subjectList)=c("subject")
cat(sprintf("Read %d subjects from the CONTROL subject bucket list\n", dim(ctrl.subjectList)[1]))

cat("Reading list of MDDs subject in the bucket file: ", mdd.subjectOrder.filename, "\n")
mdd.subjectList=read.csv(mdd.subjectOrder.filename, header=FALSE)
colnames(mdd.subjectList)=c("subject")
cat(sprintf("Read %d subjects from the MDD subject bucket list\n", dim(mdd.subjectList)[1]))

subjectOrder=data.frame("subject"=c(as.vector(mdd.subjectList$subject), as.vector(ctrl.subjectList$subject)))
colnames(subjectOrder)=c("subject")
subjectOrder$subject=as.vector(subjectOrder$subject)

## demographicsFilename=file.path(admin.data.dir, "data_entry_current_021113.csv")
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))
demographics$Grp=droplevels(demographics$Grp)

##stop()

allSubjects=data.frame(
    "subject" = subjectOrder,
    "Group"    = demographics[match(gsub("_A[0-9]?", "", subjectOrder$subject, fixed=FALSE), demographics$ID), "Grp"]
    )
colnames(allSubjects)=c("subject", "Group")

nSubjects=dim(allSubjects)[1]

e.norm.matrix=matrix(NA, nrow=nSubjects*2, ncol=1)
subjectCount=1
task.dirs=c("rsfcPreprocessed", "estopPreprocessed")
pb <- txtProgressBar(min = 0, max = nSubjects*2, style = 3)

for ( task.dir in task.dirs) {
    for (subject in allSubjects$subject) {
        if (subject == "169/300_A") {
            subject="300_A"
        }
        if (task.dir == "rsfcPreprocessed" )
            task="funcon"
        else
            task="ESTOP"
        preprocessed.data.dir=file.path(data.dir, subject, task.dir, "tmp")

        owd=getwd()
        
        e.norm.command="1d_tool.py"
        e.norm.command.arguments=sprintf("-infile %s%s.zp.float.despike_tsh_vr_motion.tcat.1D -collapse_cols euclidean_norm -write e.norm.1D -overwrite", subject, task)

        setwd(preprocessed.data.dir)
        system2(e.norm.command, e.norm.command.arguments)

        e.norm.avg=mean(scan(file="e.norm.1D", quiet=TRUE))
        
        ## print(e.norm.avg)
        setwd(owd)
        e.norm.matrix[subjectCount, 1]=e.norm.avg
        setTxtProgressBar(pb, subjectCount)
        subjectCount=subjectCount+1
    }
}
cat("\n")

enorm.df=data.frame(rbind(allSubjects, allSubjects),
    "task"=gsub("Preprocessed", "", rep(task.dirs, each=nSubjects), fixed=TRUE),
    "e.norm"=e.norm.matrix[, 1])

write.csv(enorm.df, file="enorms.csv", quote=FALSE, row.names=FALSE)
