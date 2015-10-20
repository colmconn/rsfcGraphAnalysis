rm(list=ls())
graphics.off()
library(nlme)
library(ggplot2)

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

motionTableHeader=c("roll", "pitch", "yaw", "dS",  "dL",  "dP")

groups="mddAndCtrl"


group.vbm.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/vbm.subject.list.from.matthew.n114"))
    
ctrl.subjectOrder.filename=file.path(group.vbm.dir, "rs114a.ncl.txt")
mdd.subjectOrder.filename=file.path(group.vbm.dir,  "rs114a.mdd.txt")

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
## allSubjects[allSubjects$subject=="300_A", "Group"]="MDD"
##allSubjects$subject=droplevels(allSubjects$subject)

nSubjects=dim(allSubjects)[1]

censor.matrix=matrix(0, nrow=nSubjects, ncol=1)
censor.matrix.rownames=c()

subjectCount=1
for (subject in allSubjects$subject) {
    if (subject == "169/300_A") {
        subject="300_A"
    }
    censorFile=file.path(data.dir, subject, "rsfcPreprocessed", paste(subject, "pm.censor.1D", sep="."))
    cleanedEPI=file.path(data.dir, subject, "rsfcPreprocessed", sprintf("%s.pm.cleanEPI+orig.HEAD", subject))
    
    if (file.exists(censorFile) & file.exists(cleanedEPI)) {
        ## cat("*** Attempting to read", censorFile, "\n")
        censorTable=read.table(censorFile, col.names=c("censor"), header=FALSE)
        numberOfCensoredVolumes=sum(1-censorTable[,1])
        totalNumberOfVolumes=length(censorTable[,1])
        censor.matrix[subjectCount, 1] = numberOfCensoredVolumes
        censor.matrix.rownames=c(censor.matrix.rownames, subject)
        
    } else {
        ## cat("*** No such file", censorFile, "\n")
        censor.matrix[subjectCount, 1] = NA
        censor.matrix.rownames=c(censor.matrix.rownames, subject)
    }
    subjectCount=subjectCount+1
}
colnames(censor.matrix)="numberOfCensoredVolumes"
#censor.matrix.rownames=gsub("_A", "", censor.matrix.rownames, fixed=TRUE)
#rownames(censor.matrix)=censor.matrix.rownames


censor.complete.df=droplevels(as.data.frame(cbind(allSubjects, censor.matrix)))


cat("Number of subjects per group\n")

print(addmargins(table(censor.complete.df$Group)))

censor.test=wilcox.test(censor.complete.df[censor.complete.df$Group=="MDD", "numberOfCensoredVolumes"], censor.complete.df[censor.complete.df$Group=="NCL", "numberOfCensoredVolumes"])
print(censor.test)
## stop()

nSubjects=dim(censor.complete.df)[1]
##allSubjects$subjects=data.frame(subject = droplevels(allSubjects$subjects))

cat ("Motion Analyses\n")
for ( func in c("min", "mean", "max")) { 

    motion.matrix=matrix(0, nrow=nSubjects, ncol=6)
    motion.matrix.rownames=c()
    
    subjectCount=1
    for (subject in censor.complete.df$subject) {
        ## motionFile=file.path(data.dir, paste(subject, "A", sep="_"), "functional", "motion_demean.1D")
        if (subject == "169/300_A") {
            subject="300_A"
        }
        motionFile=file.path(data.dir, subject, "rsfcPreprocessed", "tmp", sprintf("%sfuncon.zp.float.despike_tsh_vr_motion.tcat.1D", subject))
        if (file.exists(motionFile)) {
            ## cat("*** Attempting to read", motionFile, "\n")
            motionTable=read.table(motionFile, header=FALSE)
            
            ## change mean here to max to get the maximum movement excursion
            ## instead of average
            motionExcursion=apply(motionTable, 2, func)
            motion.matrix[subjectCount, ] = motionExcursion
            motion.matrix.rownames=c(motion.matrix.rownames, subject)
            
        } else {
            ## cat("*** No such file", motionFile, "\n")
        }
        subjectCount=subjectCount+1
    }
    
    colnames(motion.matrix)=motionTableHeader
    motion.matrix.rownames=gsub("_A", "", motion.matrix.rownames, fixed=TRUE)
    rownames(motion.matrix)=motion.matrix.rownames
    
    motion.df=as.data.frame(motion.matrix)
    motion.df$subject=as.factor(rownames(motion.matrix))
    motion.df$Group=demographics[match(motion.matrix.rownames, demographics$ID), "Grp"]
    ##motion.df$BDI.II=demographics[match(motion.matrix.rownames, demographics$ID), "BDI.II"]
    motion.df$Gender=demographics[match(motion.matrix.rownames, demographics$ID), "Gender"]
    ##motion.df$inFinalAnalysis=ifelse(is.na(match(motion.df$subject, gsub("_A", "", subjectsWithTooFewTimePoints$subject, fixed=TRUE))), "included", "excluded")
    ##motion.df$inFinalAnalysis=as.factor(motion.df$inFinalAnalysis)
    
    ## now fix up the 300 subject cos somehow it's got 2 IDs.
    ##motion.df[motion.df$subject==300, "Group"]="MDD"
    ##motion.df[motion.df$subject==300, "Gender"]="F"
    ##includedInFinalAnalysis.df=motion.df[motion.df$inFinalAnalysis=="included", ]
    motion.complete.df=motion.df
    
    ## only the subjects included in the final analysis
    cat (sprintf("*** Only the subjects included in the final analysis: %s motion\n", func))
    ## Y=motion.matrix[motion.df$inFinalAnalysis=="included", 1:6]
    ## model.lme=lme(Y ~ Group, random = ~ 1 | subject, data=includedInFinalAnalysis.df)
    ## model.lme.anova=anova(model.lme)
    ## print(model.lme.anova)
    
    ## model2.lme=lme(cbind(roll, pitch, yaw, dS,  dL, dP) ~ Group, random = ~ 1 | subject, data=motion.complete.df)
    model2.lme=lme(roll ~ Group, random = ~ 1 | subject, data=motion.complete.df)    
    model2.lme.anova=anova(model2.lme)
    print(model2.lme.anova)
    
    ## cat ("*** All subjects\n")
    ## ##motion.df=motion.df[, c("subject", "Group", "Gender", "inFinalAnalysis")]
    ## model.lme=lme(cbind(roll, pitch, yaw, dS,  dL, dP) ~ Group*inFinalAnalysis, random = ~ 1 | subject, data=motion.df) 
    ## model.lme.anova=anova(model.lme)
    ## print(model.lme.anova)
}
