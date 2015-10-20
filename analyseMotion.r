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

####################################################################################################
### This variable controls whether the analysis should be done for the
### VBM analysis (TRUE) of the RSFC analysis (FALSE)
####################################################################################################
vbmAnalysis=FALSE

if (! vbmAnalysis) {

    seed.name="L_BLA.weight.3mm"

    ## this file stores the order of the subjects in each of the following BRIK files
    ctrl.subjectOrder.filename=file.path(group.data.dir, paste("subjectOrder.ctrlOnly", seed.name, "csv", sep="."))
    mdd.subjectOrder.filename=file.path(group.data.dir, paste("subjectOrder.mddOnly",seed.name, "csv", sep="."))
    
} else {
    ## analysis of the VBM subjects
    
    group.vbm.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/vbm.subject.list.from.matthew.n114"))
    
    ctrl.subjectOrder.filename=file.path(group.vbm.dir, "ncl.subjectlist.txt")
    mdd.subjectOrder.filename=file.path(group.vbm.dir,  "mdd.subjectlist.txt")
}

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

mddSubjects=read.table(file.path(config.data.dir, "clean.mdd.subjectList.txt"), header=FALSE)
controlSubjects=read.table(file.path(config.data.dir, "clean.ncl.subjectList.txt"), header=FALSE)
bothGroups=rbind(mddSubjects, controlSubjects)
allSubjects=data.frame(
    "subject" = bothGroups,
    "Group"    = demographics[match(gsub("_A[0-9]?", "", bothGroups$V1, fixed=FALSE), demographics$ID), "Grp"]
    )
colnames(allSubjects)=c("subject", "Group")
allSubjects[allSubjects$subject=="300_A", "Group"]="MDD"
allSubjects$subject=droplevels(allSubjects$subject)

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
## stop()

excessiveMotionThreshold=0.2
excessiveMotionThresholdPercentage=excessiveMotionThreshold * 100

censor.df=as.data.frame(cbind(allSubjects, censor.matrix))
censor.df$drop=ifelse(is.na(censor.df$numberOfCensoredVolumes), "no data", ifelse(censor.df$numberOfCensoredVolumes > excessiveMotionThreshold*totalNumberOfVolumes, "drop", "keep"))
censor.df$drop=as.factor(censor.df$drop)
censor.df=droplevels(censor.df)

censor.df$inSubjectOrder = censor.df$subject %in% subjectOrder$subject

censor.complete.df=censor.df[censor.df$drop == "keep" , ]

cat("Number of subjects per group\n")

print(table(censor.complete.df$Group))

cat("The following subjects should be dropped do to excessive motion/censoring:\n")
print(rownames(censor.df)[censor.df$drop=="drop"])
print(as.vector(censor.df[censor.df$drop=="drop", "Group"]))

cat ("*** Only the subjects that can be included in the final analysis\n")
### cat (" *** NOTE ***\nThis table does not account for subjects dropped because of missing data.\nTo get the final N per group for your analysis you should use analyseGroupStatsOnly.r and use the tables there\n*** NOTE ***\n")
print(with(censor.df, addmargins(table(Group, drop))))

censor.df2=droplevels(subset(censor.df, drop %in% c("drop", "keep")))
print(with(censor.df2, addmargins(table(Group, drop))))

cat("*** Test of the whether the proportion of subjects dropped due to motion differed between the groups\n")
print(prop.test(with(censor.df2, table(Group, drop))))

cat("\n*** NOTE ***\nStatistics refer only to those subjects who are not dropped due to excessive motion and outlier count\n*** NOTE ***\n")

censor.test=wilcox.test(censor.complete.df[censor.complete.df$Group=="MDD", "numberOfCensoredVolumes"], censor.complete.df[censor.complete.df$Group=="NCL", "numberOfCensoredVolumes"])
print(censor.test)

my.base.size=14
my_theme =
    theme_bw(base_size =  my.base.size) +
        theme(##title = "Proportion of risky choices",
            legend.position="bottom",
            ##panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            ##axis.title=element_blank(),\
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            ##axis.title.x = element_text(size=my.base.size, vjust=0),
            ##axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
            plot.title=element_text(size=my.base.size*1.2, vjust=1))

censor.boxplot=ggplot(censor.df, aes(x=Group, numberOfCensoredVolumes, fill=Group)) +
    geom_boxplot() +
    xlab("Group") +
    ylab("Number of volumes censored") +
    ggtitle("Volumes censored in the RSFC analysis") +
    scale_fill_brewer(palette="Set1") +
    my_theme

## dev.new(); print(censor.boxplot)
censor.barplot=ggplot(censor.df, aes(x=subject, fill=Group, y=numberOfCensoredVolumes)) +
    geom_bar(stat="identity") +
        scale_fill_brewer(palette="Set1") +
            geom_hline(yintercept=excessiveMotionThreshold * totalNumberOfVolumes, color="black") +
                annotate("text", x=excessiveMotionThresholdPercentage, y=(excessiveMotionThreshold*totalNumberOfVolumes) + 4,
                         label=sprintf("%d%% (n=%d)",
                             round(excessiveMotionThresholdPercentage, 0),
                             round(excessiveMotionThreshold*totalNumberOfVolumes, 0))) +
                    my_theme +
                        labs(x="Subject", y="Number of volumes censored") +
                            ggtitle("Volumes censored in the RSFC analysis")

dev.new(); print(censor.barplot)


nSubjects=dim(censor.complete.df)[1]
##allSubjects$subjects=data.frame(subject = droplevels(allSubjects$subjects))

cat("Number of subjects per group\n")

print(table(censor.complete.df$Group))

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
