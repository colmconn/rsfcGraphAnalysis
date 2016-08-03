rm(list=ls())
graphics.off()

filterOutStimulusColumn <- function (inDataFrame) {
    df=inDataFrame
    if (dim(inDataFrame)[2] > 1) {
        cat("*** Filtering the data frame because it has more then 1 column\n")
        df=data.frame("subject"=unique(inDataFrame[, 1]))
    }
    return (df)
}


build.tsnr.filenames <- function (in.subjects, in.seed.name) {

    filenames=sapply(in.subjects,
        function(x) {
            ## /data/sanDiego/rsfcGraphAnalysis/data/105_A/rsfcPreprocessed/105_A.MNI152_T1_3mm_brain_mask.tsnr.1D
            
            sprintf("%s/%s/rsfcPreprocessed/%s.%s.tsnr.1D", data.dir, x, x, in.seed.name)
        })
    
    return(filenames)
}

## http://wiki.stdout.org/rcookbook/Graphs/Plotting%20means%20and%20error%20bars%20(ggplot2)/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    ## New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    ## This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                       c(N    = length2(xx[,col], na.rm=na.rm),
                         mean = mean (xx[,col], na.rm=na.rm),
                         median = median (xx[,col], na.rm=na.rm),
                         IQR = IQR (xx[,col], na.rm=na.rm),                         
                         mad = mad (xx[,col], na.rm=na.rm),                       
                         sd   = sd   (xx[,col], na.rm=na.rm),
                         min  = min  (xx[,col], na.rm=na.rm),
                         max  = max  (xx[,col], na.rm=na.rm),                       
                         nacount  = sum  (is.na((xx[,col])))
                         )
                   },
                   measurevar,
                   na.rm
                   )
    
    ## Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  ## Calculate standard error of the mean
    
    ## Confidence interval multiplier for standard error
    ## Calculate t-statistic for confidence interval: 
    ## e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}

####################################################################################################
## END OF FUNCTIONS
####################################################################################################

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=file.path(data.dir, "admin")
group.data.dir=file.path(data.dir, "Group.data")


demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")

subjectOrder.filename=file.path(group.data.dir, "subjectOrder.mddAndCtrl.L_whole_amygdala.3mm.csv")

subject.order=filterOutStimulusColumn(read.csv(subjectOrder.filename, header=TRUE))
colnames(subject.order)=c("subject")
subject.order$ID=gsub("([0-9]{3})_([ABCDE])", "\\1", subject.order$subject, fixed=FALSE)
cat(sprintf("Read %d subjects from the subject bucket list\n", dim(subject.order)[1]))

cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ""))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

mgd=mgd=cbind(subject.order, demographics[match(subject.order$ID, demographics$ID), c("Grp", "Gender")])
rownames(mgd)=NULL

seeds=c("L_whole_amygdala.3mm",  "MNI152_T1_3mm_brain_mask",  "R_whole_amygdala.3mm")

for (seed in seeds) {
    tsnr.filenames=build.tsnr.filenames(mgd$subject, seed)
    mgd[, seed]=vapply(tsnr.filenames, function (x) scan(x, quiet=TRUE), numeric(1))
}


sm.df=list()
for (seed in seeds) {
    cat("*** Summary for", seed, "\n")

    sm.df[[seed]]=summarySE(mgd, measure=seed, groupvars="Grp", na.rm=TRUE)
    print(sm.df[[seed]])
}


library(reshape2)
library(nlme)
## mgd.melted=melt(mgd, id.vars=c("subject", "Grp", "ID", "Gender"), measure.vars=seeds, variable.name="seed", value.name="tsnr")
mgd.melted=melt(mgd, id.vars=c("subject", "ID" ), measure.vars=seeds, variable.name="seed", value.name="tsnr")

## model.formula=as.formula("tsnr~Grp*seed")
model.formula=as.formula("tsnr~seed")
random.formula=as.formula("~ 1 | subject")
model.lme=lme(fixed=model.formula, data=mgd.melted, random=random.formula)

model.anova=anova(model.lme)

cat("\n### LME summary\n")
print(summary(model.lme))
cat("\n### ANOVA summary\n")
print(model.anova)
library(effects)
print(Effect("seed", model.lme))

sm.df.combined=do.call(rbind, lapply(names(sm.df),
    function (x) {
        sm.df[[x]]$seed=rep(x, 2);
        names(sm.df[[x]])[3] = "tsnr"
        sm.df[[x]]
    }))

library(ggplot2)
graph=ggplot(sm.df.combined, aes(x=seed, y=tsnr, color=Grp, fill=Grp))
graph=graph + geom_bar(position=position_dodge(), stat="identity") 
graph=graph + geom_errorbar(aes(ymin=tsnr-se, ymax=tsnr+se), width=.2,
    position=position_dodge(.9), color="black")
graph=graph + labs(x="Seed", y="TSNR")
print(graph)
