rm(list=ls())
graphics.off()

splitSubjectOrderIntoIdAndTimepoint <- function(inSubjectOrderTable) {

    new.subject.order.table=cbind(inSubjectOrderTable, as.data.frame(t(as.data.frame(strsplit(as.character(inSubjectOrderTable$subject), "_", fixed=TRUE)))))
    rownames(new.subject.order.table)=NULL
    colnames(new.subject.order.table)=c("id", "subject", "timepoint")

    return(new.subject.order.table)
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

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

scripts.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))
data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=normalizePath(file.path(data.dir, "admin"))
config.data.dir=normalizePath(file.path(data.dir, "config"))
group.data.dir=normalizePath(file.path(data.dir, "Group.data"))
group.results.dir=normalizePath(file.path(data.dir, "Group.results.RADS.Total.Tscore.diff.withAandC"))
seeds.data.dir=normalizePath(file.path(data.dir, "seeds"))

demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "", "Test_Subject"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

rads.time.diff.filename=file.path(admin.data.dir, "RADS.Total.Tscore.change.csv")
rads.time.diff=read.csv(rads.time.diff.filename, header=TRUE)

subjects=do.call(rbind, lapply(
    file.path(config.data.dir, c("mdd.subjectList.with.RADS.Total.Tscore.AandC.txt",
                                 "ncl.subjectList.with.RADS.Total.Tscore.AandC.txt")),
    function(X) {
        data.frame(read.table(X, header=FALSE, sep=""))
    }
))
colnames(subjects)=c("subject")
subjects=splitSubjectOrderIntoIdAndTimepoint(subjects)
subjects=droplevels(
    cbind(subjects,
          rads.time.diff[match(subjects$subject, rads.time.diff$SubjNum), ],
          demographics[match(subjects$subject, demographics$ID), c("Grp")]          
          ))
colnames(subjects)=c(colnames(subjects)[1:8], "group" )

print(addmargins(table(subjects$group)))

bp.graph=ggplot(subjects, aes(x=group, y=RADS.Total.Tscore.diff)) +
    geom_boxplot() +
    geom_jitter(alpha=0.5) +
        labs(x="Group", y="Reynolds Adolescent Depression Scale Total (Standardized; A to C Change)") +
            scale_fill_brewer(palette="Set1")

print(bp.graph)
image.file.name=file.path(group.results.dir, "rads.tscore.diff.fromAtoC.pdf")
cat("*** Writing", image.file.name, "\n")
ggsave(image.file.name, bp.graph)
print(summarySE(subjects, measurevar="RADS.Total.Tscore.diff", groupvars="group"))





