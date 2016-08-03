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
group.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.data.CDRS.t.score.scaled.diff"))
group.results.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.results.CDRS.t.score.scaled.diff"))
scriptsDir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))

mdds.only.filename=file.path(group.results.dir, "roiStats.regression.fwhm4.2.restingstate.mddOnly.R_whole_amygdala.3mm.and.CDRS.t.score.scaled.diff.txt")
mdds.only.df=read.table(mdds.only.filename, header=TRUE)
mdds.only.df=mdds.only.df[, -1]

mdds.only.subject.order.filename=file.path(group.data.dir, "subjectOrder.mddOnly.R_whole_amygdala.3mm.csv")
mdds.only.subject.order=read.table(mdds.only.subject.order.filename, header=TRUE)
mdds.only.df=cbind(mdds.only.subject.order, mdds.only.df)

mdds.only.df$Sub.brick=NULL
names(mdds.only.df)=c("subject", "R.Insula")
mdds.only.df$Group=rep("MDD", dim(mdds.only.df)[1])
mdds.only.df$Study.ID=gsub("_A", "", mdds.only.df$subject, fixed=TRUE)

change.score.filename=file.path(admin.data.dir, "new.mdd.CDRS.t.score.scaled.diff.change.score.csv")
change.scores=read.csv(change.score.filename, header=TRUE)
change.scores$Grp=NULL
change.scores$age.in.years=NULL
mdds.only.df=merge(mdds.only.df, change.scores, by.x="Study.ID", by.y="ID")

mdds.only.df$Study.ID=NULL
names(mdds.only.df)[1]="Study.ID"

mdds.only.df$Symptom.Change = rep("No Change", dim(mdds.only.df)[1])
mdds.only.df[mdds.only.df$CDRS.t.score.scaled.diff < 0, "Symptom.Change"] = "Improved"
mdds.only.df[mdds.only.df$CDRS.t.score.scaled.diff > 0, "Symptom.Change"] =  "Worsen"

cat("*** Table of symptom change over time\n")
print(addmargins(table(mdds.only.df$Symptom.Change)))

## now dump the Symptom.Change column because it will cause problems
## later when trying to merge the MDD and HCL data frames as the
## latter has no such column
mdds.only.df$Symptom.Change=NULL


## now filter to include only those MDDS whose CDRSR scores improved
## over time
mdds.only.df=subset(mdds.only.df, CDRS.t.score.scaled.diff < 0)


ctrls.only.filename=file.path(group.results.dir, "r_amy_r_insula_ctrlsOnly.1D")
ctrls.only.df=read.table(ctrls.only.filename, header=TRUE)
ctrls.only.df=ctrls.only.df[, -1]

ctrls.only.df$Sub.brick=gsub("0\\[(.*)\\]", "\\1", ctrls.only.df$Sub.brick)
names(ctrls.only.df)=c("Study.ID", "R.Insula")
ctrls.only.df$Group=rep("HCL", dim(ctrls.only.df)[1])

ctrls.only.df$CDRS.t.score.diff=rep(NA, dim(ctrls.only.df)[1])
ctrls.only.df$CDRS.t.score.scaled.diff=rep(NA, dim(ctrls.only.df)[1])

complete.df=rbind(mdds.only.df, ctrls.only.df)

library(ggplot2)
my.base.size=14
my.theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        legend.position="none",
        ## legend.position="bottom",        
        ## panel.grid.major = element_blank(),
        ## panel.grid.minor = element_blank(),

        ##remove the panel border
        ## panel.border = element_blank(),

        ## add back the axis lines
        axis.line=element_line(colour = "grey50"),
        
        ##axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))


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

my.dodge=position_dodge(.2)

data.sumamry.df=summarySE(complete.df, measure="R.Insula", groupvars=c("Group"), na.rm=TRUE)
print(data.sumamry.df)
graph=ggplot(complete.df, aes_string(x="Group", y="R.Insula", color="Group", fill="Group", shape="Group"))
graph=graph + stat_summary(fun.y = mean, geom="point", color="black", size=2, position=my.dodge)
graph=graph + geom_errorbar(data=data.sumamry.df, aes_string(ymin=paste("R.Insula", "se", sep="-"), ymax=paste("R.Insula", "se", sep="+")), width=.2, position=my.dodge)
graph=graph + geom_jitter(size=2)
graph=graph + labs(x="Group", y="Right Amygdala - Right Insula\nConnectivity")
graph=graph + scale_color_brewer(palette="Set1")
graph=graph + my.theme

print(graph)
imageFilename=file.path(group.results.dir, "r_amy_r_insula.mdd.and.hcl.connectivity.pdf")
cat(paste("*** Creating", imageFilename, "\n"))
ggsave(imageFilename, graph, width=4, height=3.5)

print(with(complete.df, t.test(R.Insula ~ Group)))
