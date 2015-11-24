rm(list=ls())
graphics.off()

library(compute.es)
library(reshape)


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

makeTableString <- function(inGroup, inMean, inSe, inMin, inMax, inNaCount, inMissingData=TRUE) {
    ##  st=paste(round(inMean, 1), " / ", round(inSe, 1),    

    st=paste(round(inMean, 1), " ± ", round(inSe, 1),
        " (", round(inMin, 1), "-", round(inMax, 1), ")", ifelse(inMissingData & inNaCount > 0, paste(" [", inNaCount, "]", sep=""), ""), sep="")
    return(st)
}

fixDates <- function (inData) {
    ## this complicated looking regexp stuff cleans up the years in DOB
    ## and date with 4 digits to be just 2 digits
    ## month day year
    inData$date=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$date)
    
    ## now convert to year/month/day
    inData$date=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$date)
    
    inData$date=as.Date(inData$date, "%y/%m/%d")

    return(inData)
}


scan.dates=read.csv("../data/Group.results/scanDates.csv", na.strings=c("ND", "NF"))

scan.dates.from.data.spreadsheet=read.csv("../data/Group.results/scanDatesFromDataSpreadsheet.csv", na.strings=c("n/a"))

scan.dates$date=as.character(scan.dates$date)
scan.dates$date=as.Date(scan.dates$date, format="%Y%m%d")
scan.dates$timepoint=as.character(scan.dates$timepoint)
## scan.dates$subject=as.factor(scan.dates$subject)


colnames(scan.dates.from.data.spreadsheet)=c("subject", "timepoint", "date")
scan.dates.from.data.spreadsheet$date=as.character(scan.dates.from.data.spreadsheet$date)
scan.dates.from.data.spreadsheet=fixDates(scan.dates.from.data.spreadsheet)
scan.dates.from.data.spreadsheet$timepoint=as.character(scan.dates.from.data.spreadsheet$timepoint)

## can't use the merge command for this as the dates in the spreadsheet
## do not match those in the DICOM files. So we consider the DICOMs
## canonical, uhhhhh
## scan.dates=merge(scan.dates, scan.dates.from.data.spreadsheet)

cat("*** Looking up missing scan dates in the scan.dates.from.data.spreadsheet\n")
for (ii in 1: dim(scan.dates)[1]) {
    if (is.na(scan.dates[ii, "date"])) { 
        ## try to lookup the subject and timepoint in the data from
        ## the spreadsheet maintained by the RAs

        cat("*** Looking for", scan.dates[ii, "subject"], "at timepoint", scan.dates[ii, "timepoint"], "in scan.dates.from.data.spreadsheet. Got:")
        scan.dates[ii, "date"]=scan.dates.from.data.spreadsheet[
                      scan.dates.from.data.spreadsheet$subject==scan.dates[ii, "subject"] &
                          scan.dates.from.data.spreadsheet$timepoint==scan.dates[ii, "timepoint"], "date"]
        print(scan.dates[ii, "date"])
    }
}


scan.dates.wide=reshape(scan.dates,
    timevar="timepoint",
    idvar=c("subject", "group"),
    direction="wide")
rownames(scan.dates.wide) = NULL


scan.dates.wide$date.diff=difftime(scan.dates.wide$date.C, scan.dates.wide$date.A, units="days")
scan.dates.wide$date.diff=as.numeric(scan.dates.wide$date.diff)

variable="date.diff"
group.variable="group"
group1="mdd"
group2="ncl"
sm.df=summarySE(scan.dates.wide, measure=variable, groupvars="group", na.rm=TRUE)
print(sm.df)

group1.string=makeTableString(sm.df[1, 1], sm.df[1, variable],  sm.df[1, "sd"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
group2.string=makeTableString(sm.df[2, 1], sm.df[2, variable],  sm.df[2, "sd"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)

dates.t.test=t.test (date.diff ~ group, data=scan.dates.wide, na.action="na.omit")
print(dates.t.test)

es=tes(dates.t.test$statistic, sm.df[1, "N"], sm.df[2, "N"], verbose=FALSE)
## upper bound on the 95% Confidence interval for the effect size
es.ci.lb=round(es$l.g, 2)
## upper bound on the 95% Confidence interval for the effect size                        
es.ci.ub=round(es$u.g, 2)
var.effect.size=round(es$g, 2)

cat("\n")
cat(sprintf("%s,%s,t(%0.2f) = %0.2f,%0.2f,g=%0.2f (%0.2f; %0.2f)\n", group1.string, group2.string, dates.t.test$parameter, dates.t.test$statistic, dates.t.test$p.value, var.effect.size, es.ci.lb,es.ci.ub ))
