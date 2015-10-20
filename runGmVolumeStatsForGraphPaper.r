#!/usr/bin/env Rscript

rm(list=ls())
graphics.off()

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

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

fixSubjectIds <- function (inDataFrame, inColumnName=NULL) {
    if ( is.null(inColumnName) ) {
        stop("Got a NULL column name in fixSubjectIds. Stopping\n")
    }

    inDataFrame[, inColumnName]=gsub("_[ABCD]", "", as.character(inDataFrame[, inColumnName]), fixed=FALSE)
    subjectIds=as.character(inDataFrame[, inColumnName])
    if (length(pmatch(subjectIds, "300")) > 0 ) {
        subjectIds=sub("^300", "169/300", subjectIds, fixed=FALSE)
        subjectIds=sub("^0", "", subjectIds, fixed=FALSE)
    }
    return (subjectIds)
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
####################################################################################################
####################################################################################################

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
standard.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/standard"))
masks.dir=file.path(standard.dir, "yeo17liberal")
scriptsDir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))
group.vbm.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/vbm.subject.list.from.matthew.n111"))

## names of files important for the analysis below
gm.mod.merg.file=file.path(group.vbm.dir, "stats", "GM_mod_merg_s3.nii.gz")
anat.file.list=file.path(group.vbm.dir, "filelist.csv")
list.of.mask.files=file.path(masks.dir, "rs111_vbm_masks.txt")
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
wasiFilename=file.path(admin.data.dir, "WASI.csv")

## read the demographics and WASI data
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ""))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

cat("*** Reading", wasiFilename, "\n")
wasi=read.csv(wasiFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
cat(sprintf("*** Read WASI-II data for %s unique subjects\n",  length(unique(wasi$SubID))))


## check that files exist before tryign to proceeed with the analysis
if (! file.exists(gm.mod.merg.file) ) {
    stop("*** No such file", gm.mod.merg.file, "Cannot continue.\n")
}

if ( ! file.exists(anat.file.list) ) {
    stop("*** No such file", anat.file.list, "Cannot continue.\n")
}

if (! file.exists(list.of.mask.files) ) {
    stop("*** No such file", list.of.mask.files, "Cannot continue.\n")
}         

masks=scan(list.of.mask.files, what=character())
anatFiles=scan(anat.file.list, what=character())
subjects.to.groups=data.frame(do.call(rbind, lapply(strsplit(anatFiles, ".", fixed=TRUE), "[", 1:2)))
colnames(subjects.to.groups)=c("Group", "subjectNumber")
subjects.to.groups$subjectNumber=as.factor(gsub("_[ABCD]", "", as.character(subjects.to.groups$subjectNumber), fixed=FALSE))

for (mask in masks) {
    if (file.exists(mask)) {
        roistats.command=sprintf("3dROIstats -nobriklab -mask %s %s", mask, gm.mod.merg.file)
        roistats.command.output=system(roistats.command, intern=TRUE)

        ## skip the first row as it will be the header and we will
        ## need to handle it separately to ensure that the column
        ## names are named properly and that the numeric columns
        ## are not transformed to factors by the call to
        ## data.frame
        ## roistats=data.frame(do.call(rbind, strsplit(roistats.command.output[-1], "\t", fixed=TRUE)))
        ## colnames(roistats)=sapply(strsplit(roistats.command.output[1], "\t", fixed=TRUE), trim)
        ## now ensure that all columns not named File are numeric and not factors
        ## numeric.columns=grep("File", colnames(roistats), invert=TRUE)
        ## if(length(numeric.columns) > 0) {
        ##     roistats[, numeric.columns] = lapply(roistats[, numeric.columns], function (xx) { as.numeric(as.character(xx)) } )
        ## }


        roistats=read.table(textConnection(roistats.command.output), header=TRUE)
        ## dump the File column as it's not really useful
        roistats$File=NULL
        roistats=computeAge(fixDates(
            cbind(
                subjects.to.groups,
                roistats,
                demographics    [match(subjects.to.groups$subjectNumber, demographics$ID), c("Gender", "DOB", "MRI")],
                wasi            [match(subjects.to.groups$subjectNumber, wasi$SubID), c("Verbal", "Performance", "Full")]
            )))
                
        n.levels=length(levels(roistats$Group))
        if (n.levels == 2) {
            g.one=levels(roistats$Group)[1]
            g.two=levels(roistats$Group)[2]

            cat("####################################################################################################\n")
            cat("### Mask:", mask, "\n")
            cat(sprintf("*** Comparing the volume of GM in Mean_1 in %s (g.one) to %s (g.two)\n", g.one, g.two))
            ## print(wilcox.test(x = roistats[roistats$Group==g.one, "Mean_1"], y=roistats[roistats$Group==g.two, "Mean_1"]))
            ## demean age.in.years
            ## roistats$age.in.years=scale(roistats$age.in.years, center=TRUE, scale=FALSE)
            print(summary(aov(Mean_1 ~ Group + Gender + age.in.years + Full, data=roistats)))

            print(summarySE(roistats, measure="Mean_1", groupvars="Group", na.rm=TRUE))
        } else {
            stop("*** Sorry cant handle more than two groups\n")
        }
    } else {
        warning("*** Skipping non-existant mask file: ", mask)
    }
}


