rm(list=ls())
graphics.off()


## Reads the seed file and does the (crude) equivalent of BASH variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character(), quiet=TRUE)
    table=gsub("$DATA", seeds.data.dir, table, fixed=TRUE)

    return (table)
}

## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed
getSeedName <- function(inSeedPath){
    name=basename(inSeedPath)
    if (grepl("\\.nii", name)) {
        return(gsub("\\.nii.*", "", name))
    } else if (grepl("\\+tlrc", name)) {
        return(gsub("\\+tlrc.*", "", name))
    } else {
        return (name)
    }
}

readSubjectOrderTable <- function (inFilename) {
    cat("*** Reading", subjectOrderFilename, "\n")
    subjectOrder=read.csv(subjectOrderFilename, header=T)

    return(subjectOrder)
}

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
}

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
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
##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################


if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

## ncpus=1

scripts.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts")
standard.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/standard")
data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")

admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")

## group.data.dir=file.path(data.dir, "Group.data.acc.amygdala.paper")
## group.results.dir=file.path(data.dir, "Group.results.acc.amygdala.paper")
## analysis="acc.amygdala.paper"

## group.data.dir=file.path(data.dir, "Group.data.kaiser.amygdala.paper")
## group.results.dir=file.path(data.dir, "Group.results.kaiser.amygdala.paper")
## analysis="kaiser.amygdala.paper"



## group.data.dir=file.path(data.dir, "Group.data.acc.graph.paper")
## group.results.dir=file.path(data.dir, "Group.results.acc.graph.paper")
## analysis="acc.graph.paper"

##group.data.dir=file.path(data.dir, "Group.data.kaiser.graph.paper")
##group.results.dir=file.path(data.dir, "Group.results.kaiser.graph.paper")
##analysis="kaiser.graph.paper"




group.data.dir=file.path(data.dir, "Group.data.afni")
group.results.dir=file.path(data.dir, "Group.results.afni")
analysis=""

seeds.data.dir=file.path(data.dir, "seeds")

demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)

wasiFilename=file.path(admin.data.dir, "WASI.csv")
wasi=readCsvFile(wasiFilename, inSubjectColumnName="SubID")

## seeds=readSeedsFile(file.path(config.data.dir, "juelich_amygdala_seeds_weights.txt"))

## seeds=readSeedsFile(file.path(config.data.dir, "kaiser_supplemental_seeds.txt"))

## seeds=readSeedsFile(file.path(config.data.dir, "short_ACC_seed_list.txt"))

seeds=readSeedsFile(file.path(config.data.dir, "caez_juelich_whole_amygdala_seeds.txt"))

## wasi.column.names=c("Verbal", "Performance", "Full")
## only use the WASI Full score as a covariate as perfromance and
## verbal are highly correlated
wasi.column.names=c("Full")

## covariate.column.names=c("age.in.years")
covariate.column.names=c("Full")
## covariate.column.names=c("Full", "age.in.years", "Gender")

## atAndNat = suicide attempters and non-attempters
grouping="mddAndCtrl"
## grouping="atAndNat"

## new.ttest.arguments=""
## new.ttest.arguments="-resid %s \\\n-CLUSTSIM \\\n-prefix_clustsim resid.CStemp.%s \\\n-tempdir tmp.resid.clustsim.%s"
new.ttest.arguments="-resid %s \\\n-CLUSTSIM \\\n-prefix_clustsim resid.CStemp.%s"



for (seed in seeds) {
    seedName=getSeedName(seed)
    
    cat("####################################################################################################\n")
    cat(sprintf("*** Creating covariate file for the %s seed of the %s grouping\n", seedName, grouping))
    
    ## this file stores the order of the subjects in each of the following BRIK files
    subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", grouping, seedName, "csv", sep="."))
    if (file.exists(subjectOrderFilename) ) {
        subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))
    } else {
        cat("*** Could not find a subject order file. Trying to load subject lists from", admin.data.dir, "\n")

        attempters=read.table(file.path(config.data.dir, "mdd.at.pilot.txt"))
        nonAttempters=read.table(file.path(config.data.dir, "mdd.nat.pilot.txt"))

        ## subjectOrder=data.frame(
        ##     "subject"=c(as.character(attempters[, 1]), as.character(nonAttempters[, 1])))

        subjectOrder=data.frame(
                c(as.character(attempters[, 1]), as.character(nonAttempters[, 1])),
                c(rep("Attempter", times=length(attempters[, 1])),  rep("Non-attempter", times=length(nonAttempters[, 1])))
            )
        colnames(subjectOrder)=c("subject", "Group")
    }

    bucketListFilename=file.path(group.data.dir, paste("bucketList", grouping, seedName, "txt", sep="."))
    if (file.exists(bucketListFilename) ) {
        cat("*** Reading", bucketListFilename, "\n")
        bucketList=read.table(bucketListFilename, header=FALSE)
        ## mgd$InputFile=sub("z-score", "z-score.masked", mgd$InputFile, fixed=TRUE)
    } else {
        cat("*** The bucket list file does not exist. Making list of files from subject order list\n")
        bucketList=data.frame("InputFile" = sapply(subjectOrder$subject, function (xx) {
            file.path(data.dir, sprintf("%s/rsfc/%s/%s.z-score.masked+tlrc.HEAD", xx, seedName, seedName))
        }) )

        ## check that these files actually exist and quit 
        for (file in bucketList$InputFile) {
            if ( ! file.exists (file) ) {
                cat(sprintf("*** No such file : %s\n", file))
            }
        }
        
        ## colnames(bucketlist)="InputFile"
        ## print(bucketList)
        ## nwo we can fix the subjectOrder table as to generate the
        ## bucketFileList we need the _A and non-fixed 169/300 subject
        ## to get the correct file names
        subjectOrder=fixSubjectOrderTable(subjectOrder)
        ## print(subjectOrder)
    }

    if (! any(grepl("Group", colnames(subjectOrder))) ) {
        mgd=cbind(subjectOrder,
            bucketList,
            demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "DOB", "MRI", "Gender")],
            wasi        [match(subjectOrder$subject, wasi$SubID), wasi.column.names]
                  )
        colnames(mgd)=c("subject", "InputFile", "Group", "DOB", "MRI", "Gender", wasi.column.names)
    } else {
        mgd=cbind(subjectOrder,
            bucketList,
            demographics[match(subjectOrder$subject, demographics$ID), c("DOB", "MRI", "Gender")],
            wasi        [match(subjectOrder$subject, wasi$SubID), wasi.column.names]
                  )
        colnames(mgd)=c("subject", "Group", "InputFile", "DOB", "MRI", "Gender", wasi.column.names)
    }

    mgd[mgd$subject=="378", "Gender"]="F"
    mgd$subject=paste(mgd$subject, "A", sep="_")


    mgd=fixDates(mgd)
    mgd=computeAge(mgd)
    
    ## now center the numeric covariates

    for (col in covariate.column.names) {
        if ( is.numeric(mgd[, col]) ) {
            cat("*** Mean centering", col, "\n")            
            mgd[, col] = scale(mgd[, col], center=TRUE, scale=FALSE)
        } else {
            ##cat(sprintf("*** Skipping centering for %s: Not a numeric column\n", col))
            cat(sprintf("*** WARNING: Converting %s to underlying numeric representation of the R factor and mean centering that\n", col))
            mgd[, col] = scale(as.numeric(mgd[, col]), center=TRUE, scale=FALSE)            
        }
    }

    rownames(mgd)=NULL

    mgd=droplevels(mgd)

    ## now reorder the columns of mgd to InputFile is last and only
    ## pick those columns that we want to correct for in the t-tests
    mgd=mgd[, c("subject", "Group", covariate.column.names, "InputFile")]

    covariates.filename=file.path(group.data.dir, gsub("..", ".", paste("3dttest.covariates", grouping, seedName, analysis, "txt", sep="."), fixed=TRUE))
    cat("*** Writing covariates file to:", covariates.filename, "\n")
    ## the ordering of the columns in the the write command below is
    ## important. It must be subject <covariates> InputFile.

    write.table(mgd[, c("subject", covariate.column.names)], file=covariates.filename, quote=FALSE, col.names=TRUE, row.names=FALSE, eol="\n")

    ## cat("*** The data table is as follows:\n")
    ## print(head(mgd))

    ## mask=sprintf("../mask.grey.%s.union.masked+tlrc.HEAD", grouping)
    ## mask=sprintf("../mask.grey.%s.union.masked+tlrc.HEAD", "mddAndCtrl")
    mask=file.path(standard.dir, "/MNI_caez_N27_brain.3mm+tlrc.HEAD")
    prefix=sprintf("new.ttest.%s.%s.covaried", grouping, seedName)
    resid.prefix=sprintf("new.ttest.%s.%s.covaried.resid", grouping, seedName)    

    ## setA.group="Attempter"
    setA.group=levels(mgd$Group)[1]
    setA.files=paste(mgd[mgd$Group==setA.group, "InputFile"], collapse=" \\\n")
    setA.labels.and.files = paste(apply(mgd[mgd$Group==setA.group, c("subject", "InputFile")], 1, paste, collapse=" "), collapse=" \\\n")
    
    ## setB.group="Non-attempter"
    setB.group=levels(mgd$Group)[2]
    setB.files=paste(mgd[mgd$Group==setB.group, "InputFile"], collapse=" \\\n")
    setB.labels.and.files = paste(apply(mgd[mgd$Group==setB.group, c("subject", "InputFile")], 1, paste, collapse=" "), collapse=" \\\n")
    
    three.d.ttest.command = "3dttest++"
    ## use center NONE here as the covariates are already mean
    ## centered before they were written to the covariates file
    if (nchar(new.ttest.arguments) > 0) {
        ## sprintf(new.ttest.arguments, resid.prefix, seedName, seedName),
        ## three.d.ttest.arguments = sprintf("%s \\\n-mask %s \\\n-prefix %s \\\n-center NONE \\\n-setA %s %s \\\n-setB %s %s \\\n-covariates %s",
        three.d.ttest.arguments = sprintf("%s \\\n-mask %s \\\n-prefix %s \\\n-center NONE \\\n-setA %s %s \\\n-setB %s %s \\\n-covariates %s",            
            sprintf(new.ttest.arguments, resid.prefix, seedName),            
            mask, prefix,
            setA.group, setA.labels.and.files,
            setB.group, setB.labels.and.files,
            covariates.filename)
    } else {
        three.d.ttest.arguments = sprintf("-mask %s \\\n-prefix %s \\\n-center NONE \\\n-setA %s %s \\\n-setB %s %s \\\n-covariates %s",
            mask, prefix, setA.group, setA.labels.and.files, setB.group, setB.labels.and.files, covariates.filename)
    }
    
    three.d.ttest.command.script.filename=file.path(scripts.dir, gsub("..", ".", sprintf("afni-07-ttest.withCovariates.%s.%s.%s.sh", grouping, seedName, analysis), fixed=TRUE))
    cat("*** Writing the 3dttest++ command to:", three.d.ttest.command.script.filename, "\n")

    ## full.three.d.ttest.command=sprintf("cd %s ; mkdir tmp.resid.clustsim.%s ; %s %s", group.results.dir, seedName, three.d.ttest.command, three.d.ttest.arguments)
    full.three.d.ttest.command=sprintf("#!/bin.bash\n\nunset AFNI_COMPRESSOR \nexport OMP_NUM_THREADS=40\nmkdir -p %s/new.ttest.%s \ncd %s/new.ttest.%s \n%s %s",
        group.results.dir, seedName,
        group.results.dir, seedName,        
        three.d.ttest.command, three.d.ttest.arguments)    
    cat (full.three.d.ttest.command, file=three.d.ttest.command.script.filename)
    ## cat(full.three.d.ttest.command, "\n\n")
    ## stop()

}
