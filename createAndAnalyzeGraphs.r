rm(list=ls())
graphics.off()

library(ppcor)


####################################################################################################
### Function definitions
####################################################################################################

loadRegionOfInterestTimeseries <- function (inSubjects) {
    ## stores the names of the percentage change file names to be read    
    fnames=sapply(inSubjects,
        function(x) {
            rsfcFile=file.path(rsfc.data.dir, sprintf("%s_%s/%s_%s.rsfc.1D", x, timePoint, x, timePoint))
            return(rsfcFile)
        })
    
    ## check that these files actually exist and quit 
    for (file in fnames) {
        if ( ! file.exists (file) ) {
            stop(sprintf("*** No such file : %s\n", file))
        }
    }
    
    roiVolumes=do.call(rbind, lapply(fnames,
        function(x) {
            subject=sub(".*/([0-9]*)_[ABCDEF].*", "\\1", x, fixed=FALSE)
            rsfcTimeseries=read.table(x)
            ## now filter out only the subcortical clusters we're interest in, i.e., the left and right amygdalaa and hippocampii
            rsfcTimeseries$subject=rep(subject, dim(rsfcTimeseries)[1])

            return(rsfcTimeseries)
        }
        ))
    return (roiVolumes)
}



####################################################################################################
### Variable definitions
####################################################################################################


if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

rsfc.admin.data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/admin")  
rsfc.data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")  
rsfc.results.dir=file.path(rsfc.data.dir, "results")

if ( ! file.exists(rsfc.results.dir) ) {
    cat("*** Creating", rsfc.results.dir, "\n")
    dir.create(rsfc.results.dir)
}

groups="mddAndCtrl"
timePoint="A"

####################################################################################################
### Variable definitions
####################################################################################################

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(rsfc.admin.data.dir, "0-data_entry_current_10152013.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

subjects=c("105", "106")

roi.timeseries=loadRegionOfInterestTimeseries(subjects)
if ( any(grep("300", roi.timeseries$subject)) ) {
    roi.timeseries$subject[which (roi.timeseries$subject=="300")]="169/300"
}
