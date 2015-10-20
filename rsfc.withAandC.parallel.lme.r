rm(list=ls())
gc()
##the name of this script
script.name=parent.frame(2)$ofile
## the location (absolute path) to this script
script.location=normalizePath(dirname(parent.frame(2)$ofile))

library(nlme)
library(contrast)

verbose=F

ncpus=1
if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
    ncpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
    ncpus=8
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))    
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

scripts.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))
data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=normalizePath(file.path(data.dir, "admin"))
config.data.dir=normalizePath(file.path(data.dir, "config"))
#group.data.dir=normalizePath(file.path(data.dir, "Group.data.withAandC"))
group.data.dir=normalizePath(file.path(data.dir, "Group.data"))
seeds.data.dir=normalizePath(file.path(data.dir, "seeds"))
#group.results.dir=normalizePath(file.path(data.dir, "Group.results.withAandC"))
group.results.dir=normalizePath(file.path(data.dir, "Group.results/afsp"))

AFNI_R_DIR=Sys.getenv("AFNI_R_DIR", unset=NA)
## use the functions for loading and saving briks from the AFNI
## distribution as they can cleanly handle floats/doubles
if ( ! is.na(AFNI_R_DIR) ) {
    source(file.path(AFNI_R_DIR, "AFNIio.R"))
} else {
    Sys.setenv(PATH=system("bash -login -c 'echo $PATH'", intern=TRUE))
    source(file.path(Sys.getenv("MRI_SOFTWARE_ROOT"), "afni", "AFNIio.R"))
}

## flag to indincate whether a progress bar should be printed or one line per slice if no progress bar is to be used
useProgressBar=TRUE

## this function takes the row and column names from an anova and
## combines them (along row) to get a list of labels for later use as
## labels to the 3drefit command so the briks in the bucket are
## correctly labeled
makeBrikLabels <- function (inAnova) {
    rns=rownames(inAnova)
    cns=colnames(inAnova)
    
    labels=c()
    
    for (i in 1:length(rns) ) {
        for (j in 1:length(cns) ) {
            name=paste(gsub("[(]|[)]", "", rns[i]), cns[j], sep=".")
            labels=c(labels, name)
        }
    }
    return (labels)
}

makeContrastBrikLabels <- function(inContrastsList) {
    return (unlist(lapply(inContrastsList, function(x) { paste(x$name, c("contrast", "t-value"), sep=".") } )))
}

## This function generates lists if indices that can be applied to an
## unlisted anova to extract the elements of said list in the correct
## order for writing to an AFNI BRIK. It works for odd and even
## numbers of rows.

makeAnovaIndices <- function (inAnova, inFirstTwoOnly=FALSE) {
    nrows=length(rownames(inAnova))
    ncols=length(colnames(inAnova))
    indices=c()
    for ( r in 1:nrows ) {
        colCount=0
        for ( cl in 1:ncols ) {
            index = (r + (cl)*nrows)-nrows
            if (inFirstTwoOnly) {
                if (colCount < 2) {
                    indices=c(indices, index)
                }
            }
            else {
                indices=c(indices, index)
            }
            colCount=colCount+1
        }
    }
    return (indices)
}

makeAfniFtestBrikIds <- function(inAnova) {
    nrows=length(rownames(inAnova))
    ncols=length(colnames(inAnova))
    
    return(seq(2, nrows*ncols, by=4))
}

## this function generates a numeric list of BRIK IDs (to be used when
## setting statpars in the call to 3drefit) from an anova and the list
## of preplanned contrasts
makeAfniTtestBrikIds <- function(inAnova, inContrastsList)  {
    numberOfContrasts=length(inContrastsList)
    if ( numberOfContrasts > 0 ) {
        numberOfAnovaBriks=length(colnames(inAnova)) * length(rownames(inAnova))
        return(seq(numberOfAnovaBriks + 1, numberOfAnovaBriks+(numberOfContrasts*2), by=2))
    } else {
        return( c() )
    }
}

## makes a appropriate statpar command line argument for use with
## 3drefit. Requires an anova object and a list of BRIK IDs
## corresponding to the list of f-stat briks. The t-test brik IDs
## (from the contrasts) and the corresponding degrees of freedom are
## optional and default to NULL indicating that there are no contrast
## BRIKs in the bucket file saved by this script
makeAfniStatparArguments <- function(inAnova, inAfniFtestBrikIds, inAfniTtestBrikIds=NULL, inContrastsDf=NULL) {
    anovaIndices=makeAnovaIndices(inAnova, inFirstTwoOnly=TRUE)
    temp = unlist(inAnova)
    indices=cbind(inAfniFtestBrikIds, matrix(anovaIndices, ncol=2, byrow=TRUE))
    anovaStatpar=paste(apply(indices, 1, function(x) { sprintf("-substatpar %d fift %d %d", x[1], temp[x[2]], temp[x[3]]) }), collapse=" ")
    contrastsStatpar=""
    if (! is.null(inAfniTtestBrikIds) ) {
        indices=cbind(inAfniTtestBrikIds, inContrastsDf)
        contrastsStatpar=paste(apply(indices, 1, function(x) { sprintf ("-substatpar %d fitt %d", x[1], x[2])} ), collapse=" ")
    }
    return(paste(anovaStatpar, contrastsStatpar))
}

stopIfFormulaTermsNotInModelMatrix <- function(inModelFormula, inRandomFormula, inModelMatrix) {

    termDifference=setdiff(c(all.vars(inModelFormula), all.vars(inRandomFormula)), colnames(inModelMatrix))
    if (length(termDifference) > 0) {
        stop(paste("The following terms were in the model forumlae but are not columns of the model matrix:", termDifference))
    }   
}

## Reads the seed file and does the (crude) equivalent of BAS variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character())
    table=gsub("$DATA", seeds.data.dir, table, fixed=TRUE)

    return (table)
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

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}

readSubjectOrderTable <- function (inFilename) {
    cat("*** Reading", subjectOrderFilename, "\n")
    subjectOrder=read.csv(subjectOrderFilename, header=T)

    return(subjectOrder)
}

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("300", "169/300", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
}

splitSubjectOrderIntoIdAndTimepoint <- function(inSubjectOrderTable) {

    new.subject.order.table=cbind(inSubjectOrderTable, as.data.frame(t(as.data.frame(strsplit(as.character(inSubjectOrderTable$subject), "_", fixed=TRUE)))))
    rownames(new.subject.order.table)=NULL
    #colnames(new.subject.order.table)=c("id", "subject", "timepoint")
    colnames(new.subject.order.table)=c("id", "subject")

    return(new.subject.order.table)
}

####################################################################################################
## runLme function definition

runLme <- function (inData, inZ, inNumberOfOutputBriks, inContrasts, inAnovaIndices, inModel, inModelFormula, inRandomFormula) {
    ## cat("There will be " , inNumberOfOutputBriks, " stats briks\n")
    outStats <- vector(mode="numeric", length=inNumberOfOutputBriks)
    ##cat (inData, "\n")
    if ( ! all(inData == 0 ) ) {
        
        ## if inData is all zero then we are in a portion of the masked
        ## out data and should therefore not perform the lme on it just
        ## return a vector of zeros, if it's not then we get in this
        ## branch and can try the lme
        
        inModel$fmri<-inData
        
        ## If lme throws an exception it returns an object, assigned to
        ## mylme, which inherits from the class try-error. If mylme
        ## inherits from this class an error was generated by a particular
        ## voxel. It will default to having its corresponding voxels in
        ## the outStats array set to 0.
        
        if( inherits(
            mylme <- try(lme(fixed=inModelFormula, random=inRandomFormula, data = inModel),
                         silent=FALSE),
            "try-error") ) {
            temp <- 0
            cat (paste("Error on slice", inZ, "\n"))
        } else {
            temp <- as.vector(unlist(anova(mylme, type="marginal")))
        }
        
        if(length(temp) > 1) {
            numberOfContrasts=length(inContrasts)
            outStats[1:(inNumberOfOutputBriks-(numberOfContrasts*2))] = temp[c(inAnovaIndices)]
            
            if (numberOfContrasts > 0 ) {
                contrastsStartAt=(inNumberOfOutputBriks-(numberOfContrasts*2))+1
                for (i in seq(1, numberOfContrasts, by=1)) {
                    if (inherits(con <- try(contrast(mylme, a =inContrasts[[i]]$a, b =inContrasts[[i]]$b, type="average"), silent=FALSE), "try-error") ) {
                        con=0
                    }
                    ## con is of class "contrast" but is also of class "list" so this works as expected
                    if (is.list(con)) {
                        outStats[contrastsStartAt:(contrastsStartAt+1)] <- c(con$Contrast, con$testStat)
                    } else {
                        outStats[contrastsStartAt:(contrastsStartAt+1)] <- c(0, 0)
                    }
                    contrastsStartAt=contrastsStartAt+2
                } ## end of for (i in seq(1, numberOfContrasts, by=1))
            } ## end of if (numberOfContrasts > 0 )
        } ## end of if(is.list(temp))
    } ## end of if ( ! all(inData == 0 ) )

    ## else {
    ##   cat("Got all zeros :-(\n")
    ## }
    
    return(outStats)
}
## enable debugging of runLme
##debug(runLme)
##trace("runLme", quote(if(! all(inData == 0 ) ) browser()))
####################################################################################################

prePlannedContrasts=list()

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)

groups="mddAndCtrl"

seeds=readSeedsFile(file.path(config.data.dir, "Harvard-Oxford_amygdala_seeds.txt"))
numberOfSeeds=length(seeds)

seedCount=1
for (seed in seeds) {

    seedName=getSeedName(seed)
    
    cat("####################################################################################################\n")
    cat(sprintf("*** Running LME for the %s contrast (%02d of %02d)\n", seedName, seedCount, numberOfSeeds))

    ## setup all the filenames
    inputBucketFilename=paste("restingstate.bucket", groups, seedName, "masked+tlrc.HEAD", sep=".")
    outputBucketPrefix= paste("restingstate",        groups, seedName, "lme.bucket",       sep=".")
    
    ## this file stores the order of the subjects in each of the following BRIK files
    subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", groups, seedName, "csv", sep="."))
    subjectOrder=splitSubjectOrderIntoIdAndTimepoint(fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename)))
    
    inputBrikFilename=file.path(group.data.dir, inputBucketFilename)
        
    cat("*** Reading", inputBrikFilename, "\n")
    inputBrik=read.AFNI(inputBrikFilename, verb=verbose)
    
    clusterLogFilename=paste("restingstate.cluster", groups, seedName, "log", sep=".")

    ## dim is the number of parameters from the lme that will be stored as
    ## subbriks and exported back for use in AFNI
    dimX=inputBrik$dim[1]
    dimY=inputBrik$dim[2]
    dimZ=inputBrik$dim[3]
    
    ## this stored the total number of subbriks for all subjects for all types (i.e. wins, losses, and ties)
    numberOfBriks=inputBrik$dim[4]
    
    ## now merge the subject order and the demographics file so we have
    ## the table completed for each subbrik in in the bucket file
    
    ## the demographics file will contain only one row for each subject
    ## but the subjectOrder file contains one row for each subject for
    ## each period/block type in the experiment and hence will be much
    ## longer, this allows us to merge the two tables so we can complete
    ## the table for use in creation of the SPM. mgd stands for merged
    
    mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI")])
    ## ensure that subject is s factor
    mgd$subject=as.factor(mgd$subject)

    mgd=fixDates(mgd)
    mgd=computeAge(mgd)

    mrData = inputBrik$brk
    dim(mrData) = c(dimX, dimY, dimZ, numberOfBriks)
    
    ## now multiply by the mask to ensure that only those voxels in the
    ## mask will be included in the lme computation, all others are set to
    ## 0. The check for whether a voxel is in the mask of not is to be
    ## found in runLme at if ( ! all(inData == 0 ) )
    ##mrData = array(apply(mrData, 4, function (x) x * maskBrik$brk[,,,1]), dim=c(dimX, dimY, dimZ, numberOfBriks))

# modelFormula=as.formula("fmri ~ group * timepoint")
    modelFormula=as.formula("fmri ~ group + gender + age")
    modelFormula=as.formula("fmri ~ group")
    randomFormula = as.formula("random = ~ 1 | subject")
    # model = data.frame(
    #   "subject"   = as.vector(mgd$subject),
    #   "group"     = as.vector(mgd$Grp),
    #   "timepoint" = as.vector(mgd$timepoint)
    #)
    model = data.frame(
    "subject"   = as.vector(mgd$subject),
    "group"     = as.vector(mgd$Grp),
    "gender"    = as.vector(mgd$Gender),
    "age"       = as.vector(mgd$age.in.years)
    )
    
    model = data.frame(
    "subject"   = as.vector(mgd$subject),
    "group"     = as.vector(mgd$Grp)
    )
    
    contrs=prePlannedContrasts

    ## print(head(model))
    cat("*** Using the following formula for the model  effects:", capture.output(print(modelFormula)),  "\n")
    cat("*** Using the following formula for the random effects:", capture.output(print(randomFormula)), "\n")

    ##stop("Check the mgd data frame\n")
    
    nContrasts=length(contrs)
    
    ## if we got this far we should be able to run a single voxel to work
    ## out the DOF for the various stats, generate labels and indices for
    ## use with the unlisted anova in the runLme function
    i = 21
    j = 21
    k = 32
    model$fmri = mrData[i, j, k, ]

    stopIfFormulaTermsNotInModelMatrix(modelFormula, randomFormula, model)

    if( inherits(tempLme <- try(lme(fixed=modelFormula, random=randomFormula, data = model),
                                silent=FALSE),
                 "try-error") ) {
        tempAnova <- 0
        stop("Got an exception trying to setup the tempLme variable. Cannot continue beyond this point. Stopping.")
    } else {
        tempAnova <- anova(tempLme, type="marginal")
        cat("The temporary ANOVA is\n")
        print(tempAnova)
    }
    contr.df <- vector(mode="numeric", length=nContrasts)
    if (nContrasts > 0 ) {
        for (i in seq(1, nContrasts, by=1)) {
            con1=try(contrast(tempLme, a =contrs[[i]]$a, b =contrs[[i]]$b, type="average"))
            if (length(con1) > 1) {
                contr.df[i] = con1$df
            }
        }
    }

    ## stop("Did the tmeplme work?\n")  
    
    ## These are the labels that will be attached to the subbriks in the
    ## output stats bucket file. The labels will be dictated by the model
    ## formula and importantly the order in which the coefficients form
    ## the ANOVA matrix are concatenated
    baseStatsBrikLabels=makeBrikLabels(anova(tempLme))
    contrastsStatsBrikLabels=makeContrastBrikLabels(contrs)
    outputStatsBrikLabels=c(baseStatsBrikLabels, contrastsStatsBrikLabels)
    
    ## total number of briks in the output file, including everything
    ## from the anovas and contrasts
    numberOfOutputBriks=length(outputStatsBrikLabels)
    anovaIndices=makeAnovaIndices(anova(tempLme))
    
    ## ## outStats <- vector(mode="numeric", length=numberOfOutputBriks)
    ## ## if (is.list(tempAnova)) {
    ## ##     outStats[1:(numberOfOutputBriks-(nContrasts*2))] = as.vector(unlist(tempAnova))[c(anovaIndices)]
    ## ##     if (nContrasts > 0 ) {
    ## ##         contrastsStartAt=(numberOfOutputBriks-(nContrasts*2))+1
    ## ##         for (i in seq(1, nContrasts, by=1)) {
    ## ##             cat(sprintf("i=%d, contrastsStartAt=%d\n", i, contrastsStartAt))
    ## ##             if (inherits(con <- try(contrast(tempLme, a =contrs[[i]]$a, b =contrs[[i]]$b, type="average"), silent=FALSE), "try-error") ) {
    ## ##                 con=0
    ## ##             }
    ## ##             print(con)
    ## ##             ## con is of class "contrast" but is also of class "list" so this works as expected
    ## ##             if (is.list(con)) {
    ## ##                 outStats[contrastsStartAt:(contrastsStartAt+1)] <- c(con$Contrast, con$testStat)
    ## ##             } else {
    ## ##                 outStats[contrastsStartAt:(contrastsStartAt+1)] <- c(0, 0)
    ## ##             }
    ## ##             contrastsStartAt=contrastsStartAt+2
    ## ##         } ## end of for (i in seq(1, numberOfContrasts, by=1)) 
    ## ##     } ## end of if (numberOfContrasts > 0 )
    ## ## } ## end of if (is.list(tempAnova)) {

    
    ##  stop("Stopping. Everything ok so far?")

    
    Stats = array(0, c(dimX, dimY, dimZ, numberOfOutputBriks))
    cat(paste("Starting at", date(), "\n"))
    startTime=proc.time()
    if (ncpus > 1 ) {
        ## multiple cpus
        library(snow)
        
        ## outfile="" should result in slave output being printed on the
        ## controlling tty device, i.e. on the master
        cat("Making a cluster with" , ncpus, "cpus\n")
        cluster = makeCluster(ncpus, type = "SOCK", outfile=file.path(group.results.dir, clusterLogFilename))
        clusterEvalQ(cluster, library(nlme));
        clusterEvalQ(cluster, library(contrast))
        ## may need to add contrasts here later
        if (useProgressBar) {
            pb <- txtProgressBar(min = 0, max = dimZ, style = 3)
        }
        for ( kk in 32:32 ) {
            if (useProgressBar) {
                setTxtProgressBar(pb, kk)
            } else {
                cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
            }
            a=parApply(cluster, mrData[ , , kk, ],  c(1, 2), runLme, inZ = kk, inNumberOfOutputBriks=numberOfOutputBriks, inContrasts=contrs,
                inAnovaIndices=anovaIndices, inModel=model, inModelFormula=modelFormula, inRandomFormula=randomFormula)
            ##cat(dim(a), "\n")
            Stats[ , , kk, ] = aperm(a, c(2, 3, 1))
        }
        stopCluster(cluster)
    } else { ## serial / single CPU
        if (useProgressBar) {
            pb <- txtProgressBar(min = 0, max = dimZ, style = 3)
        }
        for ( kk in 1:dimZ ) {
            if (useProgressBar) {
                setTxtProgressBar(pb, kk)
            } else {
                cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
            }
            a=apply(mrData[ , , kk, ],  c(1, 2), runLme, inZ = kk, inNumberOfOutputBriks=numberOfOutputBriks, inContrasts=contrs,
                inAnovaIndices=anovaIndices, inModel=model, inModelFormula=modelFormula, inRandomFormula=randomFormula)
            ##cat(dim(a), "\n")
            Stats[ , , kk, ] = aperm(a, c(2, 3, 1))
        }
    } ## single cpu
    
    if (useProgressBar) {
        close(pb)
    }
    
    cat(paste("Ended at", date(), "\n"))
    cat("Time consumed\n")
    print(proc.time() - startTime)
    
    ##stop()
    
    lmeOutBrikFilename=paste(outputBucketPrefix, ".", format(Sys.time(), "%Y%m%d-%H%M%Z"), view.AFNI.name(inputBrikFilename), sep="")
    lmeOutBrikFqfn=file.path(group.results.dir, lmeOutBrikFilename)
    hostname=system('hostname', intern=T)
    user=Sys.getenv("USER")
    cat("*** Writing bucket file ", lmeOutBrikFqfn, "\n")
    write.AFNI(lmeOutBrikFqfn,
               Stats, verb=verbose,
               ##label = baselineBrik$head$DATASET_NAME,
               label=outputStatsBrikLabels,
               note = paste(
                   paste(paste("[", user, "@", hostname, ": ",  date(), "]", sep=""), file.path(script.location, script.name)),
                   paste("Model  formula:", gsub("~", "TILDE",  capture.output(print(modelFormula)))[1], sep=" "),
                   paste("Random formula:", gsub("~", "TILDE ", capture.output(print(randomFormula)))[1], sep=" "), sep="\n"),
               origin = inputBrik$origin,
               delta = inputBrik$delta,
               orient= inputBrik$orient)
    
    
    statpar = "3drefit"
    if ( is.list(tempAnova) ) {
        cat ("Making statpar arguments\n")
        
        afniFtestBrikIds=makeAfniFtestBrikIds(tempAnova)
        afniTtestBrikIds=makeAfniTtestBrikIds(tempAnova, contrs)
        statparArguments=makeAfniStatparArguments(tempAnova, afniFtestBrikIds, afniTtestBrikIds, contr.df)
        statpar = paste(statpar, statparArguments)    
    }
    statpar = paste("(cd",  group.results.dir, ";", statpar, " -view tlrc -space MNI -addFDR -newid ", lmeOutBrikFilename, ")")
    cat(statpar, "\n")
    system(statpar)
    seedCount=seedCount + 1
} ## end of  for (seed in seeds)
