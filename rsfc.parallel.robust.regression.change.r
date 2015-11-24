#!/usr/bin/Rscript

rm(list=ls())

library(MASS)
library(getopt)
## library(compiler)

source("scoreMasc.r")

AFNI_R_DIR=Sys.getenv("AFNI_R_DIR", unset=NA)

## use the functions for loading and saving briks from the AFNI
## distribution as they can cleanly handle floats/doubles
if ( ! is.na(AFNI_R_DIR) ) {
    source(file.path(AFNI_R_DIR, "AFNIio.R"))
} else {
    stop("Couldn't find AFNI_R_DIR in environment. This points to the location from which to load functions for reading and writing AFNI BRIKS. Stopping!")
}

##########################################################################################################################################################################
### START OF FUNCTIONS ###################################################################################################################################################
##########################################################################################################################################################################


## this function takes the row and column names from an rlm
## coefficient matrix and combines these (along row) to get a list of
## labels for later use as labels to the 3drefit command so the briks
## in the bucket are correctly labeled. Spaces are replaced with . and
## ( or ) are deleted.
makeBrikLabels <- function (inRlmCoef, inBoot=FALSE) {
    rns=rownames(inRlmCoef)
    cns=colnames(inRlmCoef)
    
    labels=c()
    
    for (i in 1:length(rns) ) {
        for (j in 1:length(cns) ) {
            name=gsub("..", ".", gsub(" ", ".", paste(gsub("[(]|[)]", "", rns[i]), cns[j], sep=".")), fixed=TRUE)
            labels=c(labels, name)
        }
    }
    if (inBoot) {
        for (i in 1:length(rns) ) {
            for (j in 1:length(bootstrappingBrikLabelSuffxes)) {
                name=gsub("..", ".", gsub(" ", ".", paste(gsub("[(]|[)]", "", rns[i]), bootstrappingBrikLabelSuffxes[j], sep=".")), fixed=TRUE)
                labels=c(labels, name)
            }
        }
    }

    outList=list("numberOfLabels"=length(labels),
        "bootstrapLabelsStartAt"=ifelse(inBoot, (length(rns)*length(cns)) + 1, NA),
        "labels"=labels)

    return (outList)
}

## make a list of the AFNI brik indices that correspond to the t-stat
## from the matrix of coefficients produced by rlm. This take into
## account that AFN inumbers from 0 not 1
makeAfniTtestBrikIds <- function(inRlmCoef, inBoot=FALSE) {
    nrows=length(rownames(inRlmCoef))
    ncols=length(colnames(inRlmCoef))
    
    indices=seq(3, nrows*ncols, by=3)

    if (inBoot) {
        indices=c(indices, seq(indices[length(indices)]+2, (nrows*ncols) + (nrows*length(bootstrappingBrikLabelSuffxes)), by=4))
    }
    
    return (indices)
}

## given a list of AFNI brik indices and a degrees of freedom (inDf)
## make an appropriate list of statpar arguments to add to a 3drefit
## command line
makeAfniStatparArguments <- function(inDf, inAfniFtestBrikIds) {
    return( paste(sapply(inAfniFtestBrikIds, function(x) { sprintf("-substatpar %d fitt %d", x[1], inDf) }), collapse=" ") )
}

cleanRegressionVariable <- function(inName) {
    ## Replace 2 or more consequtive . with one . and remove any
    ## trailing . from the name using the gsub function
    return (gsub("\\.{2,}", ".", gsub("\\.$", "", inName)))
}

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

printExcludedAndIncludedSubjects <- function(inData) {

    variable.columns=-c(grep("Grp|subject", colnames(inData), fixed=FALSE))
    nas = any(is.na(inData[, variable.columns]))
    not.nas = ! nas
    if (nas) {
        cat (sprintf("*** The following subjects have missing data for %s\n", rvVariable))
        cat(c("---", as.vector(inData[nas, "subject"])), "\n")
        cat(c("---", as.vector(inData[nas, "Grp"])), "\n")
        ## cat("*** Total number of subjects to be excluded/included in the regression analysis:\n")
        ## print(addmargins(table(ifelse(is.na(inData[, rvVariable]), "Exclude", "Inlcude"))))
    }
    cat (sprintf("*** The following subjects will be included in the regression analysis\n"))
    cat(c("***", as.vector(inData[not.nas, "subject"])), "\n")
    cat(c("***", as.vector(inData[not.nas, "Grp"])), "\n")
}

## this is a very trimmed down version of the runRegression function
## below. It is only for use with bootstrapping
bootRegression <- function(inData, inIndices, inModelFormula, inMaxIt=50, inNumberOfStatsBriks) {
    outStats <- vector(mode="numeric", length=inNumberOfStatsBriks)
    if ( ! inherits(myrlm <- try(rlm(inModelFormula, data=inData[inIndices, ], maxit=inMaxIt), silent=FALSE),
                    "try-error") ) {
        ## cat ("length(coefficients(myrlm)) is ", coefficients(myrlm), "\n")
        return(coefficients(myrlm))
    } else {
        cat("Got an exception\n")
        return(outStats)
    }
}

## bootRegression <- cmpfun(bootRegression.slow)

runRegression <- function (inData, inNumberOfStatsBriks, inModel, inModelFormula, inMaxIt=50, inBoot=FALSE, inR=25, inBootstrapStatsStartAt=NA) {
    outStats <- vector(mode="numeric", length=inNumberOfStatsBriks)
    if ( ! all(inData == 0 ) ) {

        ## if inData is all zero then we are in a portion of the masked
        ## out data and should therefore not perform the rlm on it just
        ## return a vector of zeros, if it's not then we get in this
        ## branch and can try the rlm
        
        inModel$mri<-inData
        
        myrlm <- rlm(inModelFormula, data = inModel, maxit=inMaxIt)
        
        outStats[1] = mean(inModel$mri)
        ## print(summary(myrlm))
        ## > coef(summary(model))
        ## Value Std. Error    t value
        ## (Intercept) 89.09177404  8.0221830 11.1056770
        ## Age          0.04087267  0.1146234  0.3565821
        ## Educ         0.71157766  0.5370092  1.3250753
        ## as.vector(t(coef(summary(model))))
        ## [1] 89.09177404  8.02218304 11.10567705  0.04087267  0.11462344  0.35658210  0.71157766  0.53700922  1.32507531
        ##outStats[2:length(outStats)]=as.vector(t(coef(summary(myrlm))))

        if (inBoot) {
            if (is.na(inBootstrapStatsStartAt)) {
                stop("***ERROR in runRegression: inBootstrapStatsStartAt has not been set. It is currently NA. Cannot continue. Stopping\n")
            }

            ## cat ("\nRLM Indices: ", 2:(inBootstrapStatsStartAt-1), "\n")
            ## cat ("coef(summary(myrlm)): ", as.vector(t(coef(summary(myrlm)))), "\n")
            ## cat ("Length of stats vector is", length(as.vector(t(coef(summary(myrlm))))), "\n")
            ## cat ("inBootstrapStatsStartAt is ", inBootstrapStatsStartAt, "\n")

            outStats[2:(inBootstrapStatsStartAt-1)]=as.vector(t(coef(summary(myrlm))))

            ## cat ("inside inBoot outStats is: ", outStats, "\n")
            ## cat ("number of stats briks should be: ", length(2:(inBootstrapStatsStartAt-1)), "\n")
            
            bootStats=boot(inModel, bootRegression, R=inR, inModelFormula=inModelFormula, inMaxIt=inMaxIt, inNumberOfStatsBriks=length(2:(inBootstrapStatsStartAt-1)))
            ## cat("bootStats is\n")
            ## print(bootStats)
            ## print(class(bootStats))
            ## print(is.vector(bootStats))
            if (is(bootStats, "boot")) {
                bootBias=apply(bootStats$t, 2, mean) - bootStats$t0
                ## these are vectors with as many columns as there are terms in
                ## the regression model. Don't forget to include the intercept
                ## when you're trying to mentalize this
                bootSeCoeff=apply(bootStats$t, 2, sd)
                bootTCoeff=bootStats$t0 / bootSeCoeff
                ## bootCi=matrix(0, nrow=length(bootTCoeff), ncol=2)

                ## cat("1 bootBias is",    bootBias, "\n")
                ## cat("2 bootSdCoeff is", bootSdCoeff, "\n")
                ## cat("3 bootTCoeff is",  bootTCoeff, "\n")
                ## cat("4 bootCi is:",     bootCi, "\n")                                
                ## Outputbootstrapstats=vector(mode="numeric", length=length(bootTCoeff) * length(bootstrappingBrikLabelSuffxes))
                outputBootstrapStats=c()
                for (termIndex in seq(1, length(bootTCoeff))) {
                    ## the normal element of the CI value contains 3 elements: 1) the CI
                    ## level (in this case 0.95), 2) the lower bound on the CI, 3) the
                    ## upper bound on the CI
                    ci=boot.ci(bootStats, conf = c(0.95), type = c("norm"), index = termIndex)$normal[2:3]
                    ## cat("5 ci is:", ci, "\n")

                    ## bootCi[termIndex, ]=ci

                    ## cat("6 bootCi is", bootCi, "\n")                                                    
                    ## outputBootstrapStats=c(outputBootstrapStats, bootBias[termIndex], bootTCoeff[termIndex], bootCi[termIndex, ])
                    outputBootstrapStats=c(outputBootstrapStats, bootBias[termIndex], bootTCoeff[termIndex], ci)                    
                    ## cat ("7 outputBootstrapStats now: ", outputBootstrapStats, "\n")
                }
                outStats[inBootstrapStatsStartAt:inNumberOfStatsBriks]=outputBootstrapStats
            }
        } ## end of if (inBoot) {
        else {
            outStats[2:length(outStats)]=as.vector(t(coef(summary(myrlm))))
        }
    } ## end of if ( ! all(inData == 0 ) ) {
    ## cat ("8 outStats is: ", outStats, "\n")  
    ## if ( ! all(outStats == 0 ) )
    ##     stop()

    return(outStats)
}

## runRegression <- cmpfun(runRegression.slow)

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}

help <- function(){

}

checkCommandLineArguments <- function (in.opt) {
    ## if help was asked for print a friendly message
    ## and exit with a non-zero error code
    if ( !is.null(in.opt$help) ) {
        cat(getopt(spec, usage=TRUE));
        q(status=1);
    }

    if (is.null(in.opt$change)) {
        cat("A change file is required\n")
        q(status=1);
    }

    if (substr(in.opt$change, 1, 1) == "." || substr(in.opt$change, 1, 1) == "/") {
        cat ("The change file cannot begin with a / or .\n",
             "It must not be either an absolute of relative pathname.\n",
             "Rather it must be the name of a file that is assumed to live in the data/admin folder.\n", sep="")
        q(status=1)
    }

    if ( is.null(in.opt$bootstrap)) {
        in.opt$bootstrap=FALSE
    }

    if (in.opt$bootstrap) {
        ## load the boot strapping library only if bootstrapping is
        ## requested
        library(boot)
    }
    
    if (is.null(in.opt$variable)) {
        cat("A variable name is required.\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1);
    }

    if (is.null(in.opt$seeds)) {
        cat("A file name containing the list of seeds to use is required.\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1);
    }

    if (is.null(in.opt$resamples)) {
        in.opt$resamples=25
    }

    if ( is.null(in.opt$progress)) {
        in.opt$progress=FALSE
    }

    if ( is.null(in.opt$cleaned)) {
        in.opt$cleaned=FALSE
    }

    return(in.opt)
}

printOptionsSummary <- function () {
    cat("*** Change file is set to:", change.score.filename, "\n")
    if (opt$bootstrap) {
        cat("*** Performing bootstrapping of the regression using", opt$resamples, "resamples\n")
    } else {
        cat("*** Bootstrapping will NOT be performed\n")
    }
    if (opt$cleaned) {
        cat("*** Using cleaned datasets\n")
    }
    cat("*** The variable name is:", opt$variable, "\n")
    cat("*** Seeds will be read from:", seeds.filename, "\n")
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

## enable debugging of runRlm
##debug(runRegression)
##trace("runRegression", quote(if(! all(inData == 0 ) ) browser()))
####################################################################################################


## set.seed(100)

## print(Sys.info())
## if (! is.na(Sys.info()["NSLOTS"] )) {
##     ncpus=as.integer(Sys.info()["NSLOTS"])
##     cat(paste("*** Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "allocated by SGE\n"))        
## } else

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
    ncpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
    cat(paste("*** Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
    ncpus=8
    cat(paste("*** Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))    
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

## ncpus=1

scripts.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts")
data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
group.data.dir=file.path(data.dir, "Group.data")
group.results.dir=file.path(data.dir, "Group.results")
seeds.data.dir=file.path(data.dir, "seeds")

## should the AFNI read/write routines produce verobse output
verbose=FALSE

## this indicates whether bootstrapping should be performed
## doBootstrapping=FALSE

bootstrappingBrikLabelSuffxes=c("bias", "t.value", "ciLower", "ciUpper")

################################################################################
NO_ARGUMENT="0"
REQUIRED_ARGUMENT="1"
OPTIONAL_ARGUMENT="2"

## process command line arguments
spec = matrix(c(
    'help',          'h', NO_ARGUMENT,       "logical",
    'change',        'c', REQUIRED_ARGUMENT, "character",
    'bootstrap',     'b', NO_ARGUMENT,       "logical",
    "variable",      'v', REQUIRED_ARGUMENT, "character",
    "seeds",         's', REQUIRED_ARGUMENT, "character",
    "progress",      'p', NO_ARGUMENT,       "logical",
    "resamples",     'r', REQUIRED_ARGUMENT, "integer",
    "cleaned",       'e', NO_ARGUMENT,       "logical"
), byrow=TRUE, ncol=4)

if (interactive()) {
    cat("*** Setting interactive options\n")
    ## these are default arguments that are useful for testing
    ## purposes.
    ##CDRS.t.score.scaled.diff.change.score.csv 
    ## args=c(
    ##     ## "-b", "-r", "5000", ## "-e",
    ##     "-c", "new.mdd.CDRS.t.score.scaled.diff.change.score.csv",
    ##     "-v", "CDRS.t.score.scaled.diff",
    ##     "-s", "juelich_amygdala_seeds_weights.txt")
    ## opt = getopt(spec, opt=args)
} else {
    opt = getopt(spec)
}

opt=checkCommandLineArguments(opt)

## flag to indincate whether a progress bar should be printed or one
## line per slice if no progress bar is to be used
useProgressBar=interactive() || opt$progress

change.score.filename=file.path(admin.data.dir, opt$change)
seeds.filename=file.path(config.data.dir, opt$seeds)

printOptionsSummary()

## stop("Check the output of the getopt processing\n")

change=readCsvFile(change.score.filename, "ID")
groups="mddOnly"

seeds=readSeedsFile(seeds.filename)

numberOfSeeds=length(seeds)

totalNumberOfRegressions=numberOfSeeds

regressionCount=1

rvVariable=opt$variable

if (grepl("diff|rstandard", rvVariable, fixed=FALSE)) {
    ## the temporal difference regressions have seperate directories
    group.data.dir=file.path(data.dir, paste("Group.data", rvVariable, "withAandC", sep="."))
    group.results.dir=file.path(data.dir, paste("Group.results", rvVariable, "withAandC", sep="."))
} else if (grepl("both", rvVariable, fixed=TRUE) ) {
    group.data.dir=file.path(data.dir, paste("Group.data", rvVariable, sep="."))
    group.results.dir=file.path(data.dir, paste("Group.results", rvVariable, sep="."))
} else {
    ## the temporal difference regressions have seperate directories
    group.data.dir=file.path(data.dir, paste("Group.data", rvVariable, "withAandC", sep="."))
    group.results.dir=file.path(data.dir, paste("Group.results", rvVariable, "withAandC", sep="."))
}

if (opt$cleaned) {
    group.data.dir=file.path(data.dir, paste("Group.data", rvVariable, "withAandC", "cleaned", sep="."))
    group.results.dir=file.path(data.dir, paste("Group.results", rvVariable, "withAandC", "cleaned", sep="."))
}

if ( ! dir.exists( group.results.dir) ) {
    stop(sprintf("*** The group results directory (%s) does not exist. Cannot continue\n", group.results.dir))
}
## turn on compilation of all loops before they are run for the first time
## enableJIT(3)

for (seed in seeds) {
    seedName=getSeedName(seed)
    
    cat("####################################################################################################\n")
    cat(sprintf("*** Running RLM for the %s seed of the %s group (%02d of %02d regressions)\n", seedName, groups, regressionCount, totalNumberOfRegressions))
    
    ## setup all the filenames
    inputBucketFilename=paste("restingstate.bucket",                  groups, seedName,             "masked+tlrc.HEAD",  sep=".")
    if (opt$bootstrap) {
        outputBucketPrefix= paste("restingstate.reversed.formula",    groups, seedName, rvVariable, "booted.rlm.bucket", sep=".")
    } else {
        outputBucketPrefix= paste("restingstate.reversed.formula",    groups, seedName, rvVariable, "rlm.bucket",        sep=".")
    }

    clusterLogFilename=paste("restingstate.cluster", groups, seedName, "log", sep=".")
    
    ## this file stores the order of the subjects in each of the following BRIK files
    subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", groups, seedName, "csv", sep="."))
    subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))
    
    inputBrikFilename=file.path(group.data.dir, inputBucketFilename)
    
    cat("*** Reading", inputBrikFilename, "\n")
    inputBrik=read.AFNI(inputBrikFilename, verb=verbose)

    ## dim is the number of parameters from the rlm that will be stored as
    ## subbriks and exported back for use in AFNI
    dimX=inputBrik$dim[1]
    dimY=inputBrik$dim[2]
    dimZ=inputBrik$dim[3]
    
    ## this stored the total number of subbriks for all subjects for all types (i.e. wins, losses, and ties)
    numberOfBriks=inputBrik$dim[4]
    
    mgd=merge(subjectOrder, change, by.x="subject", by.y="ID", sort=FALSE)

    if (! all ( mgd$subject == subjectOrder$subject) ) {
        print(mgd$subject)
        print(subjectOrder$subject)
        print(change)
        print(dim(mgd))
        print(dim(subjectOrder))
        print(dim(change))
        stop("*** The order of the subjects in the mgd data frame is not the same as the order of subjects in subjecOrder.\n",
             "*** Cannot continue.\n")
    }
    
    if (dim(mgd)[1] != dim(subjectOrder)[1] ) {
        cat("*** The number of subjects in the merged data frame is not the same as the number of subjects in the subjectOrder file.\n")
        cat("*** Cannot continue\n")
        q(status=1)
    }

    mgd$subject=as.factor(mgd$subject)
    mgd=droplevels(mgd)
    cat("mgd:\n")
    print(head(mgd))
    
    mrData = inputBrik$brk
    dim(mrData) = c(dimX, dimY, dimZ, numberOfBriks)
    
    ## stop("Check mgd data frame\n")
    
    ## set up the data frame with the dependant variables we will want
    ## to include in the regression formula

    cat("*** Setting up model data frame\n")
    if (grepl("both.scaled", rvVariable, fixed=TRUE)) {
        formula.variable=sub(".both.scaled", "", rvVariable)

        print( paste(formula.variable, c("A.scaled", "C.scaled"), sep="."))
        model = data.frame(
            as.vector(mgd[, "Grp"]),
            as.vector(mgd[, "subject"]),            
            as.vector(mgd[, c("age.in.years.scaled", paste(formula.variable, c("A.scaled", "C.scaled"), sep="."))]))
        colnames(model) = c("Grp", "subject", "age.in.years.scaled", paste(formula.variable, c("A.scaled", "C.scaled"), sep="."))
    } else if (grepl("both", rvVariable, fixed=TRUE)) {
        formula.variable=sub(".both", "", rvVariable)
        model = data.frame(
            as.vector(mgd[, "Grp"]),
            as.vector(mgd[, "subject"]),            
            as.vector(mgd[, c("age.in.years", paste(formula.variable, c("A", "C"), sep="."))]))
        colnames(model) = c("Grp", "subject", "age.in.years", paste(formula.variable, c("A", "C"), sep="."))
    } else {
        formula.variable=rvVariable
        model = data.frame(
            as.vector(mgd[, "Grp"]),
            as.vector(mgd[, "subject"]),            
            as.vector(mgd[, rvVariable]),
            as.vector(mgd[, "age.in.years"]))
        colnames(model) = c("Grp", "subject", paste(formula.variable, c("A.scaled", "C.scaled"), sep="."), "age.in.years")
    }

    printExcludedAndIncludedSubjects(mgd)
    print(head(model))
    
    ## stop("Check model data frame\n")
    
    ## set up the formula for the model that we want to regress
    ## modelFormula  = as.formula(paste("mri ~", rvVariable, "+ age.in.years", sep=" "))
    cat("*** Setting up model formula\n")
    if (grepl("both.scaled", rvVariable, fixed=TRUE)) {
        modelFormula  = as.formula(
            paste(paste(formula.variable, "C.scaled", sep="."), "~ mri + ", paste(formula.variable, "A.scaled", sep="."), "+ age.in.years.scaled", sep=" "))
    } else if (grepl("both", rvVariable, fixed=TRUE)) {
        modelFormula  = as.formula(
            paste(paste(formula.variable, "C", sep="."), "~ mri + ", paste(formula.variable, "A", sep="."), "+ age.in.years", sep=" "))
    } else {
        modelFormula  = as.formula(paste(formula.variable, "~", "mri + age.in.years", sep=" "))
    }

    cat("*** The model formula is: ")
    print(modelFormula)
    
    ## if we got this far we should be able to run a single voxel to work
    ## out the DOF for the various stats, generate labels and indices for
    ## use with the unlisted rlm in the runRlm function
    i = 21
    j = 21
    k = 32
    model$mri = mrData[i, j, k, ]
    if( inherits(tempRlm <- try(rlm(modelFormula, data = model),
                                silent=FALSE),
                 "try-error") ) {
        tempRlmCoef=0
        tempRlmDf=0
        traceback()
        stop("Got an exception trying to setup the tempRlmCoef and tempRlmDf variables. Cannot continue beyond this point. Stopping.")
    } else {
        s=summary(tempRlm)
        cat("*** The temporary RLM summary is:\n")
        print(s)
        tempRlmCoef = coef(s)
        tempRlmDf=s$df
    }
    
    ## stop("Check RLM summary\n")
    
    ## These are the labels that will be attached to the subbriks in the
    ## output stats bucket file. The labels will be dictated by the
    ## model formula and importantly the order in which the coefficients
    ## from the RLM coefficient matrix are concatenated
    ## the makeBrikLabels
    ## function does not include the mean of the fMRI at the start of
    ## the list of labels so include it here with the c()
    d=makeBrikLabels(tempRlmCoef, inBoot=opt$bootstrap)
    outputStatsBrikLabels=c("Mean", d[["labels"]])
    
    ## we add 1 to both numberOfStatsBriks and numberOfStatsBriks
    ## becasue makeBrikLabels does not take into account the
    ## factthat we will add an additional subbrik (the mean) to
    ## the output stats. Hence the output of makeBrikLabels is
    ## always 1 too small
    
    ## the number of stats subbriks to write out. This is dictated by the
    ## model Formula, changes to it likely imply changes to this number
    ##numberOfStatsBriks=length(outputStatsBrikLabels)
    numberOfStatsBriks=d[["numberOfLabels"]] + 1
    
    bootstrapStatsStartAt=0
    if (opt$bootstrap) {
        bootstrapStatsStartAt=d[["bootstrapLabelsStartAt"]] + 1
    }
    
    ##stop("Stopping")
    
    maxIter=100
    resamples=opt$resamples
    Stats = array(0, c(dimX, dimY, dimZ, numberOfStatsBriks))
    cat(paste("Starting at", date(), "\n"))
    startTime=proc.time()
    if (useProgressBar) {
        pb <- txtProgressBar(min = 0, max = dimZ, style = 3)
    }
    if (ncpus > 1 ) {
        ## multiple cpus
        library(snow)
        
        ## cluster = makeCluster(ncpus, type = "SOCK")
        cluster = makeCluster(ncpus, type = "SOCK", outfile=file.path(group.results.dir, clusterLogFilename))
        clusterEvalQ(cluster, library(MASS))
        clusterEvalQ(cluster, library(boot))
        ## clusterEvalQ(cluster, library(compiler))        
        ## clusterEvalQ(cluster, enableJIT(3))     
        ## export the bootRegression function so that it can be sees
        ## on slave nodes in the cluster
        ## clusterExport(cluster, c("bootRegression", "bootRegression.slow"))
        
        ##function (inData, inNumberOfStatsBriks, inModel, inModelFormula, inMaxIt=50, inBoot=FALSE, inR=25, inBootstrapStatsStartAt=NA) {
        for ( kk in 1:dimZ) {
            if (useProgressBar) {
                setTxtProgressBar(pb, kk)
            } else {
                cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
            }
            Stats[ , , kk, ] = aperm(parApply(cluster, mrData[ , , kk, ],  c(1, 2), runRegression,
                     inNumberOfStatsBriks=numberOfStatsBriks, inModel=model, inModelFormula=modelFormula, inMaxIt=maxIter,
                     inBoot=opt$bootstrap, inR=resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
        }
        stopCluster(cluster)
    } else {
        for ( kk in 1:dimZ ) {
        ## for ( kk in 32:32 ) {            
            ## single cpu
            if (useProgressBar) {
                setTxtProgressBar(pb, kk)
            } else {
                cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
            }
            Stats[ , , kk, ] = aperm(apply(mrData[ , , kk, ],  c(1, 2), runRegression,
                     inNumberOfStatsBriks=numberOfStatsBriks, inModel=model, inModelFormula=modelFormula, inMaxIt=maxIter,
                     inBoot=opt$bootstrap, inR=resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
        }
    }
    
    if (useProgressBar) {
        close(pb)
    }
    
    cat(paste("Ended at", date(), "\n"))
    cat("Time consumed\n")
    print(proc.time() - startTime)
    
    ##stop()
    
    rlmOutBrikFilename=paste(outputBucketPrefix, ".", format(Sys.time(), "%Y%m%d-%H%M%Z"), view.AFNI.name(inputBrikFilename), sep="")
    rlmOutBrikFqfn=file.path(group.results.dir, rlmOutBrikFilename)
    hostname=system('hostname', intern=T)
    user=Sys.getenv("USER")
    cat("*** Writing bucket file ", rlmOutBrikFqfn, "\n")
    
    write.AFNI(rlmOutBrikFqfn,
               Stats, verb=verbose,
               ##label = baselineBrik$head$DATASET_NAME,
               label=outputStatsBrikLabels,
               note = paste(
                   paste(paste("[", user, "@", hostname, ": ",  date(), "]", sep=""), normalizePath(get_Rscript_filename())),
                   paste("Model  formula:", gsub("~", "TILDE ", capture.output(print(modelFormula)))[1], sep=" ")),
               defhead=inputBrik)
               ## origin = inputBrik$origin,
               ## delta = inputBrik$delta,
               ## orient= inputBrik$orient)
    
    statpar = "3drefit"
    
    if ( is.matrix(tempRlmCoef) ) {
        cat ("Making statpar arguments\n")
        
        afniTtestBrikIds=makeAfniTtestBrikIds(tempRlmCoef)
        statparArguments=makeAfniStatparArguments(tempRlmDf[2], afniTtestBrikIds)
        statpar = paste(statpar, statparArguments)    
    }
    statpar = paste("(cd",  group.results.dir, ";", statpar, " -view tlrc -space MNI -addFDR -newid ", rlmOutBrikFilename, ")")
    cat(statpar, "\n")
    system(statpar)
    regressionCount=regressionCount+1
} ## end of  for (seed in seeds)


