rm(list=ls())

##the name of this script
script.name=parent.frame(2)$ofile
## the location (absolute path) to this script
script.location=normalizePath(dirname(parent.frame(2)$ofile))

library(MASS)
library(getopt)

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

## Reads the seed file and does the (crude) equivalent of BAS variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character())
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

computeMascScore <- function (inData) {
    inData.dim=dim(inData)
    for (r in seq(1, inData.dim[1]) ) {
        ## cat("##################################################\n")
        subjectNumber=inData[r, "subject"]
        gender=inData[r, "Gender"]
        age=round(inData[r, "age.in.years"], 0)
        old.masc.tscore=inData[r, "MASC.tscore"]

        new.masc.tscore=scoreMasc(gender, age, inData[r, "MASC.total"])
        if (is.na(new.masc.tscore) ) {
            warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, inData[r, "MASC.total"]))
        }
        
        inData[r, "MASC.tscore"]=new.masc.tscore

        ## cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f new MASC tscore=%0.0f\n",
        ##             r, subjectNumber, gender, age, inData[r, "MASC.total"], old.masc.tscore, new.masc.tscore))

    }
    return (inData)
}


printExcludedAndIncludedSubjects <- function(inData) {
    if (any(is.na(inData))) {
        cat (sprintf("*** The following subjects have missing data for %s (%s)\n", rvName, rvVariable))
        cat(as.vector(inData[is.na(inData[, rvVariable]), "subject"]), "\n")
        cat(as.vector(inData[is.na(inData[, rvVariable]), "Grp"]), "\n")
        cat("*** Total number of subjects to be excluded/included in the regression analysis:\n")
        print(addmargins(table(ifelse(is.na(inData[, rvVariable]), "Exclude", "Inlcude"))))
    }
    cat (sprintf("*** The following subjects will be included in the regression analysis\n"))
    cat(as.vector(inData[! is.na(inData[, rvVariable]), "subject"]), "\n")
    cat(as.vector(inData[! is.na(inData[, rvVariable]), "Grp"]), "\n")
}

## this is a very trimmed down version of the runRegression function
## below. It is only for use with bootstrapping
bootRegression <- function(inData, inIndices, inModelFormula, inMaxIt=50, inNumberOfStatsBriks) {
    outStats <- vector(mode="numeric", length=inNumberOfStatsBriks)
    if ( ! inherits(myrlm <- try(rlm(inModelFormula, data=inData[inIndices, ], maxit=inMaxIt), silent=FALSE),
                    "try-error") ) {
        cat ("length(coefficients(myrlm)) is ", coefficients(myrlm), "\n")
        return(coefficients(myrlm))
    } else {
        cat("Got an exception\n")
        return(outStats)
    }
}

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
            ## outStats[2:(inBootstrapStatsStartAt-1)]=as.vector(t(coef(summary(myrlm))))
            ## cat ("inside inBoot outStats is: ", outStats, "\n")
            ## cat ("number of stats briks should be: ", length(2:(inBootstrapStatsStartAt-1)), "\n")
            
            bootStats=boot(inModel, bootRegression, R=inR, inModelFormula=inModelFormula, inMaxIt=inMaxIt, inNumberOfStatsBriks=length(2:(inBootstrapStatsStartAt-1)))
            if (is.vector(bootStats)) {
                bootBias=apply(bootStats$t, 2, mean) - bootStats$t0
                ## these are vectors with as many columns as there are terms in
                ## the regression model. Don't forget to include the intercept
                ## when you're trying to mentalize this
                bootSdCoeff=apply(bootStats$t, 2, sd)
                bootTCoeff=bootStats$t0 / bootSdCoeff
                bootCi=matrix(0, nrow=length(bootTCoeff), ncol=2)
                outputBootstrapStats=c()
                for (termIndex in seq(1, length(bootTCoeff))) {
                    ## the normal element of the CI value contains 3 elements: 1) the CI
                    ## level (in this case 0.95), 2) the lower bound on the CI, 3) the
                    ## upper bound on the CI
                    ci=boot.ci(bootStats, conf = c(0.95), type = c("norm"), index = termIndex)$normal[2:3]
                    bootCi[termIndex, ]=ci
                    outputBootstrapStats=c(outputBootstrapStats, bootBias[termIndex], bootTCoeff[termIndex], bootCi[termIndex, ])
                    ##cat ("outputBootstrapStats now: ", outputBootstrapStats, "\n")
                }
                outStats[inBootstrapStatsStartAt:inNumberOfStatsBriks]=outputBootstrapStats
            }
        } ## end of if (inBoot) {
        else {
            outStats[2:length(outStats)]=as.vector(t(coef(summary(myrlm))))
        }
    } ## end of if ( ! all(inData == 0 ) ) {
    ##cat ("outStats is: ", outStats, "\n")  
    return(outStats)
}


readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

## enable debugging of runRlm
##debug(runRegression)
##trace("runRegression", quote(if(! all(inData == 0 ) ) browser()))
####################################################################################################

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
group.data.dir=normalizePath(file.path(data.dir, "Group.data"))
group.results.dir=normalizePath(file.path(data.dir, "Group.results"))
seeds.data.dir=normalizePath(file.path(data.dir, "seeds"))

## flag to indincate whether a progress bar should be printed or one
## line per slice if no progress bar is to be used
useProgressBar=TRUE

## should the AFNI read/write routines produce verobse output
verbose=FALSE

## this indicates whether bootstrapping should be performed
doBootstrapping=FALSE

bootstrappingBrikLabelSuffxes=c("bias", "t.value", "ciLower", "ciUpper")

################################################################################

## this file stores the demographic/neuropsych/neuromed information
## for each subject. These two are later merged (the mgd variable)
## so that only some of the the info from demogrpahics (only from
## those subjects included in the subjectOrder file, which should
## match exactly the order of the inputBrik) is added to the data
## from subjectOrder

##demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)

cdrsr.change.score.filename=file.path(admin.data.dir, "cdrsr.change.csv")
cdrsr.change=readCsvFile(cdrsr.change.score.filename, "SubjID")

masc.tscore.change.filename=file.path(admin.data.dir, "MASC.tscore.change.csv")
masc.change=readCsvFile(masc.tscore.change.filename, "SubjNum")

bdi.change.filename=file.path(admin.data.dir, "BDI.change.csv")
bdi.change=readCsvFile(bdi.change.filename, "SubjNum")

cgas.change.filename=file.path(admin.data.dir, "CGAS.change.csv")
cgas.change=readCsvFile(cgas.change.filename, "SubjNum")

cdi.change.filename=file.path(admin.data.dir, "CDI.change.csv")
cdi.change=readCsvFile(cdi.change.filename, "SubjNum")

rads.change.filename=file.path(admin.data.dir, "RADS.Total.Tscore.change.csv")
rads.change=readCsvFile(rads.change.filename, "SubjNum")

regressionVariables=list(
    list(variable="CDRSR.diff",            name="Children's Depression Rating Scale (A to C Change)"),
    list(variable="MASC.tscore.diff",      name="Multidimensional Anxiety Scale for Children (Standardized; A to C Change)"),
    ## list(variable="BDI.diff",             name="Beck Depression Inventory II (A to C Change)"),
    list(variable="CGAS.diff",             name="Children's Global Assessment Scale"),
    ## list(variable="CDI.diff",                  name="Children's Depression Inventory (A to C Change)"),
    list(variable="RADS.Total.Tscore.diff", name="Reynolds Adolescent Depression Scale Total (Standardized; A to C Change)")
    
    ##list(variable="MASC.tscore",    name="Multidimensional Anxiety Scale for Children (Standardized)")
    #list(variable="CDRS.tscore",    name="Children's Depression Rating Scale (Standardized)")
    ##list(variable="BDI.II",         name="Beck Depression Inventory II")
    ## list(variable="RADS.Total.Tscore",    name="Reynolds Adolescent Depression Scale Total (Standardized)")
    )

## extract the list of variable names from the regressionVariables list rvs=regression variables
rvs=unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))]
## select only the columns we want to perform regressions on from the demographics data frame
m=match(unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))], colnames(demographics))
## remove NAs caused by the BDI and POMS not coming from the demographics frame
m=m[!is.na(m)]
## m now contains only the column numbers of the neuropsych variables we'll regress the %cs against

groups="mddOnly"

seeds=readSeedsFile(file.path(config.data.dir, "juelich_amygdala_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "juelich_bla_amygdala_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "Harvard-Oxford_amygdala_seeds.txt"))

numberOfSeeds=length(seeds)

totalNumberOfRegressions=numberOfSeeds * length(regressionVariables)

regressionCount=1
for ( regressionVariableCount in 1:length(regressionVariables ) ) {

    ## this is the version of the regressionVariable id that will be
    ## used in the model matrix and its accompaning formula and the
    ## output brik filename
    ##rvName=cleanRegressionVariable( regressionVariables[[regressionVariableCount]]$variable)
    rvVariable=regressionVariables[[regressionVariableCount]]$variable
    rvName=regressionVariables[[regressionVariableCount]]$name

    if (grepl("diff", rvVariable, fixed=TRUE)) {
        ## the temporal difference regressions have seperate directories 
        group.data.dir=normalizePath(file.path(data.dir, paste("Group.data", rvVariable, "withAandC", sep=".")))
        group.results.dir=normalizePath(file.path(data.dir, paste("Group.results", rvVariable, "withAandC", sep=".")))
    }
    
    for (seed in seeds) {
        seedName=getSeedName(seed)
        
        cat("####################################################################################################\n")
        cat(sprintf("*** Running RLM for the %s seed of the %s group (%02d of %02d regressions)\n", seedName, groups, regressionCount, totalNumberOfRegressions))
        
        ## setup all the filenames
        inputBucketFilename=paste("restingstate.bucket", groups, seedName,             "masked+tlrc.HEAD", sep=".")
        outputBucketPrefix= paste("restingstate",        groups, seedName, rvVariable, "rlm.bucket",       sep=".")
        
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
        
        ## now merge the subject order and the demographics file so we have
        ## the table completed for each subbrik in in the bucket file
        
        ## the demographics file will contain only one row for each subject
        ## but the subjectOrder file contains one row for each subject for
        ## each period type and hence will be much longer, this allows us to
        ## merge the two tables so we can complete the table for use in
        ## creation of the SPM. mgd stands for merged

        ## we branch on the regreeeion variable name because the MASC
        ## tscores used are not in the demographics df but are
        ## computed a few lines below from the MASC.total. By the time
        ## the model df is created the MASC.tscore variable will be in
        ## the mgd df as a result of calling computeMascScore
        if ( rvVariable == "MASC.tscore" ) {
            mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "MASC.total")])
        } else if (grepl("diff", rvVariable) ) {
            mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI")])
            if (rvVariable == "CDRSR.diff" ) {
                mgd=cbind(mgd, cdrsr.change[match(mgd$subject, cdrsr.change$SubjID), rvVariable])
                colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], rvVariable)

            } else if (rvVariable == "MASC.tscore.diff") {
                mgd=cbind(mgd, masc.change[match(mgd$subject, masc.change$SubjNum), rvVariable])
                colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], rvVariable)

            } else if (rvVariable == "BDI.diff") {
                mgd=cbind(mgd, bdi.change[match(mgd$subject, bdi.change$SubjNum), rvVariable])
                colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], rvVariable)

            } else if (rvVariable == "CGAS.diff") {
                mgd=cbind(mgd, cgas.change[match(mgd$subject, cgas.change$SubjNum), rvVariable])
                colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], rvVariable)

            } else if (rvVariable == "CDI.diff") {
                mgd=cbind(mgd, cdi.change[match(mgd$subject, cdi.change$SubjNum), rvVariable])
                colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], rvVariable)

            } else if (rvVariable == "RADS.Total.Tscore.diff") {
                mgd=cbind(mgd, rads.change[match(mgd$subject, rads.change$SubjNum), rvVariable])
                colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], rvVariable)
                
            } else {
                stop("Unknown rvVariable ", rvVariable, " when creating mgd data frame. Cannot continue.")
            }
        } else {
            mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", rvVariable)])
        }
        ## ensure that subject is s factor
        mgd$subject=as.factor(mgd$subject)
        mgd=droplevels(mgd)
        
        if ( rvVariable == "MASC.tscore" ) {
            mgd=fixDates(mgd)
            mgd=computeAge(mgd)
            
            mgd=computeMascScore(mgd)
        }

        print(head(mgd))

        mrData = inputBrik$brk
        dim(mrData) = c(dimX, dimY, dimZ, numberOfBriks)

        ## stop("Check mgd data frame\n")
        
        ## set up the data frame with the dependant variables we will want
        ## to include in the regression formula

        model = data.frame(
            as.vector(mgd[, "Grp"]),
            as.vector(mgd[, "subject"]),            
            as.vector(mgd[, rvVariable])
            )
        colnames(model) = c("Grp", "subject", rvVariable)

        printExcludedAndIncludedSubjects(mgd)
        ## print(head(model))
        
        ## stop("Check model data frame\n")
        
        ## set up the formula for the model that we want to regress
        modelFormula  = as.formula(paste("mri ~", rvVariable, sep=" "))
        
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

        ##stop("Check RLM summary\n")
        
        ## These are the labels that will be attached to the subbriks in the
        ## output stats bucket file. The labels will be dictated by the
        ## model formula and importantly the order in which the coefficients
        ## from the RLM coefficient matrix are concatenated
        ## the makeBrikLabels
        ## function does not include the mean of the fMRI at the start of
        ## the list of labels so include it here with the c()
        d=makeBrikLabels(tempRlmCoef, inBoot=doBootstrapping)
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
        if (doBootstrapping) {
            bootstrapStatsStartAt=d[["bootstrapLabelsStartAt"]] + 1
        }
        
        ##stop("Stopping")
        
        maxIter=100
        resamples=25
        Stats = array(0, c(dimX, dimY, dimZ, numberOfStatsBriks))
        cat(paste("Starting at", date(), "\n"))
        startTime=proc.time()
        if (useProgressBar) {
            pb <- txtProgressBar(min = 0, max = dimZ, style = 3)
        }
        if (ncpus > 1 ) {
            ## multiple cpus
            library(snow)
            
            cluster = makeCluster(ncpus, type = "SOCK")
            clusterEvalQ(cluster, library(MASS));

            ##function (inData, inNumberOfStatsBriks, inModel, inModelFormula, inMaxIt=50, inBoot=FALSE, inR=25, inBootstrapStatsStartAt=NA) {
            for ( kk in 1:dimZ) {
                if (useProgressBar) {
                    setTxtProgressBar(pb, kk)
                } else {
                    cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
                }
                Stats[ , , kk, ] = aperm(parApply(cluster, mrData[ , , kk, ],  c(1, 2), runRegression,
                         inNumberOfStatsBriks=numberOfStatsBriks, inModel=model, inModelFormula=modelFormula, inMaxIt=maxIter,
                         inBoot=doBootstrapping, inR=resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
            }
            stopCluster(cluster)
        } else {
            for ( kk in 1:dimZ ) {              
                ## single cpu
                if (useProgressBar) {
                    setTxtProgressBar(pb, kk)
                } else {
                    cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
                }
                Stats[ , , kk, ] = aperm(apply(mrData[ , , kk, ],  c(1, 2), runRegression,
                         inNumberOfStatsBriks=numberOfStatsBriks, inModel=model, inModelFormula=modelFormula, inMaxIt=maxIter,
                         inBoot=doBootstrapping, inR=resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
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
                       paste(paste("[", user, "@", hostname, ": ",  date(), "]", sep=""), file.path(script.location, script.name)),
                       paste("Model  formula:", gsub("~", "TILDE ", capture.output(print(modelFormula)))[1], sep=" ")),
                   origin = inputBrik$origin,
                   delta = inputBrik$delta,
                   orient= inputBrik$orient)

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
} ## end of for (term in regressionterms) {


