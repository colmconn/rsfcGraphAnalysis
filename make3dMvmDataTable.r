#!/usr/bin/Rscript

rm(list=ls())

library(getopt)

source("scoreMasc.r")

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

    if (is.null(in.opt$qvariables)) {
        cat("*** No covariate variable names were supplied.\n")
    } else {
        in.opt$qvariables=unlist(strsplit(in.opt$qvariables, "\\s", perl=TRUE, fixed=FALSE))
        ## print(in.opt$qvariables)
    }

    if (is.null(in.opt$bsvariables)) {
        cat("*** No between subjects variables were specified.\n")
    } else {
        in.opt$bsvariables==unlist(strsplit(in.opt$bsvariables, "\\s", perl=TRUE, fixed=FALSE))
        ## print(in.opt$bsvariables)
    }
    
    if (is.null(in.opt$wsvariables)) {
        cat("*** No within subject variable names were supplied.\n")
    } else {
        in.opt$wsvariables=unlist(strsplit(in.opt$wsvariables, "\\s", perl=TRUE, fixed=FALSE))
        ## print(in.opt$wsvariables)        
    }
    if (is.null(in.opt$wsvariables) && is.null(in.opt$bsvariables) && is.null(in.opt$wsvariables)) {
        cat("*** You must specify between, within and/or quantiative variables.\n")
    }
    
    if (is.null(in.opt$mask)) {
        cat("*** No mask was specified.\n")
        if ( ! interactive() )
            q(status=1)
    }

    if (is.null(in.opt$mask)) {
        cat("*** No mask was specified.\n")
    } else {
        mask.filename=file.path(group.results.dir, in.opt$mask)
        if ( ! file.exists( mask.filename)) {
            cat("*** No such file:", mask.filename, "\n")
            if ( ! interactive() )
                q(status=1)
        }
    }

    if (is.null(in.opt$seeds)) {
        cat("*** A file name containing the list of seeds to use is required.\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1)
    }

    if (is.null(in.opt$g2f)) {
        cat("*** A mapping of group names to file names is required.\n")
        cat(getopt(spec, usage=TRUE));    
        if ( !interactive() )
            q(status=1)
    } else {
        groupings=unlist(strsplit(in.opt$g2f, "\\s", perl=TRUE, fixed=FALSE))
        ## print(groupings)
        ## print(length(groupings))
        if (length(groupings) < 2) {
            cat("*** You must provide more than one (1) group to file name mapping\n")
            cat(getopt(spec, usage=TRUE));    
            if ( !interactive() )
                q(status=1) 
        } else {
            for (ii in 1:length(groupings) ) {
                ## cat("Groupings: ", groupings[ii], "\n")
                g2f=unlist(strsplit(groupings[ii], ":", perl=FALSE, fixed=TRUE))
                ## print(g2f)
                in.opt$groups.to.files.list[[g2f[1]]] = file.path(config.data.dir, g2f[2])
                ## print(groups.to.files.list)
                if ( ! file.exists(in.opt$groups.to.files.list[[ g2f[1] ]]) ) {
                    cat(sprintf("*** For the %s group, the subject list file %s does not exist\n", g2f[1], in.opt$groups.to.files.list[[ g2f[1] ]], "\n"))
                    if ( !interactive() )
                        q(status=1) 
                }
            }
        }
        ## cat("*** At end\n")
        ## print(in.opt$groups.to.files.list)
    }
    
    return(in.opt)
}

printDirectorySummary <- function () {

    cat("***\n")
    cat("*** Directories settings\n")
    cat("*** ====================\n")    
    for ( var in c(
        "root.dir",
        "scripts.dir",
        "data.dir",
        "admin.data.dir",
        "config.data.dir",
        "group.data.dir",
        "group.results.dir",
        "seeds.data.dir") ) {

        cat(sprintf("*** %s -> %s\n", var, eval(parse(text=var))))
    }
    cat("***\n")    
}


printOptionsSummary <- function (in.opt) {

    cat("*** Summary of ommand line arguments.\n")
    cat("*** =================================\n")
    cat("*** Prefix will be:", in.opt$prefix, "\n")
    cat("*** Within subject variable(s) will be:",  paste(in.opt$wsvariables, collapse=", "), "\n")
    cat("*** Between subject variable(s) will be:", paste(in.opt$bsvariables,  collapse=", "), "\n")
    cat("*** The following variable(s) will be used as covariates in the data table:\n")
    cat(paste("*** Quantative variable", 1:length(in.opt$qvariables), ":", in.opt$qvariables, collapse="\n"), sep="\n")
    if (in.opt$center) {
        cat("*** The quantative covariates listed above will be mean centered.\n")
    }

    cat("*** Mapping of group names to files containing list of subjects.\n")
    cat("*** -------------------------------------------------------------\n")
    cat("*** Group -> Subject list file\n")
    ## print(names(in.opt$groups.to.files.list))
    for (nn in names(in.opt$groups.to.files.list) ) {
        cat(sprintf("*** %s -> %s\n", nn,in.opt$groups.to.files.list[[ nn ]] ))
    }
    
    cat("*** Seeds will be read from:", file.path(config.data.dir, in.opt$seeds), "\n")
    cat("*** Mask will be:", file.path(group.results.dir, in.opt$mask), "\n")
    cat("***\n")    
}


readSubjectListFiles <- function( in.opt ) {

    subjectOrder=do.call( rbind,
        lapply(names(in.opt$groups.to.files.list), function (nn) { 
            so=read.table(in.opt$groups.to.files.list[[ nn ]])
            df=data.frame(rep(nn, 1:length(so)), so)
        }
               ))
    colnames(subjectOrder)=c("Group", "subject")

    return(subjectOrder)
}

splitSubjectId <- function ( in.subjectOrder )  {

    in.subjectOrder=data.frame(in.subjectOrder, do.call('rbind', strsplit(as.character( in.subjectOrder$subject ), '_', fixed=TRUE)))
    colnames(in.subjectOrder)=c("Group", "Subj", "subject", "timepoint")

    return(in.subjectOrder)
}

fixSubjectOrderTable <- function ( in.subjectOrderTable ) {

    in.subjectOrderTable$subject = as.character(in.subjectOrderTable$subject)
    in.subjectOrderTable[in.subjectOrderTable$subject=="300", "subject"] = "169/300"

    in.subjectOrderTable$subject=as.factor(in.subjectOrderTable$subject)

    return(in.subjectOrderTable)
}

makeInputFilename <- function (  in.subjectOrder, in.seed.name )  {

    df=cbind(in.subjectOrder,
        sapply(in.subjectOrder$Subj,
               function (xx) {
                   file.path(data.dir, xx, "rsfc", in.seed.name, sprintf("%s.z-score+tlrc.HEAD", in.seed.name))
               } )
                  )
    colnames(df)=c(colnames(in.subjectOrder), "InputFile")
    df$InputFile=as.character(df$InputFile)
    return (df)
}

checkInputFileExists <- function ( in.subjectOrder ) {

    files.exist=sapply(in.subjectOrder$InputFile,
        function (xx) {
            ## print(xx)
            if ( ! file.exists(xx) ) {
                cat("*** No such file:", xx, "\n")
                return(FALSE)
            } else {
                return(TRUE)
            }
        }
                       )
    if( ! any(files.exist) && ! interactive() ) {
        cat("*** Some subjects input files do not exist. Cannot continue\n")
    }
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

stopIfVariablesNotInModelMatrix <- function(in.opt, in.model.matrix ) {

    termDifference=setdiff(c(in.opt$bsvariables, in.opt$wsvariables, in.opt$qvariables), colnames(in.model.matrix))
    if (length(termDifference) > 0) {
        stop(paste("The following were in the list of between, within, and quantative variables but are not columns of the model matrix:", termDifference))
    }   
}

checkIsNa <- function (inData, inColumns) {
    for (column in inColumns) {
        
        if (any(is.na(inData[, column]))) {
            cat ("****************************************************************************************************\n")
            cat (sprintf("*** The following subjects have NA data for %s\n", column))

            print(data.frame ("Group" = as.vector ( inData[is.na(inData[, column]), "Group"]),
                              "Subj" = as.vector ( inData[is.na(inData[, column]), "Subj"])))
            
            cat ("****************************************************************************************************\n")      
        }
    } ## end of for (column in inColumns) {
    
} ## end of checkIsNa

computeMascScore <- function (inData) {
    inData.dim=dim(inData)
    for (r in seq(1, inData.dim[1]) ) {
        ## cat("##################################################\n")
        subjectNumber=inData[r, "subject"]
        gender=inData[r, "Gender"]
        age=round(inData[r, "age.in.years"], 0)
        old.masc.tscore=inData[r, "MASC.tscore"]

        ## cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f,\n", r, subjectNumber, gender, age, inData[r, "MASC.total"], old.masc.tscore))
        
        new.masc.tscore=scoreMasc(gender, age, inData[r, "MASC.total"])
        if (is.na(new.masc.tscore) ) {
            warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, inData[r, "MASC.total"]))
        }
        
        inData[r, "MASC.tscore"]=new.masc.tscore
        
        ## cat (sprintf("Old MASC tscore=%0.0f, new MASC tscore=%0.0f\n", old.masc.tscore, new.masc.tscore))
    }
    return (inData)
}

scaleCovariates <- function ( in.data) {
    for (col in opt$qvariables) {
        if ( is.numeric(in.data[, col]) ) {
            cat("*** Mean centering quantative variable:", col, "\n")            
            in.data[, col] = scale(in.data[, col], center=TRUE, scale=FALSE)
        } else {
            cat(sprintf("*** Skipping centering for %s: Not a numeric column\n", col))
        }
    }
    return(in.data)
}

makePrefix <- function (in.prefix, in.grouping, in.seed.name) {
    
    return(sprintf("%s.%s.%s.3dmvm.bucket.%s", opt$prefix, in.grouping, in.seed.name, format(Sys.time(), "%Y%m%d-%H%M%Z") ))
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

scripts.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts")
data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
group.data.dir=file.path(data.dir, "Group.data")
group.results.dir=file.path(data.dir, "Group.results")
seeds.data.dir=file.path(data.dir, "seeds")

wasi.column.names=c("Verbal", "Performance", "Full")

################################################################################
NO_ARGUMENT="0"
REQUIRED_ARGUMENT="1"
OPTIONAL_ARGUMENT="2"

## process command line arguments
spec = matrix(c(
    'help',          'h', NO_ARGUMENT,       "logical",
    "center",        'c', NO_ARGUMENT,       "logical",
    "seeds",         's', REQUIRED_ARGUMENT, "character",
    "g2f",           'g', REQUIRED_ARGUMENT, "character",

    "qvariables",    'q', REQUIRED_ARGUMENT, "character",
    "bsvariables",   'b', REQUIRED_ARGUMENT, "character",
    "wsvariables",   'w', REQUIRED_ARGUMENT, "character",

    "prefix",        'p', REQUIRED_ARGUMENT, "character",
    "mask",          'm', REQUIRED_ARGUMENT, "character",
    "model",         'o', REQUIRED_ARGUMENT, "character"
), byrow=TRUE, ncol=4)

################################################################################
if (interactive()) {
    ##
    ## THE FOLLOWING IS FOR TESTING PURPOSES ONLY
    ## 
    ## these arguments that are useful for testing purposes only.
    ##
    ## THE FOLLOWING IS FOR TESTING PURPOSES ONLY
    ##
    
    args=c(
        ## "-b", "-r", "5000", ## "-e",
        "-c",
        "-p", "restingstate",
        "-m", "mask.grey.Attempter.and.Non-attempter.union.masked+tlrc.HEAD",
        "-b", "Group",
        "-w", "Gender age.in.years",
        "-s", "juelich_amygdala_seeds_weights.txt",
        "-g", "Attempter:mdd.at.pilot.txt Non-attempter:mdd.nat.pilot.txt",
        ## "-q", "Full RSQ.fixed MASC.tscore"
        "-q", "Full"         
    )
    opt = getopt(spec, opt=args)
} else {
    opt = getopt(spec)
}
################################################################################
printDirectorySummary()
opt=checkCommandLineArguments(opt)
printOptionsSummary(opt)
    
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)
if (any(grepl("Group", colnames(demographics)))) {
    demographics$Group = NULL
}

wasiFilename=file.path(admin.data.dir, "WASI.csv")
wasi=readCsvFile(wasiFilename, inSubjectColumnName="SubID")

seeds.filename=file.path(config.data.dir, opt$seeds)

seeds=readSeedsFile(seeds.filename)

subjectOrder = fixSubjectOrderTable(splitSubjectId(readSubjectListFiles(opt)))
## print(subjectOrder)


mgd=cbind(subjectOrder,
    demographics[match(subjectOrder$subject, demographics$ID), ],
    ## demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "DOB", "MRI", "Gender")],    
    wasi        [match(subjectOrder$subject, wasi$SubID), 2:length(colnames(wasi))]
          )

if (any(mgd$subject=="378") ) {
    mgd[mgd$subject=="378", "Gender"]="F"
    mgd$subject=paste(mgd$subject, "A", sep="_")
}

mgd=fixDates(mgd)
mgd=computeAge(mgd)

mgd=computeMascScore(mgd)
stopIfVariablesNotInModelMatrix(opt, mgd)

mgd=mgd[, c("Subj", "subject", opt$bsvariables, opt$wsvariables, opt$qvariables)]
checkIsNa(mgd, opt$qvariables)

if (opt$center) {
    mgd=scaleCovariates(mgd)
}

rownames(mgd)=NULL
for (seed in seeds) {
    seedName=getSeedName(seed)

    cat("####################################################################################################\n")
    cat(sprintf("*** Creating covariate file for the %s seed of the %s grouping\n", seedName, paste(names(opt$groups.to.files.list), collapse=", ") ))

    mgd=makeInputFilename(mgd, seedName)
    ## print(mgd)

    checkInputFileExists(mgd)

    grouping=paste(names(opt$groups.to.files.list), collapse=".and.")

    dataTableFilename=file.path(group.data.dir, paste("dataTable", grouping, seedName, "txt", sep="."))
    cat("*** Writing data table to", dataTableFilename, "\n")
    write.table(mgd[1:dim(mgd)[1]-1, c("Subj", opt$bsvariables,  opt$wsvariables,  opt$qvariables, "InputFile")], file=dataTableFilename, quote=FALSE, col.names=TRUE, row.names=FALSE, eol=" \\\n")
    write.table(mgd[  dim(mgd)[1],   c("Subj", opt$bsvariables,  opt$wsvariables,  opt$qvariables, "InputFile")], file=dataTableFilename, quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)    

    three.d.mvm.command.script.filename=file.path(scripts.dir, sprintf("07-3dmvm.withCovariates.%s.%s.sh", grouping, seedName))
    cat("*** Writing the 3dMVM command to:", three.d.mvm.command.script.filename, "\n")
    cat("*** Note you MUST edit the command in this file and make sure that the correct formulae and contrasts are specified! It will not run properly as it stands!\n")

    three.d.mvm.command=sprintf("#!/bin/bash
%s
%s
3dMVM \\
      -prefix %s %s \\
      -jobs 8 \\
      -ranEff '~1' \\
      -SS_type 3 \\
      -qVars '%s' \\
      -bsVars '%s' \\
      -wsVars '%s' \\
      -num_glt  <n> \\
      -gltLabel 1 '< add GLT here>' \\
      -gltLabel 2 '< add GLT here>' \\
      -dataTable @%s\n",
        ifelse(opt$center, "\n### Variables in the data table are already mean centered\n", "\n"),
        paste("cd", group.results.dir),
        makePrefix(opt$prefix, grouping, seedName),
        ifelse(is.null(opt$mask), "", paste("-mask", opt$mask)),
        paste(opt$qvariables, collapse=","),
        paste(opt$bsvariables, collapse=","),
        paste(opt$wsvariables, collapse=","),               
        dataTableFilename)

    cat(three.d.mvm.command, file=three.d.mvm.command.script.filename)
}
