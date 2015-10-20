#!/usr/bin/Rscript

rm(list=ls())

library(MASS)
library(getopt)

checkCommandLineArguments <- function (in.opt) {
    ## if help was asked for print a friendly message
    ## and exit with a non-zero error code
    if ( !is.null(in.opt$help) ) {
        cat(getopt(spec, usage=TRUE));
        q(status=1);
    }

    if ( is.null(in.opt$debug) ) {
        in.opt$debug=FALSE
    }

    if (in.opt$debug)
        cat("*** Debug output enabled\n")

    if (is.null(in.opt$session)) {
        cat("A session directory is required\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1);
    }

    if (is.null(in.opt$prefix)) {
        cat("A prefix name is required.\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1);
    }

    if (is.null(in.opt$LHS)) {
        cat("A file name containing data for the left handside of the regression equation is required.\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1);
    }

    if (is.null(in.opt$RHS)) {
        cat("At least one file name containing data for the right handside of the regression equation is required.\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1);
    }
    else {
        in.opt$RHS=unlist(strsplit(in.opt$RHS, "[[:space:]]+", fixed=FALSE))
    }

    if (is.null(in.opt$robustregression) && is.null(in.opt$gramschmidt) ) {
        cat("*** No cleaning methods specified on the command line. Defaulting to robust regression\n")
        in.opt$robustregression=TRUE
        in.opt$gramschmidt=FALSE
    }

    if (is.null(in.opt$robustregression) && in.opt$gramschmidt ) {
        in.opt$robustregression=FALSE
    }

    if (in.opt$robustregression && is.null(in.opt$gramschmidt) ) {
        in.opt$gramschmidt=FALSE
    }
    
    if (in.opt$robustregression && in.opt$gramschmidt ) {
        stop("*** Only one of -r|--robustregression or -g|--gramschmidt should be provided. Not both. Exiting\n")
    }
    
    ## if (in.opt$robustregression ) {
    ##     in.opt$gramschmidt=FALSE
    ## }
    
    ## if (in.opt$gramschmidt) {
    ##     in.opt$robustregression=FALSE
    ## }
     
    return(in.opt)
}

printOptionsSummary <- function () {
    cat("*** Script name:", get_Rscript_filename(), "\n")
    cat("*** Session directory is set to:", opt$session, "\n")
    cat("*** Prefix is set to:", opt$prefix, "\n")
    cat("*** The LHS of the regression equation will be read from:", opt$LHS, "\n")
    cat("*** The RHS of the regression equation will be read from:", paste(opt$RHS, collapse=" "), "\n")
    if (opt$robustregression) {
        cat("*** Using robust regression remove signal from RHS in LHS\n")
    } else if (opt$gramschmidt) {
        cat("*** Using Gram-Schmidt orthoganization to render LHS and RHS orthogonal to one another\n")
    } else {
        stop("*** Don't know what method to use to clean the data. This should not happen. Stopping.")
    }
}

getSeedNameFromFilename <- function (in.name) {

    return (sub(".ts.1D", "", basename(in.name), fixed=TRUE))

}

makeTimeseriesFilenames <- function (in.timeseries.names) {
    fnames=sapply(in.timeseries.names,
        function(x) {
            tsFile=file.path(opt$session, x)
            return(tsFile)
        })
}

checkTimeSeriesFilesExist <- function (in.filenames) {
    files.exist=vector(mode = "logical", length = length(in.filenames))
    ## check that these files actually exist and quit
    ii=1
    for (file in in.filenames) {
        if ( ! file.exists (file) ) {
            cat(sprintf("*** No such file : %s\n", file))
            files.exist[ii]=FALSE
        } else {
            cat("*** Found", file, "\n")
            files.exist[ii]=TRUE            
        }
        ii=ii+1
    }

    return(all(files.exist))
}


loadTimeseriesFiles <- function (in.timeseries.files) {
    timeseries=do.call(cbind, lapply(in.timeseries.files,
        function(x) {
            rsfcTimeseries=read.table(x)
            return(rsfcTimeseries)
        }
    ))
    colnames(timeseries) = sapply(in.timeseries.files, getSeedNameFromFilename)
    return (timeseries)
}

stopIfFormulaTermsNotInModelMatrix <- function(inModelFormula, inModelMatrix) {

    termDifference=setdiff(all.vars(inModelFormula), colnames(inModelMatrix))
    if (length(termDifference) > 0) {
        stop(paste("The following terms were in the model forumlae but are not columns of the model matrix:", termDifference))
    }   
}

makeCleanedTimerseriesFilename <- function() {
    return(file.path(opt$session, paste(getSeedNameFromFilename(opt$LHS), "cleaned.ts.1D", sep=".")))
}

savedCleanedTimeseries <- function (in.model, in.cleaned.timeseries.filename) {
    cat("*** Cleaned timeseries for", opt$LHS, "will be saved to", in.cleaned.timeseries.filename, "\n")

    if (opt$robustregression) {
        write.table(residuals(in.model), in.cleaned.timeseries.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
    } else {
        write.table(in.model, in.cleaned.timeseries.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
    }
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

NO_ARGUMENT="0"
REQUIRED_ARGUMENT="1"
OPTIONAL_ARGUMENT="2"

## --session /data/sanDiego/rsfcGraphAnalysis/data/105_A/rsfc/L_CMA.weight.3mm.cleaned
## --prefix L_CMA.weight.3mm.cleaned
## --LHS L_CMA.weight.3mm.ts.1D
## --RHS "L_SFA.weight.3mm.ts.1D L_BLA.weight.3mm.ts.1D"


## process command line arguments
spec = matrix(c(
    'help',             'h', NO_ARGUMENT,       "logical",
    'debug',            'd', NO_ARGUMENT,       "logical",    
    'session',          's', REQUIRED_ARGUMENT, "character",
    'prefix',           'p', REQUIRED_ARGUMENT, "character",
    "LHS",              'L', REQUIRED_ARGUMENT, "character",
    "RHS",              'R', REQUIRED_ARGUMENT, "character",
    "robustregression", 'r', NO_ARGUMENT,       "logical",
    "gramschmidt",      'G', NO_ARGUMENT,       "logical"  
), byrow=TRUE, ncol=4)

if (interactive()) {
    ## these are default arguments that are useful for testing
    ## purposes.
    args=c(
        "-d", "-G",
        "--session", "/data/sanDiego/rsfcGraphAnalysis/data/105_A/rsfc/L_CMA.weight.3mm.cleaned",
        "--prefix", "L_CMA.weight.3mm.cleaned",
        "--LHS", "L_CMA.weight.3mm.ts.1D",
        "--RHS", "L_SFA.weight.3mm.ts.1D L_BLA.weight.3mm.ts.1D")
    opt = getopt(spec, opt=args)
} else {
    opt = getopt(spec)
}

opt=checkCommandLineArguments(opt)
printOptionsSummary()

timeseriesFilenames=makeTimeseriesFilenames(c(opt$LHS, opt$RHS))
cat("*** LHS and RHS Timeseries will be read from the following files:\n", paste("*** +++", timeseriesFilenames, collapse="\n"), "\n", sep="")
cat("*** Checking if the timeseries files exist:\n")
if (! checkTimeSeriesFilesExist(timeseriesFilenames)) {
    cat("*** Some of the timeseries files do not exist. See messages above\n")
}

timeseries=loadTimeseriesFiles(timeseriesFilenames)
if(opt$debug) {
    cat("*** Head of the timeseries data frame:\n")
    print(head(timeseries))
}

if (opt$robustregression) {
    regression.formula=paste(getSeedNameFromFilename(opt$LHS), "~", paste(sapply(opt$RHS, getSeedNameFromFilename), collapse=" + "))
    cat("*** Regression formula is set to:", regression.formula, "\n")
    regression.formula=as.formula(regression.formula)
    stopIfFormulaTermsNotInModelMatrix(regression.formula, timeseries)

    model=rlm(regression.formula, timeseries)
    if(opt$debug){
        cat("*** Regression model and its summary:\n")
        print(model)
        print(summary(model))
    }
} else if (opt$gramschmidt) {
    timeseries.qr=qr(timeseries)
    timeseries.ortho=qr.Q(timeseries.qr)
}

cleaned.timeseries.filename=makeCleanedTimerseriesFilename()
if (opt$robustregression) {
    savedCleanedTimeseries(model, cleaned.timeseries.filename)
    if (opt$debug) {
        cat("*** Head of the residuals:\n")
        print(head(residuals(model)))
    }
} else if (opt$gramschmidt) {
    savedCleanedTimeseries(timeseries.ortho[, 1], cleaned.timeseries.filename)
    if (opt$debug) {
        cat("*** Head of the orthogonal timeseries (only the left most column is saved to the file above:\n")
        print(head(timeseries.ortho))
    }
}
