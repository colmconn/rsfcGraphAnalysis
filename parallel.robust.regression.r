#!/usr/bin/Rscript

rm(list=ls())
library(nlme)
library(MASS)
library(getopt)

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
# in the bucket are correctly labeled. Spaces are replaced with . and
## ( or ) are deleted.
make.brik.labels <- function (inRlmCoef, inBoot=FALSE) {
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
## account that AFNI numbers from 0 not 1
make.afni.ttest.brik.ids <- function(inRlmCoef, inBoot=FALSE) {
    nrows=length(rownames(inRlmCoef))
    ncols=length(colnames(inRlmCoef))
    
    indices=seq(ncols, nrows*ncols, by=ncols)

    if (inBoot) {
        indices=c(indices, seq(indices[length(indices)] + grep("t.value", bootstrappingBrikLabelSuffxes), (nrows*ncols) + (nrows*length(bootstrappingBrikLabelSuffxes)), by=length(bootstrappingBrikLabelSuffxes)))
    }
    
    return (indices)
}

## given a list of AFNI brik indices and a degrees of freedom (inDf)
## make an appropriate list of statpar arguments to add to a 3drefit
## command line
make.afni.statpar.arguments <- function(inDf, inAfniTtestBrikIds) {
    return( paste(sapply(inAfniTtestBrikIds, function(x) { sprintf("-substatpar %d fitt %d", x[1], inDf) }), collapse=" ") )
}

make.statssym.list <- function (in.dof, in.afni.ttest.brik.ids) {
    lapply(in.afni.ttest.brik.ids, function (xx) list(sb=xx, typ="Ttest", par=in.dof))
}


## this is a very trimmed down version of the run.regression function
## below. It is only for use with bootstrapping
boot.regression<- function(inData, inIndices, inModelFormula, inMaxIt=50, inNumberOfBetaValues) {
    outStats <- vector(mode="numeric", length=inNumberOfBetaValues)
    if ( ! inherits(myrlm <- try(rlm(inModelFormula, data=inData[inIndices, ], maxit=inMaxIt), silent=TRUE),
                    "try-error") ) {
        ## cat ("length(coefficients(myrlm)) is ", coefficients(myrlm), "\n")
        return(coefficients(myrlm))
    } else {
        ## cat("Got an exception\n")
        return(outStats)
    }
}

run.regression <- function (inData, inNumberOfStatsBriks, inModel, inModelFormula, inMaxIt=50, inBoot=FALSE, inR=25, inBootstrapStatsStartAt=NA) {
    ## cat("inNumberOfStatsBriks =>", inNumberOfStatsBriks, "\n")
    ## cat("inBootstrapStatsStartAt =>", inBootstrapStatsStartAt, "\n")
    out.stats <- vector(mode="numeric", length=inNumberOfStatsBriks)
    ## full.model.residuals <- vector(mode="numeric", length=dim(inModel)[1])
    if ( ! all(inData == 0 ) ) {

        ## if inData is all zero then we are in a portion of the masked
        ## out data and should therefore not perform the rlm on it just
        ## return a vector of zeros, if it's not then we get in this
        ## branch and can try the rlm
        
        inModel$mri<-inData
        if ( ! inherits(myrlm <- try(rlm(inModelFormula, data = inModel, maxit=inMaxIt), silent=TRUE),
                        "try-error") ) {
            ## full.model.residuals=residuals(myrlm)
            numberOfBetaValues=length(coefficients(myrlm))
            ## print(myrlm)
            ## print(summary(myrlm))
            out.stats[1] = mean(inModel$mri)
            ## cat("out.stats[1] = ",  out.stats[1], "\n")
            
            ## print(summary(myrlm))
            ## > coef(summary(model))
            ## Value Std. Error    t value
            ## (Intercept) 89.09177404  8.0221830 11.1056770
            ## Age          0.04087267  0.1146234  0.3565821
            ## Educ         0.71157766  0.5370092  1.3250753
            ## as.vector(t(coef(summary(model))))
            ## [1] 89.09177404  8.02218304 11.10567705  0.04087267  0.11462344  0.35658210  0.71157766  0.53700922  1.32507531

            model.coefficients=as.vector(t(coef(summary(myrlm))))
            ## cat ("model.coefficients: ", model.coefficients, "\n")
            
            ## stop("Check stuff")
            
            if (inBoot) {
                ## cat("boot out.stats indices 2:(inBootstrapStatsStartAt-1)", 2:(inBootstrapStatsStartAt-1), "\n")
                out.stats[2:(inBootstrapStatsStartAt-1)]=model.coefficients
                
                if (is.na(inBootstrapStatsStartAt)) {
                    stop("*** ERROR in run.regression: inBootstrapStatsStartAt has not been set. It is currently NA. Cannot continue. Stopping\n")
                }

                ## cat ("inside inBoot out.stats is: ", out.stats, "\n")
                ## cat ("number of stats briks should be: ", length(2:(inBootstrapStatsStartAt-1)), "\n")
                
                boot.stats=boot(inModel, boot.regression, R=inR, inModelFormula=inModelFormula, inMaxIt=inMaxIt, inNumberOfBetaValues=numberOfBetaValues)
                ## cat("boot.stats is\n")
                ## print(boot.stats)
                ## print(class(boot.stats))
                ## print(is.vector(boot.stats))
                if (is(boot.stats, "boot")) {
                    ## these are vectors with as many columns as there are terms in
                    ## the regression model. Don't forget to include the intercept
                    ## (assuming there is one in the model) when you're trying
                    ## to mentalize this
                    ## really need to add the mean bootstrapped coefficient values to the output stats vector
                    boot.beta.mean=apply(boot.stats$t, 2, mean)
                    boot.bias=boot.beta.mean - boot.stats$t0
                    boot.se.coeff=apply(boot.stats$t, 2, sd)
                    boot.t.coeff=boot.stats$t0 / boot.se.coeff

                    ## cat("boot.beta.mean is: ", boot.beta.mean, "\n")
                    ## cat("boot.bias is     : ", boot.bias, "\n")                
                    ## cat("boot.se.coeff is : ", boot.se.coeff, "\n")
                    ## cat("boot.t.coeff is  : ", boot.t.coeff, "\n")

                    start.at=inBootstrapStatsStartAt
                    for (term.index in seq(1, length(boot.t.coeff))) {
                        ## the normal element of the CI value contains 3 elements: 1) the CI
                        ## level (in this case 0.95), 2) the lower bound on the CI, 3) the
                        ## upper bound on the CI
                        ci=boot.ci(boot.stats, conf = c(0.95), type = c("norm"), index = term.index)$normal[2:3]
                        ## cat("ci for", names(boot.stats$t0)[term.index], "is:", ci, "\n")
                        ## cat("out.stats indices start.at:(start.at+3)", start.at:(start.at+3), "\n")
                        out.stats[start.at:(start.at+4)] = c(boot.beta.mean[term.index], boot.bias[term.index], boot.t.coeff[term.index], ci)
                        start.at=start.at+5
                    }
                }
            } else {
                ## cat("no boot out.stats indices 2:inNumberOfStatsBriks", 2:inNumberOfStatsBriks, "\n")
                out.stats[2:inNumberOfStatsBriks]=model.coefficients
            }
        }
    } ## end of if ( ! all(inData == 0 ) ) {

    ## cat ("out.stats is: ", out.stats, "\n")
    ## longest.label=max(sapply(output.stats.brik.labels, nchar))
    ## cc=cbind (output.stats.brik.labels, out.stats)
    ## cc=cbind(seq(1:dim(cc)[1]), cc)
    ## cat(apply(cc, 1,
    ##           function(xx) {
    ##               sprintf("%02d: %s => %0.5f",
    ##                       as.integer(xx[1]),
    ##                       format(xx[2], justify="left", width=longest.label),
    ##                       as.numeric(xx[3]), 5)
    ##           } ), sep="\n")
    
    ## if ( ! all(out.stats == 0 ) )
    ##     stop()
    ## cat(".")

    ##return(list("stats"=out.stats, "full.model.residuals"=full.model.residuals))
    return(out.stats)
}


read.data.table <- function (inFilename) {
    if (file.exists(inFilename)) {
        cat("*** Reading", inFilename, "\n")
        data.table=read.table(inFilename, header=TRUE, allowEscapes=TRUE)
    } else {
        stop(paste("*** No such file:", inFilename, "\n"))
    }
    data.table$InputFile=as.character( data.table$InputFile)
    return(data.table)
}

conditionally.make.session.dir <- function (in.opt) {
    if ( ! dir.exists(in.opt$session)) {
        if (in.opt$verbose)
            cat("*** Recursively creating", in.opt$session, "\n")
        dir.create(in.opt$session, recursive=TRUE)
    }
}


help <- function(){

}

check.command.line.arguments <- function (in.opt) {
    ## if help was asked for print a friendly message
    ## and exit with a non-zero error code
    if ( !is.null(in.opt$help) ) {
        cat(getopt(spec, usage=TRUE));
        q(status=1);
    }

    if (is.null(in.opt$datatable)) {
        cat("A data table is required\n")
        cat(getopt(spec, usage=TRUE));
        q(status=1);
    }

    if ( is.null(in.opt$bootstrap)) {
        in.opt$bootstrap=FALSE
    }
    
    if (in.opt$bootstrap) {
        ## load the boot strapping library only if bootstrapping is
        ## requested
        library(boot)
    }

    if (is.null(in.opt$formula)) {
        cat("A formula is required\n")
        cat(getopt(spec, usage=TRUE));
        q(status=1);
    }

    if (is.null(in.opt$session)) {
        cat("A session directory name (where to save results) is required.\n")
        cat(getopt(spec, usage=TRUE));    
        q(status=1);
    }

    if ( is.null(in.opt$infix)) {
        cat("An infix used in the naming of the errts and stats outpiut files is required.\n")
        cat(getopt(spec, usage=TRUE));
        q(status=1)
    }

    if (is.null(in.opt$resamples)) {
        in.opt$resamples=25
    }

    if ( is.null(in.opt$progress)) {
        in.opt$progress=FALSE
    }

    if ( is.null(in.opt$threads)) {
        in.opt$threads=1
    } else {
        if (in.opt$threads < 1) {
            cat(sprintf("*** Value for the number of threads is too low (%d). Setting to 1\n", in.opt$threads))
            in.opt$threads=1
        } else if (in.opt$threads > max.cpus) {
            cat(sprintf("*** Value for the number of threads is too high (%d). Setting to maximum number of available CPUs (%d)\n", in.opt$threads, max.cpus))
            in.opt$threads=max.cpus
        }
    }

    if ( is.null(in.opt$verbose)) {
        in.opt$verbose=FALSE
    } else {
        in.opt$verbose=TRUE
    }
    
    
    return(in.opt)
}

print.command.line.arguments.summary <- function () {
    cat("*** Data table will be read from", opt$datatable, "\n")
    if (opt$bootstrap) {
        cat("*** Performing bootstrapping of the regression using", opt$resamples, "resamples\n")
    } else {
        cat("*** Bootstrapping will NOT be performed\n")
    }
    cat("*** Model formula is:", opt$formula, "\n")
    cat("*** Results will be written to:", opt$session, "\n")
    cat("*** Stats and errts files will be named with", opt$infix, "infix\n")
    if (is.null(opt$mask)) {
        cat("*** No maksing will be performed\n")
    } else {
        cat("*** Mask will be read from", opt$mask, "\n")
    }
    cat(paste("*** Running with", opt$threads, ifelse(opt$threads == 1, "CPU", "CPUs"), "\n"))
    if (opt$verbose) {
        cat("*** Verbose messages enabled\n")
    }
    if (opt$progress) {
        cat("*** Using progress bar\n")
    }
}


stop.if.not.valid.data.table <- function(in.data.table) {
    if (names(in.data.table)[1] != "Subj" ) {
        stop(paste("*** The first column of the data table must be named 'Subj' not", names(in.data.table)[1], "\n"))
    }
    if (names(in.data.table)[length(names(in.data.table))] != "InputFile" ) {
        stop(paste("*** The last column of the data table must be named 'InputFile' not", names(in.data.table)[length(names(in.data.table))], "\n"))
    }
    do.input.files.exist=sapply(in.data.table$InputFile, file.exists)
    if(! any(do.input.files.exist)) {
        cat(paste("*** The following input files do not exist.\n",
                  paste(in.data.table$InputFile[!do.input.files.exist], collapse="\n"), "\n"))
        stop("*** Cannot continue\n")
    }
    unique.brik.file.views=unique(vapply(data.table[, "InputFile"], view.AFNI.name, ""))
    if (length(unique.brik.file.views) > 1) {
        cat("*** More than one view provided for the input files in the data table\n")
        cat(paste("*** The following views were derived from the input files in the data table:",
                  paste(unique.brik.file.views, collapse=", "), "\n"))
        stop("*** Cannot continue\n")
    }

    n.nonmandatory.columns=length(colnames(in.data.table)[-c(1, length(colnames(in.data.table)))])
    if (n.nonmandatory.columns < 1) {
        stop("*** Data table must contain at least 1 (one) non mandatory column to conduct regression\n")
    }
    if (opt$verbose) {
        cat("*** Data table passes validity checks\n")
    }
}

print.data.table.summary <- function(in.data.table) {
    cat(sprintf("*** Read data for %s subjects\n",  length(in.data.table[, "Subj"])))
    nonmandatory.columns=colnames(in.data.table)[-c(1, length(colnames(in.data.table)))]
    for (cc in nonmandatory.columns) {
        cat(sprintf("*** Data table column %s, of class %s, is summarized as follows:\n", cc, class(in.data.table[, cc])))
        print(summary(in.data.table[, cc]))

        print(addmargins(table(in.data.table[, c("Subj", cc)])))
    }
}

get.model.terms <- function(in.formula) {
    model.variables=all.vars(as.formula(in.formula))
    model.variables=model.variables[grep("mri", model.variables, fixed=TRUE, invert=TRUE)]
    return(model.variables)
}


stop.if.not.valid.formula <- function (in.formula, in.data.table) {
    if (! grepl("mri", in.formula, fixed=TRUE)) {
        stop(paste("*** Model formula does not contain the term 'mri' as either a explanatory or observed variable.",
                   "Cannot continue", sep="\n"))
    }

    model.variables=get.model.terms(in.formula)
    
    termDifference=setdiff(model.variables, colnames(in.data.table)[-c(1, length(names(in.data.table)))])

    if (length(termDifference) > 0) {
        stop(paste("The following terms were in the model forumla but are not columns of the data table:", paste(termDifference, collapse=", ")))
    }

    if (opt$verbose) {
        cat("*** Formula passes validity checks\n")
    }
}


read.input.briks <- function (in.data.table, in.opt) {
    cat("*** Reading each subjects' input data file. This may take some time\n")
    ## culled from 3dLME.R
    ## inputBrik=unlist(lapply(lapply(in.data.table[, "InputFile"], read.AFNI, verb=in.opt$verbose, meth="clib", forcedset = TRUE), '[[', 1))
    input.briks=lapply(in.data.table[, "InputFile"], read.AFNI, verb=in.opt$verbose, meth="clib", forcedset = TRUE)
    cat("*** Done\n")
    all.input.brik.dims=do.call(rbind, lapply(input.briks, function(ib) ib$dim))

    ## the folowing line of code is useful for testing whether the
    ## dimension checking code below works or not
    ## all.input.brik.dims[1, 1]=all.input.brik.dims[1, 1]+1
    ## now check that the dimensions of all of the briks read in are the same. If they are not we have a problem
    if ( ! all(apply(all.input.brik.dims, 2, function (xx) { abs(max(xx) - min(xx)) < .Machine$double.eps }))) {
        dims.table=cbind(in.data.table[ , c("Subj", "InputFile")], all.input.brik.dims)
        colnames(dims.table)=c("Subj", "InputFile", "dim.x", "dim.y", "dim.z", "dim.t")
        cat("*** The dimensions of all of all of the input briks is not the same\n")
        print(dims.table)
        stop("*** Cannot continue\n")
    } else {
        input.brik.dims=unique(all.input.brik.dims)
        input.brik.dims[4]=length(in.data.table$InputFile)
        mr.data=unlist(lapply(input.briks, '[[', 1))
        dim(mr.data)=input.brik.dims
    }

    return(mr.data)
}

read.first.inputfile <- function(in.data.table) {
    cat("*** Reading first subject's input data file\n")
    ## culled from 3dLME.R
    ## inputBrik=unlist(lapply(lapply(in.data.table[, "InputFile"], read.AFNI, verb=in.opt$verbose, meth="clib", forcedset = TRUE), '[[', 1))
    first.brik=read.AFNI(in.data.table[1, "InputFile"], verb=opt$verbose, meth="clib", forcedset = TRUE)

    return(first.brik)
}

read.mask <- function() {
    cat("*** Reading mask file\n")
    ## culled from 3dLME.R
    ## inputBrik=unlist(lapply(lapply(in.data.table[, "InputFile"], read.AFNI, verb=in.opt$verbose, meth="clib", forcedset = TRUE), '[[', 1))
    mask=read.AFNI(opt$mask, verb=opt$verbose, meth="clib", forcedset = TRUE)

    return(mask$brk)
}

binarize.mask <- function (in.mask, in.mr.data) {
    if (opt$verbose) {
            cat("*** Binarizing mask\n")
    }
    bb=as.numeric(abs(in.mask) > 0)
    dim(bb)=dim(in.mr.data)[1:3]
    return (bb)
}

check.mask.and.input.brik.dimensions.match <- function(in.mask, in.mr.data) {

    input.and.mask.dims=rbind(dim(in.mask)[1:3], dim(in.mr.data)[1:3])
    ## print(input.and.mask.dims)
    ## the following line is useful to check the code below whether
    ## checks the dimensions of the mask and mr.data match
    ## input.and.mask.dims[1, 1] = input.and.mask.dims[1, 1] + 2
    ## print(input.and.mask.dims)    
    rownames(input.and.mask.dims)=c("Input Briks'", "Mask brik")
    colnames(input.and.mask.dims)=c("x", "y", "z")
    if (! all(apply(input.and.mask.dims, 2, function (xx) { abs(max(xx) - min(xx)) < .Machine$double.eps }))) {
        cat("*** The dimensions of the input briks and the mask brik are not the same\n")
        print(input.and.mask.dims)
        stop("*** Cannot continue\n")
    } else {
        if (opt$verbose) {
            cat("*** Dimensions of input briks and mask brik match\n")
        }
    }
}

test.model <-function(in.model.formula, in.model, in.mr.data, in.opt) {
    ## if we got this far we should be able to run a single voxel to work
    ## out the DOF for the various stats, generate labels and indices for
    ## use with the unlisted rlm in the runRlm function
    i = 21
    j = 21
    k = 31

    if (opt$verbose)
        cat(sprintf("*** Testing model at voxel (i, j, k)=(%d, %d, %d)\n", i, j, k))
    in.model$mri = mr.data[i, j, k, ]
    ## mdl=rlm(in.model.formula, data = in.model)
    ## if (opt$verbose)
    ##     print(mdl)

    ## print(in.model)
    ## lme.model=lme(fixed=in.model.formula, random=as.formula("~ 1 | Subj"),
    ##               data=in.model)
    ## print(lme.model)
    ## stop()
    
    if( inherits(test.rlm <- try(rlm(in.model.formula, data = in.model),
                                silent=FALSE),
                 "try-error") ) {
        test.rlm.coeff=0
        test.rlm.dof=0
        traceback()
        stop("*** Got an exception trying to setup the test.rlm.coeff and test.rlm.dof variables. Cannot continue beyond this point. Stopping.")
    } else {
        s=summary(test.rlm)
        if (opt$verbose) {
            cat("*** The test robust linear model summary is:\n")
            print(s)
        }
        test.rlm.coeff = coef(s)
        test.rlm.dof=s$df
        ## rlm.residuals=residuals(test.rlm)
        ## if (in.opt$verbose) {
        ##     cat("*** The test robust linear model residuals are:\n")
        ##     print(rlm.residuals)
        ## }
        if (opt$verbose)
            cat("*** Success! Model appears to be solvable\n")
    }

    return(list("test.rlm.coeff"=test.rlm.coeff, "test.rlm.dof"=test.rlm.dof))##, "test.rlm.residuals"=rlm.residuals))
           
}


##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

## enable debugging of runRlm
##debug(run.regression)
##trace("run.regression", quote(if(! all(inData == 0 ) ) browser()))
####################################################################################################

if ( Sys.info()["sysname"] == "Darwin" ) {
    max.cpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
} else if ( Sys.info()["sysname"] == "Linux" ) {
    max.cpus=as.integer(system("cat /proc/cpuinfo |grep -c processor", intern=T))
} else {
    max.cpus=1
}

bootstrappingBrikLabelSuffxes=c("booted.mean.beta", "bias", "booted.t.value", "ciLower", "ciUpper")

################################################################################
NO_ARGUMENT="0"
REQUIRED_ARGUMENT="1"
OPTIONAL_ARGUMENT="2"

## process command line arguments
command.line.argument.specification = matrix(c(
    'bootstrap',	'b',	NO_ARGUMENT,		'logical',
    'datatable',	'd',	REQUIRED_ARGUMENT,	'character',
    'formula',		'f',	REQUIRED_ARGUMENT,	'character',
    'help',		'h',	NO_ARGUMENT,		'logical',
    'infix',		'i',	REQUIRED_ARGUMENT,	'character',
    'progress',		'p',	NO_ARGUMENT,		'logical',
    'resamples',	'r',	REQUIRED_ARGUMENT,	'integer',
    'session',		's',	REQUIRED_ARGUMENT,	'character',
    'threads',		't',	REQUIRED_ARGUMENT,	'integer',
    'mask',             'm',    REQUIRED_ARGUMENT,	'character',
    'verbose',		'v',	NO_ARGUMENT,		'logical'
), byrow=TRUE, ncol=4)

if (interactive()) {
    cat("*** Setting interactive options\n")

    args=c(
        "-v", "-r", "100", "-t", "2", "-t", "1", "-p",
        "-f", "CDRS.t.score.scaled.diff ~ mri + age.in.years",
        "-d", "/data/sanDiego/machLearnT1Analysis/data/config/data.table.mdd.CDRS.t.score.L_whole_amygdala.3mm.tab",
        "-s", "/data/sanDiego/machLearnT1Analysis/data/group.results",
        "-i", "restingstate.mddOnly.L_whole_amygdala3.mm.CDRS.t.score.scaled.diff"
    )
    ##the name of this script
    script.name=parent.frame(2)$ofile
    ## the location (absolute path) to this script
    script.location=normalizePath(dirname(parent.frame(2)$ofile))
    
    command.line=paste(file.path(script.location, script.name), paste(args, collapse=" "))
    opt = getopt(command.line.argument.specification, opt=args)
} else {
    args.start.at=grep("--args", commandArgs(), fixed=TRUE)
    command.line=paste(get_Rscript_filename(), paste(commandArgs()[-c(1:args.start.at)], collapse=" "))
    opt = getopt(command.line.argument.specification)
}

opt=check.command.line.arguments(opt)
print.command.line.arguments.summary()

## stop("Check the output of the getopt processing\n")

data.table=read.data.table(opt$datatable)
## data.table$InputFile=file.path("/Volumes", data.table$InputFile)
stop.if.not.valid.data.table(data.table)
print.data.table.summary(data.table)
stop.if.not.valid.formula(opt$formula, data.table)
conditionally.make.session.dir(opt)

## setup all the filenames
stats.bucket.prefix=paste("stats",    opt$infix,        sep=".")
errts.bucket.prefix=paste("errts",    opt$infix,        sep=".")
cluster.log.filename=paste("cluster", opt$infix, "log", sep=".")


mr.data = read.input.briks(data.table, opt)

if (! is.null(opt$mask)) {
    mask=read.mask()
    ## print(dim(mask))
    ## print(dim(mr.data))
    check.mask.and.input.brik.dimensions.match(mask, mr.data)
    mask=binarize.mask(mask, mr.data)
    ## apply the mask to the input data
    if (opt$verbose) 
        cat("*** Applying mask to input data\n")
    mr.data = array(apply(mr.data, 4, function(x) x*mask), dim=dim(mr.data))
    if (opt$verbose) 
        cat("*** Done\n")
    ## print(mr.data[ , , 21, 1])
    ## stop()
}

## stop("Check the input briks\n")

## set up the data frame with the dependant variables we will want
## to include in the regression formula

cat("*** Setting up model data frame\n")
model=data.table[, c("Subj", get.model.terms(opt$formula))]

cat("*** The head of the model data frame is as follows:\n")
print(head(model))

## stop("Check model data frame\n")

model.formula=as.formula(opt$formula)

rlm.values=test.model(model.formula, model, mr.data, opt)
test.rlm.coeff=rlm.values[["test.rlm.coeff"]]
test.rlm.dof=rlm.values[["test.rlm.dof"]]
## number.of.residuals=length(rlm.values[["test.rlm.residuals"]])

## stop("Check the test model\n")

## These are the labels that will be attached to the subbriks in the
## output stats bucket file. The labels will be dictated by the
## model formula and importantly the order in which the coefficients
## from the RLM coefficient matrix are concatenated
## the make.brik.labels
## function does not include the mean of the fMRI at the start of
## the list of labels so include it here with the c()
sub.brik.labels=make.brik.labels(test.rlm.coeff, inBoot=opt$bootstrap)

output.stats.brik.labels=c("Mean", sub.brik.labels[["labels"]])

## we add 1 to both number.of.stats.briks and number.of.stats.briks
## becasue make.brik.labels does not take into account the
## fact that we will add an additional subbrik (the mean) to
## the output stats. Hence the output of make.brik.labels is
## always 1 too small

## the number of stats subbriks to write out. This is dictated by the
## model Formula, changes to it likely imply changes to this number
number.of.stats.briks=length(output.stats.brik.labels)

bootstrapStatsStartAt=NA
if (opt$bootstrap) {
    bootstrapStatsStartAt=sub.brik.labels[["bootstrapLabelsStartAt"]] + 1
}

## stop("Stopping")

maxIter=25
mr.data.dims=dim(mr.data)

statistics.array = array(0, c(mr.data.dims[1], mr.data.dims[2], mr.data.dims[3], number.of.stats.briks))
## residuals.array  = array(0, c(mr.data.dims[1], mr.data.dims[2], mr.data.dims[3], number.of.residuals))
cat(paste("*** Starting at", date(), "\n"))
start=Sys.time()

if (opt$threads > 1 ) {
    ## multiple cpus
    library(snow)
    seed.value=123456789
    ## cluster = makeCluster(ncpus, type = "SOCK")
    cat("*** Starting cluster... ")
    cluster = makeCluster(rep("localhost", opt$threads), type = "SOCK", outfile=file.path(opt$session, cluster.log.filename))
    clusterEvalQ(cluster, library(MASS))
    clusterEvalQ(cluster, library(boot))
    clusterExport(cluster, c("boot.regression", "run.regression"))
    clusterSetupRNG(cluster, seed=seed.value)
    cat("Done\n")
    ##function (inData, inNumberOfStatsBriks, inModel, inModelFormula, inMaxIt=50, inBoot=FALSE, inR=25, inBootstrapStatsStartAt=NA) {
    if (opt$progress) {
        pb <- txtProgressBar(min = 0, max = mr.data.dims[3], style = 3)
    }
    for ( kk in 1:mr.data.dims[3]) {
        if (opt$progress) {
            setTxtProgressBar(pb, kk)
        } else {
            cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
        }
        statistics.array[ , , kk, ] =
            aperm(parApply(cluster, mr.data[ , , kk, ],  c(1, 2), run.regression,
                           inNumberOfStatsBriks=number.of.stats.briks, inModel=model, inModelFormula=model.formula, inMaxIt=maxIter,
                           inBoot=opt$bootstrap, inR=opt$resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
    }
    stopCluster(cluster)
} else {
    if (opt$progress) {
        pb <- txtProgressBar(min = 0, max = mr.data.dims[3], style = 3)
    }
    for ( kk in 1:mr.data.dims[3] ) {
        ## for ( kk in 31:31 ) {            
        ## single cpu
        if (opt$progress) {
            setTxtProgressBar(pb, kk)
        } else {
            cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
        }
        statistics.array[ , , kk, ] =
            aperm(apply(mr.data[ , , kk, ],  c(1, 2), run.regression,
                        inNumberOfStatsBriks=number.of.stats.briks, inModel=model, inModelFormula=model.formula, inMaxIt=maxIter,
                        inBoot=opt$bootstrap, inR=opt$resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
    }
    ## do something with the resudials here
    
    ## the lines below is useful for debugging the run.regression function, particularly when using boot strapping
    ## i = 21
    ## j = 21
    ## k = 31
    ## out.stats.list = run.regression(mr.data[i, j, k, ], inNumberOfStatsBriks=number.of.stats.briks, inModel=model, inModelFormula=model.formula, inMaxIt=maxIter,
    ##     inBoot=opt$bootstrap, inR=opt$resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt)
}
cat("\n")
if (opt$progress) {
    close(pb)
}

end=Sys.time()
cat("*** Finished running regression models at:", date(), "\n")
suppressMessages(library(chron))
cat("*** Computation took", format(as.chron(end) - as.chron(start)), "\n")

## now if save the HEAD/BRIK and run 3drefit to set the statistical parameters, view, space and ID, 
if (is.matrix(test.rlm.coeff)) {

    ## this is used to get the header for saving the stats brik file
    first.brik=read.first.inputfile(data.table)

    rlm.out.brik.filename=paste(stats.bucket.prefix, ".", format(Sys.time(), "%Y%m%d-%H%M%Z"), sep="")
    rlm.out.brik.fqfn=file.path(opt$session, rlm.out.brik.filename)
    cat("*** Writing bucket file ", rlm.out.brik.fqfn, "\n")

    command.history=gsub("~", "TILDE", command.line)
    afni.ttest.brik.ids=make.afni.ttest.brik.ids(test.rlm.coeff, inBoot=opt$boot)
    statssym.list=make.statssym.list(test.rlm.dof[2], afni.ttest.brik.ids)
    write.AFNI(rlm.out.brik.fqfn,
               statistics.array,
               statsym=statssym.list,
               verb=opt$verbose,
               view=view.AFNI.name(data.table[1, "InputFile"]),
               addFDR=FALSE,
               label=output.stats.brik.labels,
               com_hist=command.history,
               defhead=first.brik)

    cat ("*** Making statpar arguments\n")
    statpar = "3drefit"
    
    statpar = paste("(cd",  opt$session, ";", statpar, "-addFDR",
        paste(rlm.out.brik.filename, view.AFNI.name(data.table[1, "InputFile"]), sep=""),
        ")")    
    cat("*** Running", statpar, "\n")
    system(statpar)
} else {
    stop("*** The test.rlm.coeff was not a matrix. Could not save the results HEAD/BRIK\n")
}
   
