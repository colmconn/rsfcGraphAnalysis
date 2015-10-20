rm(list=ls())
graphics.off()

##library(gmodels)
library(gdata)
library(ggplot2)
##library(lubridate)
library(compute.es)
library(reshape)
library(orddom)
library(car)

source("scoreMasc.r")

########################################################################################################################################################################################################
### START OF FUNCTIONS #################################################################################################################################################################################
########################################################################################################################################################################################################

stack <- function(){ 
    it <- list() 
    res <- list( 
        push=function(x){ 
            it[[length(it)+1]] <<- x 
        }, 
        pop=function(){ 
            val <- it[[length(it)]] 
            it <<- it[-length(it)] 
            return(val) 
        }, 
        value=function(){ 
            return(it) 
        } 
        ) 
    class(res) <- "stack" 
    res 
} 
print.stack <- function(x,...){ 
    print(x$value()) 
} 
push <- function(stack,obj){ 
    stack$push(obj) 
} 
pop <- function(stack){ 
    stack$pop() 
}

make.significance.indications <- function(pValues, which.pValues=c(1)) {

    Signif=symnum(pValues, corr=FALSE, na=FALSE, cutpoints = c(0,  .001,.01, .05, .1, 1),
        symbols   =  c("***", "**", "*", ".", " "))
    f=format(Signif)

    ## only return the first one as we're only interested in marking significant group effects
    return(f[which.pValues])
}

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


## Reads the seed file and does the (crude) equivalent of BAS variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character())
    table=gsub("$scriptsDir", scriptsDir, table, fixed=TRUE)

    return (table)
}

## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed
getSeedName <- function(inSeedPath){
    name=basename(inSeedPath)
    return(gsub("\\.nii.*", "", name))
}


makeTableString <- function(inGroup, inMean, inSe, inMin, inMax, inNaCount, inMissingData=TRUE) {
    ##  st=paste(round(inMean, 1), " / ", round(inSe, 1),    

    st=paste(round(inMean, 1), " ± ", round(inSe, 1),
        " (", round(inMin, 1), "-", round(inMax, 1), ")", ifelse(inMissingData & inNaCount > 0, paste(" [", inNaCount, "]", sep=""), ""), sep="")
    return(st)
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


analyse <- function(inData, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=NULL, results.stack=stack()) {

    if (length(levels(group.variable)) > 2 ) {
        stop("analyse: Cannot handle more then two levels in the grouping variable. Quitting.\n")
    }

    inData=inData[complete.cases(inData[, group.variable]), ]
    inData[, group.variable] = as.factor( inData[, group.variable])
    
    group1=levels(inData[, group.variable])[1]
    group2=levels(inData[, group.variable])[2]    
    
    cat("****************************************************************************************************\n")
    cat(sprintf("*** Using %s as the grouping variable. Group1: %s Group2: %s\n", group.variable, group1, group2))
    

    inData$Grp=drop.levels(inData[ , group.variable])
    inData$Gender=drop.levels(inData$Gender)
    ## inData$Race=drop.levels(inData$Race)
    inData$Race=as.factor(tolower(drop.levels(inData$Race)))
    
    cat("\n*** Gender\n")
    gender.table=table(inData[, c("Gender", group.variable)])
    gender.test=prop.test(gender.table)
    gender.table=addmargins(gender.table)
    print(gender.table)
    print(gender.test)

    ##ethnicity.table=table(inData$ethnicity, inData$Group)
    ## print(gender.table)

    cat("\n*** Number of participants (n)\n")
    n.table=table(inData[, c("Grp")])
    n.test=prop.test(n.table)
    n.table=addmargins(n.table)
    print(n.table)
    print(n.test)
    
    ## csvLine=paste("Number of participants in final analysis (n)", n.table[group1], n.table[group2],
    ##     sprintf("Chi(%0.2f) = %0.2f",
    ##             round(n.test$parameter, 2), round(n.test$statistic, 2)),
    ##     round(n.test$p.value, 2), "", make.significance.indications(n.test$p.value), sep=",")

    csvLine=make.proportions.test.string.for.results.stack(n.table, n.test, group1, group2, "Number of participants in final analysis (n)")
    push(results.stack, csvLine)

    ## we dont use the make.proportions.test.string.for.results.stack because it cant handle the gender split
    csvLine=paste("Gender (M / F)", paste(gender.table["M", group1], gender.table["F", group1], sep=" / "),
        paste(gender.table["M", group2], gender.table["F", group2], sep=" / "),
        sprintf("Chi(%0.2f) = %0.2f",
                round(gender.test$parameter, 2), round(gender.test$statistic, 2)),
        round(gender.test$p.value, 2), "", make.significance.indications(gender.test$p.value), sep=",")
    push(results.stack, csvLine)
    
    ## cat("\n*** Race\n")
    ## race.table=table(inData[, c("Race", "Grp")])
    ## ##ethnicity.table=table(inData$ethnicity, inData$Group)
    ## print(addmargins(race.table))
    ## print(prop.test(race.table))
    
    cat("*** Now performing tests on psychMeasures\n")
    
    for (i in 1: length(psychMeasures)) {
        variable=psychMeasures[[i]]$variable
        name=psychMeasures[[i]]$name
        
        ##cat("################################################################################\n");
        ##cat("Summary for ", name, "variable = ", variable, "\n")
        if ( ! variable %in% colnames(inData)) {
            message (sprintf("Cannot find variable (%s) in the data frame provided. Skipping\n", variable))            
        } else if (is.factor(inData[, variable])) {
            message (sprintf("%s, is a factor. Skipping\n", variable))
        } else {
            
            ## if (any(inData[, variable] < 0, na.rm=TRUE )) {
            ##     cat ("****************************************************************************************************\n")
            ##     cat (sprintf("*** The following subjects have data less than zero for %s (%s)\n", name, variable))
            ##     print(as.vector(inData[inData[, variable] < 0 & ! is.na(inData[, variable]), "subject"]))
            ##     print(as.vector(inData[inData[, variable] < 0 & ! is.na(inData[, variable]), "Grp"]))
            ##     cat ("*** Now setting these values to NA\n")
            ##     inData[inData[, variable] < 0 & !is.na(inData[, variable]), variable]=NA
            ##     inData[inData[, variable] < 0 & !is.na(inData[, variable]), variable]=NA      
            ##     cat ("****************************************************************************************************\n")      
            ## }
            
            sm.df=summarySE(inData, measure=variable, groupvars=c(group.variable), na.rm=TRUE)
            ##print(sm.df)
            ##print(rownames(sm.df))
            
            group1.string=""
            group2.string=""
            test=0
            if (isNonparametricTestVariable(variable)) {
                group1.string=makeTableString(sm.df[1, 1], inMean=sm.df[1, "median"],  sm.df[1, "IQR"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
                group2.string=makeTableString(sm.df[2, 1], inMean=sm.df[2, "median"],  sm.df[2, "IQR"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)
                
                ##name=paste(name, "†", sep="")
            } else {
                group1.string=makeTableString(sm.df[1, 1], sm.df[1, variable],  sm.df[1, "sd"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
                group2.string=makeTableString(sm.df[2, 1], sm.df[2, variable],  sm.df[2, "sd"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)
            }

            if (any(is.na(inData[, variable]))) {
                cat (sprintf("*** The following subjects have missing data for %s (%s)\n", name, variable))
                cat(as.vector(inData[is.na(inData[, variable]), "subject"]), "\n")
                cat(as.vector(inData[is.na(inData[, variable]), group.variable]), "\n")      
            }

            ## cat("*** Analyzing ", variable, "\n")
            if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {
                if (isNonparametricTestVariable(variable)) {
                    ##cat("Control length: ", length(inData[inData[, group.variable]==group1, variable]), "MDD length: ", length(inData[inData[, group.variable]==group2, variable]), "\n")
                    if( inherits(test <- try(wilcox.test(inData[inData[, group.variable]==group1, variable],
                                                         inData[inData[, group.variable]==group2, variable]),
                                             silent=FALSE),
                                 "try-error") ) {
                        test <- 0
                    }
                } else {
                    
                    if( inherits(test <- try(t.test(inData[inData[, group.variable]==group1, variable],
                                                    inData[inData[, group.variable]==group2, variable]),
                                             silent=FALSE),
                                 "try-error") ) {
                        test <- 0
                    }
                }
                ## print(test)
            } else {
                cat ("*** Insufficient opservations\n")
            } ## end of if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {

            var.statistic=""
            var.df=""
            var.pvalue=""
            var.parameter=""
            var.significance=""
            var.effect.size=""
            
            if (is.list(test)) {
                var.statistic=round(test$statistic, 2)
                var.df=ifelse(is.numeric(test$parameter), round(test$parameter, 2), "")
                var.pvalue=round(test$p.value, 2)
                var.parameter=ifelse(is.null(test$parameter), "NA", round(test$parameter, 2))
                var.significance=make.significance.indications(test$p.value)

                if (isNonparametricTestVariable(variable)) {
                    ## compute Probability of Superiority here
                    if (var.pvalue < 1.0) {

                        if ( compute.effect.size ) {
                            orddom.ps   =dmes     (na.omit(inData[inData[, group.variable]==group1, variable]), na.omit(inData[inData[, group.variable]==group2, variable]))
                            orddom.ps.ci=dmes.boot(na.omit(inData[inData[, group.variable]==group1, variable]), na.omit(inData[inData[, group.variable]==group2, variable]), theta.es="PSc")

                            var.effect.size=round(orddom.ps$PSc, 3)
                            ## upper bound on the 95% Confidence interval for the effect size
                            es.ci.lb=round(orddom.ps.ci$theta.bci.lo, 2)
                            ## upper bound on the 95% Confidence interval for the effect size                        
                            es.ci.ub=round(orddom.ps.ci$theta.bci.up, 2)
                        }
                        else {
                            ## This is the method for calculating the
                            ## Probability of Superioroty Score
                            ## mentioned in 1. Erceg-Hurn DM,
                            ## Mirosevich VM (2008) Modern robust
                            ## statistical methods: an easy way to
                            ## maximize the accuracy and power of your
                            ## research. Am Psychol 63:
                            ## 591–601. doi:10.1037/0003-066X.63.7.591.
                            ##
                            ## var.effect.size=round(
                            ##     var.statistic /
                            ##     ( length(inData[inData[, group.variable]==group1, variable]) *
                            ##       length(inData[inData[, group.variable]==group2, variable]) ), 2)
                            
                            var.effect.size=0
                            es.ci.ub=0
                            es.ci.lb=0
                        }
                        
                        st = sprintf("%s,%s,%s,W = %0.0f,%0.2f,PS=%s (%s; %s),%s",
                            name, paste(group1.string, "†", sep=""), paste(group2.string, "†", sep=""),
                            round(var.statistic, 2),
                            round(var.pvalue, 2),  var.effect.size, es.ci.lb , es.ci.ub, var.significance)
                    } else {
                        st = sprintf("%s,%s,%s,W = %0.0f,%0.2f,,%s",
                            name, paste(group1.string, "†", sep=""), paste(group2.string, "†", sep=""),
                            round(var.statistic, 2),
                            round(var.pvalue, 2), var.significance )
                    }
                } else {
                    if (var.pvalue < 1.0) {
                        if ( compute.effect.size ) {
                            es=tes(var.statistic, length(inData[inData[, group.variable]==group1, variable]), length(inData[inData[, group.variable]==group2, variable]), verbose=FALSE)
                            ## upper bound on the 95% Confidence interval for the effect size
                            es.ci.lb=round(es$l.g, 2)
                            ## upper bound on the 95% Confidence interval for the effect size                        
                            es.ci.ub=round(es$u.g, 2)
                            var.effect.size=round(es$g, 2)
                        } else {
                            var.effect.size=0
                            es.ci.ub=0
                            es.ci.lb=0
                        }
                            
                        st = sprintf("%s,%s,%s,t(%0.2f) = %0.2f,%0.2f,g=%s (%0.2f; %0.2f),%s",
                            name, group1.string, group2.string,
                            round(var.parameter, 2), round(var.statistic, 2),
                            round(var.pvalue, 2),    var.effect.size, es.ci.lb, es.ci.ub, var.significance )
                    } else {
                        st = sprintf("%s,%s,%s,t(%0.2f) = %0.2f,%0.2f,,%s",
                            name, group1.string, group2.string,
                            round(var.parameter, 2), round(var.statistic, 2),
                            round(var.pvalue, 2),    var.significance )
                    }
                }
                
            } else {
                st = sprintf("%s,%s,%s", name, group1.string, group2.string)
            } ## end of if (is.list(test)) {

            ## run levern's test and append the results

            homogeneityFormula = as.formula(paste(variable, group.variable, sep=" ~ "))
            ## print(homogeneityFormula)
            homogeneityTest=leveneTest(homogeneityFormula, data=inData)
            ## print(homogeneityTest)
            ht.string=sprintf("F(%d; %d)=%0.3f; p=%0.3f", homogeneityTest[["Df"]][1], homogeneityTest[["Df"]][2], homogeneityTest[["F value"]][1], homogeneityTest[["Pr(>F)"]][1])
            st=paste(st, ht.string, sep=",")
            ## stop("Check the homogeneity test results\n")
            
            ## st=paste(name, ctrl.string, mdd.string,
            ##     var.parameter, var.statistic, var.pvalue, var.effect.size, var.significance, sep=",")
            push(results.stack, st)
        } ## end of if (is.factor(inData[, variable])) { ... } else { ... }
    }
    if ( ! is.null(results.table.filename)) {
        cat("*** Results table is in ", results.table.filename, "\n")
        ff=file(results.table.filename, open="w", encoding="utf-8")
        sink(ff, append=FALSE)
    } else {
        cat("################################################################################\n");
        cat("Summary statistics table\n")
    }
    l=results.stack$value()
    header=sprintf("Characteristic,%s,%s,Stat.,pValue,Effect Size (95%% CI),Signif.,Levene's Test", group1, group2)
    cat(header, "\n")
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
    if ( ! is.null(results.table.filename)) {
        sink()
    }
}

isNonparametricTestVariable <- function (inVariableName) {

    ## the C_ will capture the Caprara Irritablity ans emotional
    ## susceptibility scales
    ##nonParametricRegexp="Tanner|SES|PSWQ|CGI.CGAS|CoRum|RSQ|Hand|C_|CTQ"
    ## nonParametricRegexp="Tanner|SES|PSWQ|CGI.CGAS|CoRum|RSQ|Hand|C_|num.stressful.events|num.severe.events|total.sum.stress"
    nonParametricRegexp="Tanner|SES|PSWQ|CoRum|RSQ|Hand|C_|num.stressful.events|num.severe.events|total.sum.stress"        

    if ( inVariableName == "NS" ) {
        return (TRUE)
    } else if (inVariableName == "HA" ) {
        return (TRUE)
    } else if (inVariableName == "RD" ) {
        return (TRUE)
    } else if (inVariableName == "P" ) {
        return (TRUE)
    } else if (inVariableName == "SD" ) {
        return (TRUE)
    } else if (inVariableName == "C" ) {
        return (TRUE)
    } else if (inVariableName == "ST" ) {
        return (TRUE)
    } else if (any(grep(nonParametricRegexp, inVariableName)) ) {
        return (TRUE)
    } else {
        return (FALSE)
    }
}


filterOutStimulusColumn <- function (inDataFrame) {
    df=inDataFrame
    if (dim(inDataFrame)[2] > 1) {
        cat("*** Filtering the data frame because it has more then 1 column\n")
        df=data.frame("subject"=unique(inDataFrame[, 1]))
    }
    return (df)
}

checkLessThanZero <- function (inData, inSetToNa=FALSE) {
    variablesLessThanZero=FALSE
    for (i in 1: length(psychMeasures)) {
        variable=psychMeasures[[i]]$variable
        name=psychMeasures[[i]]$name

        cat(sprintf("*** Checking whether %s (%s) is less than zero\n", name, variable))
        
        if (any(inData[, variable] < 0, na.rm=TRUE )) {
            variablesLessThanZero=TRUE
            cat ("****************************************************************************************************\n")
            cat (sprintf("*** The following subjects have data less than zero for %s (%s)\n", name, variable))
            cat(paste(as.vector(inData[inData[, variable] < 0 & ! is.na(inData[, variable]), "subject"]), collapse=" "), "\n")
            cat(paste(as.vector(inData[inData[, variable] < 0 & ! is.na(inData[, variable]), "Grp"]), collapse=" "), "\n")
            cat(paste(as.vector(inData[inData[, variable] < 0 & ! is.na(inData[, variable]), variable]), collapse=" "), "\n")
            if (inSetToNa) {
                cat ("*** Now setting these values to NA\n")
                inData[inData[, variable] < 0 & !is.na(inData[, variable]), variable]=NA
            }
            cat ("****************************************************************************************************\n")      
        }
    } ## end of for (i in 1: length(psychMeasures)) {
    if(variablesLessThanZero & ! inSetToNa) {
        stop("Found variables with values less than zero. Cannot continue\n")
    }
    return(inData)
} ## end of checkLessThanZero

checkIsNa <- function (inData, inColumns, inDrop=FALSE) {
    for (column in inColumns) {
        
        if (any(is.na(inData[, column]))) {
            cat ("****************************************************************************************************\n")
            cat (sprintf("*** The following subjects have NA data for %s\n", column))

            cat (paste (as.vector ( inData[is.na(inData[, column]), "subject"]), collapse=" "), "\n")
            cat (paste (as.vector ( inData[is.na(inData[, column]), "Grp"]), collapse=" "), "\n")            
            
            ##cat(paste(as.vector(is.na(inData[inData[, column]) & ! is.na(inData[, column]), "subject"]), collapse=" "), "\n")
            ##cat(paste(as.vector(is.na(inData[inData[, column]) & ! is.na(inData[, column]), "Grp"]), collapse=" "), "\n")
            ##cat(paste(as.vector(is.na(inData[inData[, column]) & ! is.na(inData[, column]), column]), collapse=" "), "\n")
            if (inDrop) {
                cat ("*** Nowdropping those subjects\n")
                inData = inData[ complete.cases(inData[, column]), ]
            }
            cat ("****************************************************************************************************\n")      
        }
    } ## end of for (column in inColumns) {
    
    return(inData)
} ## end of checkIsNa

checkColumnAgainstThreshold <- function (inData, inColumns, inThreshold, inDirection=c("greater", "less")) {
    inDirection=match.arg(inDirection)

    retVal=FALSE
    
    for (column in inColumns) {
        ## print(inData[, c("subject", "Grp", column)])
    
        if (inDirection == "greater") {
            overThresholdRows=inData[, column] > inThreshold
            if (any(overThresholdRows)) {
                ## print(inData[, column] > inThreshold)
                cat ("****************************************************************************************************\n")
                cat (sprintf("*** The following subjects have %s data greater than %f\n", column, inThreshold))
                
                cat (paste (as.vector ( inData[which(overThresholdRows), "subject"]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(overThresholdRows), "Grp"    ]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(overThresholdRows), column   ]), collapse=" "), "\n")                
                cat ("****************************************************************************************************\n")
                retVal=TRUE
            }
        } else if (inDirection == "less") {
            underThresholdRows=inData[, column] < inThreshold
            if (any(underThresholdRows)) {
                cat ("****************************************************************************************************\n")
                cat (sprintf("*** The following subjects have %s data less than %f\n", column, inThreshold))
                
                cat (paste (as.vector ( inData[which(underThresholdRows), "subject"]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(underThresholdRows), "Grp"    ]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(underThresholdRows), column   ]), collapse=" "), "\n")                
                cat ("****************************************************************************************************\n")                
            }
            retVal=TRUE            
        }
    }
    return(retVal)
}

checkForMedicatedSubjects <- function(inData) {
    retVal=FALSE
    ## Now remove subjects that are NOT medication-naive (Prozac, Zoloft, Celexa, Citalopram, Klondpin, Seroquil, Cymbalta and ADHD meds)
    medicatedSubjects= inData$subject %in% c("130", "132", "318", "319", "320", "322", "323", "324", "325", "329", "345", "376", "384", "333", "349")
    if (any(medicatedSubjects)) {
        cat ("****************************************************************************************************\n")
        cat ("*** The following medicated subjects were found in the analysis\n")
        cat (paste (as.vector (  inData[which(medicatedSubjects), "subject"]), collapse=" "), "\n")
        cat (paste (as.vector (  inData[which(medicatedSubjects), "Grp"]), collapse=" "), "\n")
        cat ("****************************************************************************************************\n")
        retVal=TRUE
    }

    return(retVal)
}

checkForNonMaleOrFemale <- function(inData) {
    retVal=FALSE
    tgSubjects=! inData$Gender %in% c("M", "F") 
    if (any(tgSubjects) ) {

        cat ("****************************************************************************************************\n")
        cat ("*** The following transgender subjects were found in the analysis\n")
        cat (paste (as.vector (  inData[which(tgSubjects), "subject"]), collapse=" "), "\n")
        cat (paste (as.vector (  inData[which(tgSubjects), "Grp"]), collapse=" "), "\n")
        cat (paste (as.vector (  inData[which(tgSubjects), "Gender"]), collapse=" "), "\n")        
        cat ("****************************************************************************************************\n")
        

        retVal=TRUE
    }
    
    return(retVal)
}


graphVariables <- function (inData) {
    my.base.size=14
    my.theme=
        theme_bw(base_size =  my.base.size) +
            theme(
                ##legend.position="none",
                legend.position="bottom",        
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),

                ##remove the panel border
                ##panel.border = element_blank(),

                ## add back the axis lines
                axis.line=element_line(colour = "grey50"),
                legend.text = element_text(angle=45),
                
                ##axis.title.x=element_blank(),
                axis.title.x = element_text(size=my.base.size, vjust=0),
                axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
                plot.title=element_text(size=my.base.size*1.2, vjust=1))

    imageDirectory=file.path(group.results.dir, "psychMeasureGraphs")
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }
    
    for (i in 9: length(psychMeasures)) {
        variable=psychMeasures[[i]]$variable
        name=psychMeasures[[i]]$name
        
        ss=inData[, c("Grp", variable)]
        graph = ggplot(ss, aes_string(x=variable)) +
            geom_histogram(aes(y=..density.., fill=..density..)) +
                geom_density(aes(color="red")) +
                    labs(x=name) +
                        my.theme
        imageFilename=file.path(imageDirectory, sprintf("%s.pdf", variable))
        cat(paste("*** Creating", imageFilename, "\n"))
        ggsave(imageFilename, graph)
         ##print(graph)
        ##stop("Check histogram\n")
    } ## end of for (i in 1: length(psychMeasures)) {
} ## end of graphVariables


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

make.proportions.test.string.for.results.stack <- function(in.table, in.prop.test, in.group1, in.group2, in.characteristic.string) {

    csvLine=paste(in.characteristic.string, in.table[in.group1], in.table[in.group2],
        sprintf("Chi(%0.2f) = %0.2f",
                round(in.prop.test$parameter, 2), round(in.prop.test$statistic, 2)),
        round(in.prop.test$p.value, 2), "", make.significance.indications(in.prop.test$p.value), sep=",")

    return(csvLine)
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

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
scriptsDir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))

##demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
ctqFilename=file.path(admin.data.dir, "Exisiting CTQ Scored_with RC.csv")
slesFilename=file.path(admin.data.dir, "SLES_20140820.csv")
wasiFilename=file.path(admin.data.dir, "WASI.csv")


subjectList=c("106_A", "107_A", "108_A", "109_A", "111_A", "112_A", "114_A",
"116_A", "121_A", "122_A", "123_A", "124_A", "126_A", "130_A",
"131_A", "132_A", "133_A", "134_A", "136_A", "137_A", "138_A",
"141_A", "142_A", "143_A", "144_A", "145_A", "146_A", "147_A",
"149_A", "150_A", "151_A", "152_A", "155_A", "157_A", "158_A",
"160_A", "161_A", "162_A", "164_A", "165_A", "167_A", "301_A",
"303_A", "304_A", "306_A", "310_A", "312_A", "313_A", "316_A",
"317_A", "318_A", "320_A", "321_A", "322_A", "323_A", "324_A",
"325_A", "326_A", "328_A", "331_A", "332_A", "333_A", "335_A",
"336_A", "339_A", "342_A", "343_A", "344_A", "346_A", "347_A",
"348_A", "349_A", "350_A", "351_A", "353_A", "354_A", "355_A",
"356_A", "357_A", "358_A", "359_A", "360_A", "361_A", "363_A",
"364_A", "365_A", "366_A", "367_A", "368_A", "370_A", "371_A",
"372_A", "373_A", "377_A", "379_A", "380_A", "386_A", "387_A")

## Create the subjectOrder variable used throughout the rest of the script
## the clean subject list has no timepoint suffix
subjectOrder=data.frame("subject"=c(subjectList))
colnames(subjectOrder)=c("subject")
subjectOrder$subject=as.vector(subjectOrder$subject)
## due to some sort of screw up 169/300_A is the same subject so fix it here.
badSubjectIds=which (subjectOrder$subject=="300_A")
if (length(badSubjectIds) > 0 ) {
    subjectOrder[badSubjectIds, "subject"]="169/300_A"
}
subjectOrder$subject=as.factor(gsub("_A[0-9]?", "", as.character(subjectOrder$subject), fixed=FALSE))
cat(sprintf("*** The subjectOrder data frame has %s unique subjects\n",  length(unique(subjectOrder$subject))))

###stop("Stooooooooooooooooooooooooooooooooooooooooooooooooooping")

cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ""))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

cat("*** Reading", ctqFilename, "\n")
ctq=read.csv(ctqFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
cat(sprintf("*** Read CTQ data for %s unique subjects\n",  length(unique(ctq$ID))))

cat("*** Reading", slesFilename, "\n")
sles=read.csv(slesFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
cat(sprintf("*** Read SLES data for %s unique subjects\n",  length(unique(sles$subject))))
sles=sles[, c(1:92)]

cat("*** Reading", wasiFilename, "\n")
wasi=read.csv(wasiFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
cat(sprintf("*** Read WASI-II data for %s unique subjects\n",  length(unique(wasi$SubID))))

selectedColumns=c(
    "Grp", "Gender", "DOB", "MRI", "Race", "Hand", "SES", "Tanner1", "Tanner2", "TannerAvg",
    "CGI.CGAS",                 "CDRS.raw",                  "CDRS.tscore",               "SFMQ.total",              
    "MADRS.raw",                "PSWQ", 		     "CoRum",
    "RADS.DM",                  "RADS.AN",                   "RADS.NS",                   "RADS.SC",
    "RADS.DM.Tscore",           "RADS.AN.Tscore",            "RADS.NS.Tscore",            "RADS.SC.Tscore",  	"RADS.Total.Tscore",
    "BDI.II",                   "CDI",                       "MASC.total",                "MASC.tscore",        "RSQ.fixed",
    "SCARED.som.panic",         "SCARED.gen.anx",            "SCARED.seper.anx",          "C_Irritability",	"C_EmoSuscept",
    "SCARED.soc.phobia",        "SCARED.school.phobia",      "BIS",                       "BAS.drive",                
    "BAS.funseek",              "BAS.reward.responsiveness", "PsycAche.total", "NS",	"HA",	"RD",	"P",	"SD",	"C",	"ST")

if (length(setdiff(selectedColumns, colnames(demographics))) != 0) {
    message ("*** The following column name(s) is/are causing problems: ", setdiff(selectedColumns, colnames(demographics)), "\n")
    stop("There is a mismatch between the selected columns and the columns in the demographics data frame. Stopping.")
    
}

####################################################################################################
### Set up the data frames for statistical analysis
####################################################################################################

mgd=cbind(
    subjectOrder,
    demographics    [match(subjectOrder$subject, demographics$ID), selectedColumns],
    wasi            [match(subjectOrder$subject, wasi$SubID), c("Verbal", "Performance", "Full")],
    ctq             [match(subjectOrder$subject, ctq$ID), c("CTQ_TOTAL", "CTQ_EA", "CTQ_PA", "CTQ_SA", "CTQ_EN", "CTQ_PN", "CTQ_MD")])


## add this column with the subject ID concatenated with the timepoint
## so we can pull the correct data from the SLES data frame. This is
## done because the SLES data frame contains data for all time points.
mgd$subjectWithTimePoint=as.factor(paste(mgd$subject, "A", sep="_"))
mgd=cbind(
    mgd,
    sles [match(mgd$subjectWithTimePoint, sles$subject),
          c("num.stressful.events", "num.severe.events", "total.sum.stress", "num.1.little.to.no.effect", "num.2.some.effect", "num.3.moderate.effect", "num.4.great.effect")])
## remove the now unnecessary subject with timepoint column
mgd$subjectWithTimePoint=NULL

mgd$Grp=drop.levels(mgd$Grp)
mgd$SCARED.som.panic=as.numeric(as.vector(mgd$SCARED.som.panic))
mgd$MASC.tscore=as.numeric(mgd$MASC.tscore)

mgd=fixDates(mgd)
mgd=computeAge(mgd)

mgd=computeMascScore(mgd)
rownames(mgd)=NULL

## stop("stopppppppppping")

psychMeasures=list(
    list(variable="age.in.years", name="Age at time of scan (years)"),
    ## list(variable="Hand",         name="Edinburgh Handedness Inventory"),
    list(variable="SES",          name="Hollingshead Socioeconomic Score"),
    list(variable="TannerAvg",    name="Tanner Score"),
    ## list(variable="Tanner2",    name="Tanner Score"),    

    ## list(variable="WASI.PERF",    name="Wechsler Abbreviated Scale of Intelligence (Performance)"),
    ## list(variable="WASI.Full.4",  name="Wechsler Abbreviated Scale of Intelligence"),
    
    ## use these for the seperate WASI table    
    list(variable="Verbal",    name="Wechsler Abbreviated Scale of Intelligence (Verbal)"),
    list(variable="Performance",    name="Wechsler Abbreviated Scale of Intelligence (Performance)"),    
    list(variable="Full",           name="Wechsler Abbreviated Scale of Intelligence (Full)"),
    
    
    list(variable="CGI.CGAS",     name="Children's Global Assessment Scale"),
    ## list(variable="PSWQ",         name="Penn State Worry Questionnaire"),
    ## list(variable="CoRum",        name="Corumination Questionnaire"),
    ## list(variable="RSQ.fixed",    name="Ruminative Responses Styles Questionnaire"),
    
    list(variable="CDRS.tscore",  name="Children's Depression Rating Scale (Standardized)"),
    
    ## list(variable="RADS.DM",      name="Reynolds Adolescent Depression Scale Dysphoric Mood"),
    ## list(variable="RADS.AN",      name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect"),
    ## list(variable="RADS.NS",      name="Reynolds Adolescent Depression Scale Negative Self-evaluation"),
    ## list(variable="RADS.SC",      name="Reynolds Adolescent Depression Scale Somatic Complaints"),    
    
    list(variable="RADS.DM.Tscore",      name="Reynolds Adolescent Depression Scale Dysphoric Mood (Standardized)"),
    list(variable="RADS.AN.Tscore",      name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect (Standardized)"),
    list(variable="RADS.NS.Tscore",      name="Reynolds Adolescent Depression Scale Negative Self-evaluation (Standardized)"),
    list(variable="RADS.SC.Tscore",      name="Reynolds Adolescent Depression Scale Somatic Complaints (Standardized)"),
    list(variable="RADS.Total.Tscore",   name="Reynolds Adolescent Depression Scale Total (Standardized)"),
    
    
    ## list(variable="SFMQ.total",     name="SFMQ"),
    ## list(variable="MADRS.raw",      name="Montgomery-Asberg Depression Rating Scale"),
    list(variable="BDI.II",         name="Beck Depression Inventory II"),
    list(variable="CDI",            name="Children's Depression Inventory"),
    ##list(variable="MASC.total",     name="Multidimensional Anxiety Scale for Children"),
    list(variable="MASC.tscore",    name="Multidimensional Anxiety Scale for Children (Standardized)")
    ## list(variable="C_Irritability", name="Caprara Irritability Scale"),
    ## list(variable="C_EmoSuscept",   name="Caprara Emotional Susceptibility Scale")
    
    ## list(variable="SCARED.som.panic", name="SCARED Som. Panic"),
    ## list(variable="SCARED.gen.anx", name="SCARED Ganeralized Anxiety"),
    ## list(variable="SCARED.seper.anx", name="SCARED Seperation Anxiety"),
    ## list(variable="SCARED.soc.phobia", name="SCARED Social Phobia"),

    ## list(variable="BIS", name="Behavioral Inhibition System"),
    ## list(variable="BAS.drive", name="Behavioral Approach System: Drive"),
    ## list(variable="BAS.funseek", name="Behavioral Approach System: Fun Seeking"),
    ## list(variable="BAS.reward.responsiveness", name="Behavioral Approach System: Reward Responsiveness"),


    ## ## TCI related variables
    ## list(variable="NS", name="TCI: Novelty Seeking"),
    ## list(variable="HA", name="TCI: Harm Avoidance"),
    ## list(variable="RD", name="TCI: Reward Dependence"),
    ## list(variable="P",  name="TCI: Persistence"),
    ## list(variable="SD", name="TCI: Self-Directedness"),
    ## list(variable="C",  name="TCI: Cooperativeness"),
    ## list(variable="ST", name="TCI: Self-Transcendence"),
    
    ## list(variable="PsycAche.total", name="PsycAche.total")

    ## list(variable="CTQ_TOTAL", name="Childhood Trauma Questionnaire (Total)"),

    ## list(variable="CTQ_EA", name="CTQ: Emotional Abuse"),
    ## list(variable="CTQ_PA", name="CTQ: Physical Abuse"),
    ## list(variable="CTQ_SA", name="CTQ: Sexual Abuse"),
    ## list(variable="CTQ_EN", name="CTQ: Emotional Neglect"),
    ## list(variable="CTQ_PN", name="CTQ: Physical Neglect"),
    ## list(variable="CTQ_MD", name="CTQ: Minimization/Denial"),

    ## list(variable="num.stressful.events", name="SLES: Number of stressful events"),
    ## list(variable="num.severe.events",    name="SLES: Number of severe events"),
    ## list(variable="total.sum.stress",     name="SLES: Total Sum Stress")
    )

cat("####################################################################################################\n")
cat("### Demographic and Psychiatric Statistics\n")

##sink("../data/Group.results/restingStateStatictics28Nov.txt")
## analyse(mgd, "../data/Group.results/restingStateDistributionOfGenders.png", "../data/Group.results/restingStateDistributionOfRaces.png")

##stop("check out the orddom package\n")

mgd=checkLessThanZero(mgd, inSetToNa=TRUE)
depressedControls=FALSE
underdepressed=FALSE
medicatedSubjects=FALSE

if(checkColumnAgainstThreshold(mgd[mgd$Grp=="NCL", ], "CDRS.tscore", 54, inDirection="g")) {
    cat("You have depressed control subjects included in the analysis\n")
    depressedControls=TRUE
}

if(checkColumnAgainstThreshold(mgd[mgd$Grp=="MDD", ], "CDRS.tscore", 55, inDirection="l")) {
    cat("You have underdepressed MDD subjects included in the analysis\n")
    underdepressed=TRUE
}

if(checkForMedicatedSubjects(mgd)) {
    cat("You have medicated subjects included in the analysis\n")
    medicatedSubjects=TRUE
}

## stopifnot(!depressedControls, !underdepressed, !medicatedSubjects)


## if (checkForNonMaleOrFemale(mgd)) {
##     warning("Cannot continue with subjects whose gender is neither male or female included in the analysis\n")
## }

##mgd=checkIsNa(mgd, c("num.stressful.events", "num.severe.events", "total.sum.stress"), inDrop=FALSE)
##graphVariables(mgd)
results.table.filename=file.path(group.results.dir, "dti.psychometrics.results.table.csv")
## stop()

## 378 identified as trensgender, was biologically female and not on any hormone therapy
mgd[mgd$subject=="378", "Gender"]="F"
mgd$Gender=droplevels(mgd$Gender)

if (exists("results.stack") ) {
     analyse(mgd, group.variable="Grp", compute.effect.size=FALSE, results.table.filename=NULL, results.stack=results.stack)
     ## analyse(mgd, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=results.table.filename, results.stack=results.stack)
} else {
    analyse(mgd, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=NULL)
    ## analyse(mgd, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=results.table.filename)
}

write.csv(mgd, file=file.path(group.results.dir, "mgd.dataframe.for.olga.csv"), col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)


## WASI.Full.4 outliers
##     subject Grp WASI.Full.4
## 64      164 MDD         195
## 65      165 NCL         204
## 188     417 NCL         245


## > t.test(WASi.Full.4.nooutliers[WASi.Full.4.nooutliers$Grp=="MDD", "WASI.Full.4"], WASi.Full.4.nooutliers[WASi.Full.4.nooutliers$Grp=="NCL", "WASI.Full.4"])

## 	Welch Two Sample t-test

## data:  WASi.Full.4.nooutliers[WASi.Full.4.nooutliers$Grp == "MDD", "WASI.Full.4"] and WASi.Full.4.nooutliers[WASi.Full.4.nooutliers$Grp == "NCL", "WASI.Full.4"]
## t = -3.772, df = 94.933, p-value = 0.0002817
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -14.113612  -4.379995
## sample estimates:
## mean of x mean of y 
##  100.2826  109.5294 


## ****************************************************************************************************
## *** The following subjects have CDRS.tscore data greater than 54.000000
## 119 346 370 
## NCL NCL NCL 
## 59 55 61 
## ****************************************************************************************************
## ****************************************************************************************************
## *** The following subjects have CDRS.tscore data less than 55.000000
## 118 120 
## MDD MDD 
## 53 44 
## ****************************************************************************************************
## ****************************************************************************************************
## *** The following medicated subjects were found in the analysis
## 130 132 319 320 322 323 325 329 333 376 
## MDD MDD MDD MDD MDD MDD MDD MDD MDD MDD 
## ****************************************************************************************************

## pwr.2p2n.test(h=NULL,n1=80,n2=245,sig.level=0.05,alternative="greater")
