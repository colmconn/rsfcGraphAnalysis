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
            ## legend.text = element_text(angle=45),
            
            ##axis.title.x=element_blank(),
            axis.title.x = element_text(size=my.base.size, vjust=0),
            axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
            plot.title=element_text(size=my.base.size*1.2, vjust=1))

########################################################################################################################################################################################################
### START OF FUNCTIONS #################################################################################################################################################################################
########################################################################################################################################################################################################
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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

analyse.differences.between.included.and.excluded <- function(inData, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=NULL) {
    
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

     for (i in 1: length(psychMeasures)) {
        variable=psychMeasures[[i]]$variable
        name=psychMeasures[[i]]$name
        
        cat("################################################################################\n");
        cat("*** Summary for ", name, "variable = ", variable, "\n")
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
            
            if (any(is.na(inData[, variable]))) {
                cat (sprintf("*** The following subjects have missing data for %s (%s)\n", name, variable))
                cat(as.vector(inData[is.na(inData[, variable]), "subject"]), "\n")
                cat(as.vector(inData[is.na(inData[, variable]), group.variable]), "\n")      
            }

            ## cat("*** Analyzing ", variable, "\n")
            if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {
                if (isNonparametricTestVariable(variable)) {
                    ##cat("Control length: ", length(inData[inData[, group.variable]==group1, variable]), "MDD length: ", length(inData[inData[, group.variable]==group2, variable]), "\n")
                    ## if( inherits(test <- try(wilcox.test(inData[inData[, group.variable]==group1, variable],
                    ##                                      inData[inData[, group.variable]==group2, variable]),
                    ##                          silent=FALSE),
                    ##              "try-error") ) {
                    ##     test <- 0
                    ## }
                    test <- 0
                    cat(sprintf("*** The wilcox.test hasn't been programmed yet. Can't analyze %s (%s)\n", name, variable))
                } else {
                    form=sprintf("%s ~ %s*include", variable, group.variable)
                    cat(sprintf("*** Using %s as formula for ANOVA\n", form))
                    if( inherits(test <- try(aov(as.formula(form), data=inData),
                                             silent=FALSE),
                                 "try-error") ) {
                        test <- 0
                    }
                    ## test=aov(as.formula(form), data=inData)
                }
                if (is.list(test)) {
                    print(summary(test))
                    print(TukeyHSD(test))
                }

                sm.df=summarySE(inData, measure=variable, groupvars=c(group.variable), na.rm=TRUE)
                print(sm.df)
                
                sm.df=summarySE(inData, measure=variable, groupvars=c("include"), na.rm=TRUE)
                print(sm.df)
                
                sm.df=summarySE(inData, measure=variable, groupvars=c(group.variable, "include"), na.rm=TRUE)
                print(sm.df)
                
            } else {
                cat ("*** Insufficient opservations\n")
            } ## end of if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {
            
        } ## end of if (is.factor(inData[, variable])) { ... } else { ... }
    }
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

graph.rads.by.inclusion.interaction <- function(inData) {

    imageDirectory=file.path(group.results.dir, "psychMeasureGraphs")
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }
    
    y.axis="RADS.Total.Tscore"
    data.sumamry.df=summarySE(inData, measure=y.axis, groupvars=c("Grp", "include"), na.rm=TRUE)
    print(data.sumamry.df)
    
    my.dodge=position_dodge(.2)

    interaction.graph.a=ggplot(inData, aes_string(x="include", y=y.axis, group="Grp", color="Grp", shape="Grp")) +
        geom_jitter() +
            geom_errorbar(data=data.sumamry.df, aes_string(ymin=paste(y.axis, "se", sep="-"), ymax=paste(y.axis, "se", sep="+")), width=.2, position=my.dodge) +
                stat_summary(fun.y = mean, geom="point", position=my.dodge) +
                    stat_summary(fun.y=mean, geom="line", position=my.dodge) +
                        scale_color_brewer(name="Group:", palette="Set1") +
                            scale_shape_discrete(name="Group:") +                    
                                labs(x="Inclusion", y="Reynolds Adolescent Depression Scale Total (Standardized)") +
                                    my.theme
    print(interaction.graph.a)
    imageFilename=file.path(imageDirectory, "RADS-2.by.inclusion.interaction.a.pdf")
    cat(paste("*** Creating", imageFilename, "\n"))
    ggsave(imageFilename, interaction.graph.a)
    
    interaction.graph.b=ggplot(inData, aes_string(x="Grp", y=y.axis, group="include", color="include", shape="include")) +
        geom_jitter() +
            geom_errorbar(data=data.sumamry.df, aes_string(ymin=paste(y.axis, "se", sep="-"), ymax=paste(y.axis, "se", sep="+")), width=.2, position=my.dodge) +
                stat_summary(fun.y = mean, geom="point", position=my.dodge) +
                    stat_summary(fun.y=mean, geom="line", position=my.dodge) +
                        scale_color_brewer(name="Inclusion:", palette="Set1") +
                            scale_shape_discrete(name="Inclusion:") +                    
                                labs(x="Group", y="Reynolds Adolescent Depression Scale Total (Standardized)") +
                                    my.theme
    dev.new()
    print(interaction.graph.b)
    imageFilename=file.path(imageDirectory, "RADS-2.by.inclusion.interaction.b.pdf")
    cat(paste("*** Creating", imageFilename, "\n"))
    ggsave(imageFilename, interaction.graph.b)

}
    

graphVariables <- function (inData) {


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

##original.control.subjectList.filename=file.path(config.data.dir, "estop.subjList.ncl.txt")
##original.mdd.subjectList.filename=file.path(config.data.dir, "estop.subjList.mdd.txt")

original.control.subjectList.filename=file.path(config.data.dir, "rs114a.ncl.txt")
original.mdd.subjectList.filename=file.path(config.data.dir, "rs114a.mdd.txt")


####################################################################################################
### This variable controls whether the analysis should be done for the
### VBM analysis (TRUE) of the RSFC analysis (FALSE)
####################################################################################################
vbmAnalysis=FALSE

if (! vbmAnalysis) {
    cat("####################################################################################################\n")
    cat("### Setting up files for the VBM analysis for Matthew's graph paper\n")
    cat("####################################################################################################\n")
    
    ## seed.name="L_BLA.tiff.3mm"
    seed.name="L_BLA.weight.3mm"    

    ctrl.subjectOrder.filename=file.path(group.data.dir, paste("subjectOrder.ctrlOnly", seed.name, "csv", sep="."))
    mdd.subjectOrder.filename=file.path(group.data.dir, paste("subjectOrder.mddOnly",seed.name, "csv", sep="."))

### Original subject lists
    cat ("########## Reading the original subject lists\n")

    cat("Reading original CONTROL subject list: ", original.control.subjectList.filename, "\n")
    original.control.subjectList=read.table(original.control.subjectList.filename, header=FALSE)
    colnames(original.control.subjectList)=c("subject")
    cat(sprintf("Read %d subjects from the original CONTROL subject list\n", dim(original.control.subjectList)[1]))

    cat("Reading original MDD subject list: ", original.mdd.subjectList.filename, "\n")
    original.mdd.subjectList=read.table(original.mdd.subjectList.filename, header=FALSE)
    colnames(original.mdd.subjectList)=c("subject")
    cat(sprintf("Read %d subjects from the original MDD subject list\n", dim(original.mdd.subjectList)[1]))


### Subjects actually in the bucket files
    cat ("########## Reading the bucklet lists\n")

    cat("Reading list of CONTROLs subject in the bucket file: ", ctrl.subjectOrder.filename, "\n")
    ctrl.subjectList=filterOutStimulusColumn(read.csv(ctrl.subjectOrder.filename, header=TRUE))
    colnames(ctrl.subjectList)=c("subject")
    cat(sprintf("Read %d subjects from the CONTROL subject bucket list\n", dim(ctrl.subjectList)[1]))

    cat("Reading list of MDDs subject in the bucket file: ", mdd.subjectOrder.filename, "\n")
    mdd.subjectList=filterOutStimulusColumn(read.csv(mdd.subjectOrder.filename, header=TRUE))
    colnames(mdd.subjectList)=c("subject")
    cat(sprintf("Read %d subjects from the MDD subject bucket list\n", dim(mdd.subjectList)[1]))


### Subjects that were not in the bucket files, either due to error or missing regressors etc.
    cat ("########## The following subjects were in the original subjects list but were automatically excluded from bucket generation due to missing data\n")

    nodata.ctrl.filename=file.path(group.data.dir, paste("nodata.ctrlOnly", seed.name, "txt", sep="."))
    cat("Reading missing subject information file for CONTROLs: ", nodata.ctrl.filename, "\n")
    if (file.info(nodata.ctrl.filename)$size > 0 ) {
        nodata.ctrl=filterOutStimulusColumn(read.csv(nodata.ctrl.filename, header=FALSE))
        colnames(nodata.ctrl)=c("subject")
    } else {
        nodata.ctrl=data.frame("subject"=c())
    }
    ## nodata.ctrl$subject=as.factor(gsub("_A", "", as.character(nodata.ctrl$subject), fixed=TRUE))
    cat(sprintf("Read missing subject information for %s unique CONTROL subjects\n",  length(unique(nodata.ctrl$subject))))
    cat("There was no data for the following CTRL subjects:\n")
    cat(as.vector(nodata.ctrl$subject))
    cat("\n")

    nodata.mdd.filename=file.path(group.data.dir, paste("nodata.mddOnly", seed.name, "txt", sep="."))
    cat("Reading missing subject information file for MDDs: ", nodata.mdd.filename, "\n")
    if (file.info(nodata.mdd.filename)$size > 0 ) {
        nodata.mdd=filterOutStimulusColumn(read.csv(nodata.mdd.filename, header=FALSE))
        colnames(nodata.mdd)=c("subject")
    } else {
        nodata.mdd=data.frame("subject"=c())
    }
    ##nodata.mdd$subject=as.factor(gsub("_A", "", as.character(nodata.mdd$subject), fixed=TRUE))
    cat(sprintf("Read missing subject information for %s unique MDD subjects\n",  length(unique(nodata.mdd$subject))))
    cat("There was no data for the following MDD subjects:\n")
    cat(as.vector(nodata.mdd$subject))
    cat("\n")


### Subjects that were not in the bucket files, either due to too much motion.
    cat ("########## The following subjects were in the original subjects list but were automatically excluded from bucket generation due to much motion\n")
    doNotAnalyze.ctrl.filename=file.path(group.data.dir, paste("doNotAnalyse.ctrlOnly", seed.name, "txt", sep="."))
    cat("Reading missing subject information file for CONTROLs: ", doNotAnalyze.ctrl.filename, "\n")
    if (file.info(doNotAnalyze.ctrl.filename)$size > 0 ) {
        doNotAnalyze.ctrl=filterOutStimulusColumn(read.csv(doNotAnalyze.ctrl.filename, header=FALSE))
        colnames(doNotAnalyze.ctrl)=c("subject")
    } else {
        doNotAnalyze.ctrl=data.frame("subject"=c())
    }
    ## doNotAnalyze.ctrl$subject=as.factor(gsub("_A", "", as.character(doNotAnalyze.ctrl$subject), fixed=TRUE))
    cat(sprintf("Read missing subject information for %s unique CONTROL subjects\n",  length(unique(doNotAnalyze.ctrl$subject))))
    cat("There was no data for the following CTRL subjects:\n")
    cat(as.vector(doNotAnalyze.ctrl$subject))
    cat("\n")

    doNotAnalyze.mdd.filename=file.path(group.data.dir, paste("doNotAnalyse.mddOnly", seed.name, "txt", sep="."))
    cat("Reading missing subject information file for MDDs: ", doNotAnalyze.mdd.filename, "\n")
    if (file.info(doNotAnalyze.mdd.filename)$size > 0 ) {
        doNotAnalyze.mdd=filterOutStimulusColumn(read.csv(doNotAnalyze.mdd.filename, header=FALSE))
        colnames(doNotAnalyze.mdd)=c("subject")
    } else {
        doNotAnalyze.mdd=data.frame("subject"=c())
    }
    ##doNotAnalyze.mdd$subject=as.factor(gsub("_A", "", as.character(doNotAnalyze.mdd$subject), fixed=TRUE))
    cat(sprintf("Read missing subject information for %s unique MDD subjects\n",  length(unique(doNotAnalyze.mdd$subject))))
    cat("There was no data for the following MDD subjects:\n")
    cat(as.vector(doNotAnalyze.mdd$subject))
    cat("\n")

} else {
    ## analysis of the VBM subjects

    group.vbm.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/vbm.subject.list.from.matthew.n114"))
    
    ctrl.subjectOrder.filename=file.path(group.vbm.dir, "ncl.subjectlist.txt")
    mdd.subjectOrder.filename=file.path(group.vbm.dir,  "mdd.subjectlist.txt")

    cat("Reading list of CONTROLs subject in the bucket file: ", ctrl.subjectOrder.filename, "\n")
    ctrl.subjectList=read.csv(ctrl.subjectOrder.filename, header=FALSE)
    colnames(ctrl.subjectList)=c("subject")
    cat(sprintf("Read %d subjects from the CONTROL subject bucket list\n", dim(ctrl.subjectList)[1]))
    
    cat("Reading list of MDDs subject in the bucket file: ", mdd.subjectOrder.filename, "\n")
    mdd.subjectList=read.csv(mdd.subjectOrder.filename, header=FALSE)
    colnames(mdd.subjectList)=c("subject")
    cat(sprintf("Read %d subjects from the MDD subject bucket list\n", dim(mdd.subjectList)[1]))

}


## Create the subjectOrder variable used throughout the rest of the script
## the clean subject list has no timepoint suffix
subjectOrder=data.frame("subject"=c(as.vector(mdd.subjectList$subject), as.vector(ctrl.subjectList$subject)))
colnames(subjectOrder)=c("subject")
subjectOrder$subject=as.vector(subjectOrder$subject)
## due to some sort of screw up 169/300_A is the same subject so fix it here.
badSubjectIds=which (subjectOrder$subject=="300_A")
if (length(badSubjectIds) > 0 ) {
    subjectOrder[badSubjectIds, "subject"]="169/300_A"
}
subjectOrder$subject=as.factor(gsub("_A[0-9]?", "", as.character(subjectOrder$subject), fixed=FALSE))
cat(sprintf("*** The subjectOrder data frame has %s unique subjects\n",  length(unique(subjectOrder$subject))))

##stop("Stooooooooooooooooooooooooooooooooooooooooooooooooooping")


## due to some sort of screw up 169/300_A is the same subject so fix it here.
mdd.subjectList$subject=as.vector(mdd.subjectList$subject)

badSubjectIds=which (mdd.subjectList$subject=="300_A")
if (length(badSubjectIds) > 0 ) {
    mdd.subjectList[badSubjectIds, "subject"]="169/300_A"
}
##mdd.subjectList[which (mdd.subjectList$subject=="300_A"), "subject"]="169/300_A"

mdd.subjectList$subject=as.factor(gsub("_A[0-9]?", "", as.character(mdd.subjectList$subject), fixed=FALSE))

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
    "MADRS.raw",                "PSWQ", 	             "CoRum",
    "RADS.DM",                  "RADS.AN",                   "RADS.NS",                   "RADS.SC",
    "RADS.DM.Tscore",           "RADS.AN.Tscore",            "RADS.NS.Tscore",            "RADS.SC.Tscore",  	"RADS.Total.Tscore",
    "BDI.II",                   "CDI",                       "MASC.total",                "MASC.tscore",        "RSQ.fixed",
    "SCARED.som.panic",         "SCARED.gen.anx",            "SCARED.seper.anx",          "C_Irritability",	"C_EmoSuscept",
    "SCARED.soc.phobia",        "SCARED.school.phobia",      "BIS",                       "BAS.drive",                
    "BAS.funseek",              "BAS.reward.responsiveness", "PsycAche.total", "NS",	"HA",	"RD",	"P",	"SD",	"C",	"ST")

## 2015/08/13: Don't use these from the demographics sheet
## "WASI.PERF", "WASI.Full.4",

if (length(setdiff(selectedColumns, colnames(demographics))) != 0) {
    message ("*** The following column name(s) is/are causing problems: ", setdiff(selectedColumns, colnames(demographics)), "\n")
    stop("There is a mismatch between the selected columns and the columns in the demographics data frame. Stopping.")
}

####################################################################################################
### Set up the data frames to work out how many were excluded for motion and missing data
####################################################################################################

if (! vbmAnalysis) {

    results.stack <- stack()
    
    original.subjectList=rbind(original.control.subjectList, original.mdd.subjectList)
    original.subjectList$subject=gsub("_A[0-9]?", "", as.character(original.subjectList$subject), fixed=FALSE)
    ## due to some sort of screw up 169/300_A is the same subject so fix it here.
    original.subjectList[which (original.subjectList$subject=="300"), ]="169/300"
    original.subjectList.mgd=cbind(original.subjectList, demographics[match(original.subjectList$subject, demographics$ID), c("Grp", "Gender")])
    original.subjectList.mgd$subject=as.factor(original.subjectList.mgd$subject)
    original.subjectList.mgd$Grp=drop.levels(original.subjectList.mgd$Grp)
    original.subjectList.mgd$Gender=drop.levels(original.subjectList.mgd$Gender)

    nodata.subjectList=rbind(nodata.ctrl, nodata.mdd)
    nodata.subjectList$subject=gsub("_A[0-9]?", "", as.character(nodata.subjectList$subject), fixed=FALSE)
    ## due to some sort of screw up 169/300_A is the same subject so fix it here.
    if (length(which (nodata.subjectList$subject=="300")) > 0 ) {
        nodata.subjectList[which (nodata.subjectList$subject=="300"), ]="169/300"
    }
    nodata.subjectList.mgd=cbind(nodata.subjectList, demographics[match(nodata.subjectList$subject, demographics$ID), c("Grp", "Gender")])
    nodata.subjectList.mgd$subject=as.factor(nodata.subjectList.mgd$subject)
    nodata.subjectList.mgd$Grp=drop.levels(nodata.subjectList.mgd$Grp)
    nodata.subjectList.mgd$Gender=drop.levels(nodata.subjectList.mgd$Gender)


    doNotAnalyze.subjectList=rbind(doNotAnalyze.ctrl, doNotAnalyze.mdd)
    doNotAnalyze.subjectList$subject=gsub("_A[0-9]?", "", as.character(doNotAnalyze.subjectList$subject), fixed=FALSE)
    ## due to some sort of screw up 169/300_A is the same subject so fix it here.
    if (length(which (doNotAnalyze.subjectList$subject=="300")) > 0 ) {
        doNotAnalyze.subjectList[which (doNotAnalyze.subjectList$subject=="300"), ]="169/300"
    }
    ## doNotAnalyze.subjectList.mgd=cbind(doNotAnalyze.subjectList, demographics[match(doNotAnalyze.subjectList$subject, demographics$ID), c("Grp", "Gender")])
    doNotAnalyze.subjectList.mgd=cbind(doNotAnalyze.subjectList, demographics[match(doNotAnalyze.subjectList$subject, demographics$ID), selectedColumns])    
    doNotAnalyze.subjectList.mgd$subject=as.factor(doNotAnalyze.subjectList.mgd$subject)
    doNotAnalyze.subjectList.mgd$Grp=drop.levels(doNotAnalyze.subjectList.mgd$Grp)
    doNotAnalyze.subjectList.mgd$Gender=drop.levels(doNotAnalyze.subjectList.mgd$Gender)

    cat("\n*** Table of individuals originally included in analysis\n")
    original.gender.table=table(original.subjectList.mgd[, c("Gender", "Grp")])
    ##original.gender.test=prop.test(original.gender.table)
    print(addmargins(original.gender.table))

    original.n.table=table(original.subjectList.mgd[, c("Grp")])
    original.n.test=prop.test(original.n.table)
    print(addmargins(original.n.table))
    print(original.n.test)
    csvLine=make.proportions.test.string.for.results.stack(original.n.table, original.n.test, "MDD", "NCL", "Total number of participants recruited (n)")
    push(results.stack, csvLine)
    
    cat("\n*** Table of individuals dropped because of missing data\n")
    nodata.gender.table=table(nodata.subjectList.mgd[, c("Gender", "Grp")])
    ##nodata.gender.test=prop.test(nodata.gender.table)
    print(addmargins(nodata.gender.table))

    nodata.n.table=table(nodata.subjectList.mgd[, c("Grp")])
    nodata.n.test=prop.test(nodata.n.table)
    nodata.n.table=nodata.n.table
    print(addmargins(nodata.n.table))
    print(nodata.n.test)
    csvLine=make.proportions.test.string.for.results.stack(nodata.n.table, nodata.n.test, "MDD", "NCL", "Total number of participants dropped becasue of missing data (n)")
    push(results.stack, csvLine)
    

    if (dim(doNotAnalyze.subjectList.mgd)[1] > 0) {
        cat("\n*** Table of originally included in analysis less individuals dropped because of missing data\n")
        total.subjects.less.dropped.for.missing.data.table = original.n.table - nodata.n.table
        print(addmargins(total.subjects.less.dropped.for.missing.data.table))
        small.n.test=prop.test(total.subjects.less.dropped.for.missing.data.table)
        print(small.n.test)
        csvLine=make.proportions.test.string.for.results.stack(total.subjects.less.dropped.for.missing.data.table, small.n.test, "MDD", "NCL",
            "Total number of participants in final analysis less those dropped becasue of missing data (n)")
        push(results.stack, csvLine)
        
        cat("\n*** Table of individuals dropped because of excessive motion\n")
        doNotAnalyze.gender.table=table(doNotAnalyze.subjectList.mgd[, c("Gender", "Grp")])
        ##doNotAnalyze.gender.test=prop.test(doNotAnalyze.gender.table)
        doNotAnalyze.gender.table=addmargins(doNotAnalyze.gender.table)
        print(doNotAnalyze.gender.table)

        doNotAnalyze.n.table=table(doNotAnalyze.subjectList.mgd[, c("Grp")])
        doNotAnalyze.n.test=prop.test(doNotAnalyze.n.table)
        print(addmargins(doNotAnalyze.n.table))
        print(doNotAnalyze.n.test)
        csvLine=make.proportions.test.string.for.results.stack(doNotAnalyze.n.table, doNotAnalyze.n.test, "MDD", "NCL",
            "Total number of participants dropped because of excessive motion (n)")
        push(results.stack, csvLine)
    }
}

####################################################################################################
### Set up the data frames for statistical analysis
####################################################################################################

mgd=cbind(
    subjectOrder,
    demographics    [match(subjectOrder$subject, demographics$ID), selectedColumns],
    wasi            [match(subjectOrder$subject, wasi$SubID), c("Verbal", "Performance", "Full")],
    ctq             [match(subjectOrder$subject, ctq$ID), c("CTQ_TOTAL", "CTQ_EA", "CTQ_PA", "CTQ_SA", "CTQ_EN", "CTQ_PN", "CTQ_MD")])

doNotAnalyze.mgd=cbind(
    doNotAnalyze.subjectList.mgd,
    wasi            [match(doNotAnalyze.subjectList.mgd$subject, wasi$SubID), c("Verbal", "Performance", "Full")],
    ctq             [match(doNotAnalyze.subjectList.mgd$subject, ctq$ID), c("CTQ_TOTAL", "CTQ_EA", "CTQ_PA", "CTQ_SA", "CTQ_EN", "CTQ_PN", "CTQ_MD")])


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

mgd$include=rep("include", dim(mgd)[1])


## add this column with the subject ID concatenated with the timepoint
## so we can pull the correct data from the SLES data frame. This is
## done because the SLES data frame contains data for all time points.
doNotAnalyze.mgd$subjectWithTimePoint=as.factor(paste(doNotAnalyze.mgd$subject, "A", sep="_"))
doNotAnalyze.mgd=cbind(
    doNotAnalyze.mgd,
    sles [match(doNotAnalyze.mgd$subjectWithTimePoint, sles$subject),
          c("num.stressful.events", "num.severe.events", "total.sum.stress", "num.1.little.to.no.effect", "num.2.some.effect", "num.3.moderate.effect", "num.4.great.effect")])
## remove the now unnecessary subject with timepoint column
doNotAnalyze.mgd$subjectWithTimePoint=NULL

doNotAnalyze.mgd$Grp=drop.levels(doNotAnalyze.mgd$Grp)
doNotAnalyze.mgd$SCARED.som.panic=as.numeric(as.vector(doNotAnalyze.mgd$SCARED.som.panic))
doNotAnalyze.mgd$MASC.tscore=as.numeric(doNotAnalyze.mgd$MASC.tscore)

doNotAnalyze.mgd=fixDates(doNotAnalyze.mgd)
doNotAnalyze.mgd=computeAge(doNotAnalyze.mgd)

doNotAnalyze.mgd=computeMascScore(doNotAnalyze.mgd)
rownames(doNotAnalyze.mgd)=NULL

doNotAnalyze.mgd$include=rep("exclude", dim(doNotAnalyze.mgd)[1])

all.subjects=rbind(mgd, doNotAnalyze.mgd)
all.subjects$include=as.factor(all.subjects$include)

## stop("stopppppppppping")

psychMeasures=list(
    list(variable="age.in.years", name="Age at time of scan (years)"),
    ## list(variable="Hand",         name="Edinburgh Handedness Inventory"),
    list(variable="SES",          name="Hollingshead Socioeconomic Score"),
    list(variable="TannerAvg",    name="Tanner Score"),
    ## list(variable="Tanner2",    name="Tanner Score"),    
    
    ## use these for the demographics table
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
    
    ## list(variable="RADS.DM.Tscore",      name="Reynolds Adolescent Depression Scale Dysphoric Mood (Standardized)"),
    ## list(variable="RADS.AN.Tscore",      name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect (Standardized)"),
    ## list(variable="RADS.NS.Tscore",      name="Reynolds Adolescent Depression Scale Negative Self-evaluation (Standardized)"),
    ## list(variable="RADS.SC.Tscore",      name="Reynolds Adolescent Depression Scale Somatic Complaints (Standardized)"),
     list(variable="RADS.Total.Tscore",      name="Reynolds Adolescent Depression Scale Total (Standardized)"),
    
    
    ## list(variable="SFMQ.total",     name="SFMQ"),
    ## list(variable="MADRS.raw",      name="Montgomery-Asberg Depression Rating Scale"),
    ## list(variable="BDI.II",         name="Beck Depression Inventory II"),
    ## list(variable="CDI",            name="Children's Depression Inventory"),
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
    warning("Cannot continue with those control subjects included in the analysis\n")
    depressedControls=TRUE
}

if(checkColumnAgainstThreshold(mgd[mgd$Grp=="MDD", ], "CDRS.tscore", 55, inDirection="l")) {
    warning("Cannot continue with those MDD subjects included in the analysis\n")
    underdepressed=TRUE
}

if(checkForMedicatedSubjects(mgd)) {
    warning("Cannot continue with medicated subjects included in the analysis\n")
    medicatedSubjects=TRUE
}

## stopifnot(!depressedControls, !underdepressed, !medicatedSubjects)


## if (checkForNonMaleOrFemale(mgd)) {
##     warning("Cannot continue with subjects whose gender is neither male or female included in the analysis\n")
## }

##mgd=checkIsNa(mgd, c("num.stressful.events", "num.severe.events", "total.sum.stress"), inDrop=FALSE)
##graphVariables(mgd)
results.table.filename=file.path(group.results.dir, "psychometrics.results.table.csv")
## stop()
## analyse(mgd, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=results.table.filename)

## 378 identified as trensgender, was biologically female and not on any hormone therapy
mgd[mgd$subject=="378", "Gender"]="F"
mgd$Gender=droplevels(mgd$Gender)

## if (exists("results.stack") ) {
##     analyse(mgd, group.variable="Grp", compute.effect.size=FALSE, results.table.filename=NULL, results.stack=results.stack)
##     ## analyse(mgd, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=results.table.filename, results.stack=results.stack)
## } else {
##     analyse(mgd, group.variable="Grp", compute.effect.size=FALSE, results.table.filename=NULL)
##     #analyse(mgd, group.variable="Grp", compute.effect.size=TRUE, results.table.filename=results.table.filename)
## }


analyse.differences.between.included.and.excluded(all.subjects, group.variable="Grp", compute.effect.size=FALSE, results.table.filename=NULL)
graph.rads.by.inclusion.interaction(all.subjects)

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
