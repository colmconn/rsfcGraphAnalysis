rm(list=ls())
graphics.off()

sink.reset <- function(){
    for(i in seq_len(sink.number())){
        sink(NULL)
    }
}

sink.reset()

## for the some function
library(car)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
library(grid)
library(plyr)
library(compute.es)
##source('pcor.R')

####################################################################################################
### Start of functions
####################################################################################################

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

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s,1,1)),
                             {s <- substring(s,2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

make.significance.indications <- function(pValues, which.pValues=c(1)) {

    Signif=symnum(pValues, corr=FALSE, na=FALSE, cutpoints = c(0,  .001,.01, .05, .1, 1),
        symbols   =  c("***", "**", "*", ".", " "))
    f=format(Signif)

    ## only return the first one as we're only interested in marking significant group effects
    return(f[which.pValues])
}

substituteShortLabels <- function(inLevel) {
    returnSubstitutedLabel = gsub("[0-9]+ ", "", gsub("Inf", "Inferior", gsub("Sup", "Superior", gsub("Gy", "Gyrus",
        gsub("^R", "Right", gsub("^L", "Left", inLevel, fixed=FALSE), fixed=FALSE), fixed=TRUE), fixed=TRUE), fixed=TRUE))
    
    return (returnSubstitutedLabel)
}

stderror <- function(x) sd(x)/sqrt(length(x))

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

    ## New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    ## This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                       c( N    = length2(xx[,col], na.rm=na.rm),
                         mean = mean   (xx[,col], na.rm=na.rm),
                         sd   = sd     (xx[,col], na.rm=na.rm)
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

makePublicationTable <- function(inClusterWhereAmI, inClusters, inRoistats, inRoistats.averageStatValue=NULL, inRoistats.averageContrastValue=NULL,
                                 inSummaryColumns=NULL,
                                 inStatColumnName="Default Stat Name", inContrastColumnName="Default Contrast Name", in.dof=NULL, inCom=TRUE) {
  hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
  ##print(hemisphere)
  if ( inCom ) {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
  } else {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
  }
  publication.table.header=c("Structure", colnames(locations)[-1])
  
  ## cat("Locations: Volume and coordinates\n")
  ## print(locations)

  ## cat ("Columns matching Mean: ", grep("Mean", colnames(inRoistats)), "\n")
  ## cat ("Data from the above columns:\n")
  ## print(inRoistats[, grep("Mean", colnames(inRoistats))])
  ## print(list(inRoistats$Group))
  ## agg=aggregate(inRoistats[, grep("Mean", colnames(inRoistats))], list(inRoistats$Group), mean)

  ##ddply.agg=ddply(inRoistats, .(timepoint, Group), summarise, mean)


  ## now make a summary table with the mean of each timepoint for each
  ## group. There will be one row for each timepoint x group with the
  ## means for the ROIs for each timepoint x group occupying the
  ## remaining columns
  ddply.agg=ddply(inRoistats, inSummaryColumns,
      .fun=colwise(
          .fun=function (xx) {
              c(mean=mean(xx))
          },
          ## which columns to apply the function to are listed below
          colnames(inRoistats)[grep("Mean", colnames(inRoistats))]
          )
      )

  
  ## cat("ddply.agg:\n")
  ## print(ddply.agg)
  
  ## cat("t(ddply.agg):\n")
  ## print(t(ddply.agg))

  ## now transpose the means so that there are as many rows as
  ## columns. This is done so that it can be cbinded with the ROI
  ## center of mass and volumes later
  ddply.agg.means=t(ddply.agg[grep("Mean", colnames(ddply.agg))])

  ## print ("*** Column numbers that are not Mean")
  ## print(which (! grepl("Mean", colnames(ddply.agg))))
  ## print ("*** Column names that are not Mean")  
  ## print(ddply.agg[, which (! grepl("Mean", colnames(ddply.agg)))])
  ## print("*** cnames computation")
  
  ## now make the column names
  ##
  ## the first branch handles more than two summary variables, e.g.,
  ## Group X timepoint or Group X Gender
  ##
  ## the econd branch handles only single variables as with a main
  ## effect, e.g., Group, or Gender
  if (length(which (! grepl("Mean", colnames(ddply.agg)))) > 1)
      cnames=apply(ddply.agg[, which (! grepl("Mean", colnames(ddply.agg)))], 1,
          function(xx) {
              ## xx[1] is group, xx[2] is timepoint
              return(sprintf("%s (%s)", xx[1], xx[2]))
          })
  else {
      cnames=as.character(ddply.agg[, which (! grepl("Mean", colnames(ddply.agg)))])
  }
  ## cat("cnames:\n")
  ## print(cnames)

  colnames(ddply.agg.means)=cnames
  ## cat("ddply.agg.means:\n")
  ## print(ddply.agg.means)

  mns=round(ddply.agg.means, 2)
  
  ##cat("agg: mean for each group in each ROI\n")  
  ##print(agg)
  ##mns=round(t(agg[, -1]), 2)
  ##cat("mns: transposed mean for each group in each ROI\n")    
  ##colnames(mns)=levels(agg[,1])
  ##print(mns)
  
  ## cat("inRoistats.averageStatValue\n")
  ## print (inRoistats.averageStatValue)
  if (! is.null(inRoistats.averageStatValue) ) {
      ## cat("Adding average t stats\n")
      pubTable=cbind(locations, round(t(inRoistats.averageStatValue), 2))
      ## print(pubTable)

      if ( ! is.null(in.dof)) {
          ## print((inRoistats.averageStatValue * 2) / sqrt(in.dof))
          ## print("before")
          ## print(pubTable)
          ## ##stop()
          pubTable=cbind(pubTable, round(t((inRoistats.averageStatValue * 2) / sqrt(in.dof)), 3))
          ## print("after")
          ## print(pubTable)
          
          publication.table.header=c(publication.table.header, "Cohen's d")
      }
  } else {
      pubTable=cbind(locations, mns)
  }

  if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageContrastValue) ) {
      ## cat("Adding average coefficient values\n")      
      pubTable=cbind(pubTable, round(t(inRoistats.averageContrastValue), 2))
      ## print(pubTable)      
  }

  pubTable=cbind(pubTable, mns)
  if (! is.null(inRoistats.averageStatValue) ) {
      colnames(pubTable)=
          c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", inStatColumnName, colnames(mns))
  } else {
      colnames(pubTable)=c("Structure", colnames(pubTable)[-1])
  }

  if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageContrastValue) ) {
      colnames(pubTable)=
          c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", inStatColumnName, inContrastColumnName, colnames(mns))
  }
  
  rownames(pubTable)=NULL
  ## print(pubTable)
  return(pubTable)
}

savePublicationTable <- function (inPublicationTable, inPublicationTableFilename, append=TRUE) {
    cat("*** Writing publication table to", inPublicationTableFilename, "\n")
    if ( append ) {
        write.table(inPublicationTable, file=inPublicationTableFilename, quote=F, col.names=TRUE, row.names=FALSE, sep=",", append=TRUE)
    } else {
        write.table(inPublicationTable, file=inPublicationTableFilename, quote=F, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
    }
    cat("\n", file=inPublicationTableFilename, append=TRUE)        
}

readStatsTable <- function (inFilename) {

    cat("*** Reading" , inFilename, "\n")
    statsTable=read.table(inFilename, header=T, sep="")
    ## dump the first column as it's only the file name
    statsTable=statsTable[, -1]
    return(statsTable)
}

readClustersTable <- function (inFilename){
    cat("*** Reading", file.path(inFilename), "\n")
    clusters=read.table(file.path(inFilename))
    colnames(clusters) = clust.header
    return (clusters)
}

readClusterLocationsTable <- function (inFilename) {
    cat("*** Reading cluster locations from", inFilename, "\n")
    ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
    clusterWhereAmI=gsub(" +", " ", scan(file=inFilename, what='character', sep=','))

    return (clusterWhereAmI)
}

readDataTable <- function (inFilename) {
    cat("*** Reading", inFilename, "\n")
    dataTable=read.table(inFilename, header=T, allowEscapes=TRUE)

    return(dataTable)
}

fixDataTable <- function (inDataTable) {
    inDataTable$Subj=gsub("S\\.([0-9]{3})", "\\1", as.character(inDataTable$Subj), fixed=FALSE)
    if (any(grepl("300", inDataTable$Subj, fixed=TRUE))) {
        inDataTable$subject=gsub("300", "169/300", as.character(inDataTable$Subj), fixed=TRUE)
        inDataTable$subject=as.factor(inDataTable$Subj)
    }
    return(inDataTable)
}

splitDataTableIntoIdAndTimepoint <- function(inDataTable) {

    new.subject.order.table=cbind(inDataTable, as.data.frame(t(as.data.frame(strsplit(as.character(inDataTable$subject), "_", fixed=TRUE)))))
    rownames(new.subject.order.table)=NULL
    colnames(new.subject.order.table)=c("id", "subject", "timepoint")

    return(new.subject.order.table)
}

## Reads the seed file and does the (crude) equivalent of BAS variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character())
    table=gsub("$DATA", seeds.data.dir, table, fixed=TRUE)

    return (table)
}


readDegreesOfFreedom <- function(inFilename) {
    cat("*** Reading", inFilename, "\n")
    dof=scan(inFilename, what=integer(), quiet=TRUE)

    return(dof)
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

is.f.stat <- function (inStatLabel) {
    return(grepl("F$", inStatLabel, fixed=FALSE))
}

is.z.stat <- function (inStatLabel) {
    return(grepl("Z$", inStatLabel, fixed=FALSE))
}

makeGraphTitle <- function(inLevel) {

    ## graphTitle= gsub("^[0-9]+[ ]*", "", gsub("Inf", "Inferior", gsub("Gy", "Gyrus", gsub("Sup", "Superior", inLevel))))
    graphTitle= gsub("Inf", "Inferior", gsub("Gy", "Gyrus", gsub("Sup", "Superior", inLevel)))    

    return(graphTitle)
}

run.correlations <- function (in.stats.labels, in.results.stack, in.method="pearson") {

    for (seed in seeds) {
        seedName=getSeedName(seed)

        for (statLabel in in.stats.labels) {
            
            cat("####################################################################################################\n")
            cat(sprintf("*** Graphing ROIs from the %s seed for the %s groups\n", seedName, groups))
            
            infix=sprintf("fwhm%s.%s.%s.%s.%s.%s", usedFwhm, task, groups, analysis, seedName, statLabel)

            roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))            
            if (is.f.stat(statLabel) ) {
                roistats.averageFvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageFvalue.txt", infix))
            } else if (is.z.stat(statLabel) ) {
                roistats.averageZvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageZValue.txt", infix))
                roistats.averageContrastValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageContrastValue.txt", infix))
            }
            
            if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the contrast in each ROI,
            ## you do not what to graph this

                roistats=readStatsTable(roistats.filename)
                ## print(roistats)
                ## print(class(roistats))
                ## stop()
                roistats$Sub.brick=NULL
                if (is.f.stat(statLabel) ) {
                    roistats.averageFvalue=readStatsTable(roistats.averageFvalue.filename)
                    roistats.averageFvalue$Sub.brick=NULL                
                } else if (is.z.stat(statLabel) ) {
                    roistats.averageZvalue=readStatsTable(roistats.averageZvalue.filename)
                    roistats.averageContrastValue=readStatsTable(roistats.averageContrastValue.filename)
                    roistats.averageZvalue$Sub.brick=NULL
                    roistats.averageContrastValue$Sub.brick=NULL                
                }
            
                clusterCount=length(grep("Mean", colnames(roistats)))
                if (clusterCount > 0 ) {
                    cat(sprintf("*** %d clusters found in %s\n", clusterCount, roistats.filename))
                
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
                
                    clustersFilename=file.path(group.results.dir, sprintf("clust.%s.txt", infix))
                    clusters=readClustersTable(clustersFilename)
                
                    ## this table contains the locations, as text, of the clusters and is the output of a perl script
                    clusterLocationsFilename=file.path(group.results.dir, sprintf("clusterLocations.%s.csv", infix))
                    clusterWhereAmI=readClusterLocationsTable(clusterLocationsFilename)
                
                    ## this file stores the order of the subjects in each of the following BRIK files
                    ##                     dataTable.MDD.and.NCL.L_whole_amygdala.3mm.any.timepoint.txt
                    dataTableFilename=file.path(group.data.dir, paste("dataTable.rsfc", groups, seedName, analysis, "txt", sep="."))
                    dataTable=fixDataTable(readDataTable(dataTableFilename))
                    ## print(dataTable)
                    ## print(head(demographics))
                    
                    ## stop()
                    cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(dataTable$Subj))))

                    for (cor.var.itr in seq_along(regressionVariables))  {

                        regression.variable=regressionVariables[[cor.var.itr]]$variable
                        regression.name=regressionVariables[[cor.var.itr]]$name                        
                        cat("*** Running correlations for on", regression.variable, "\n")

                        ## The use of the match and interaction comes
                        ## from:
                        ## http://stackoverflow.com/questions/6880450/match-two-columns-with-two-other-columns
                        ## this could probably be done more
                        ## efficiently using the data tables instead
                        ## of data frames
                        
                        mgd=cbind(dataTable, roistats, demographics[ match(
                                                           interaction( dataTable$Subj, dataTable$timepoint),
                                                           interaction( demographics$ID, demographics$timepoint)
                                                       ),
                                                      c("Gender", regression.variable)]) 

                        
                        ## mgd=cbind(dataTable, roistats, demographics[match(dataTable$Subj, demographics$ID), c("Gender")])
                        ## colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], "Gender")
                        ## print(mgd)
                        ## stop()
                        if (any(mgd$Subj=="378") ) {
                            mgd[mgd$Subj=="378", "Gender"]="F"
                        }
                        rownames(mgd)=NULL
                        ## ensure that subject is s factor
                        mgd$Subj=as.factor(mgd$Subj)
                        mgd=droplevels(mgd)
                        
                        ## print(clusterWhereAmI)
                        ## print(clusters)
                        ## print(roistats)
                        ## print(mgd)
                        ##cat("*** Some of the mgd data frame\n")
                        ##print(some(mgd))
                        ## stop("Check the mgd data frame\n")

                        if (is.f.stat(statLabel) ) {
                            for (mean.col.itr in seq_len(clusterCount)) {

                                mean.column = paste("Mean", mean.col.itr, sep="_")
                                
                                n.timepoints=length(levels(mgd$timepoint))
                                if (n.timepoints > 2) {
                                    cat("*** Number of timepoints (", n.timepoints, ") is more than 2 timepoints\n")
                                    mgd.dim.before.complete.cases=dim(mgd)
                                    mgd.complete=mgd[complete.cases(mgd), ]
                                    
                                    
                                    reg.formula=as.formula(paste(regression.variable, "~", mean.column, "*timepoint"))
                                    mdl = lm(reg.formula, data=mgd.complete)
                                    smry=summary(mdl)
                                    ## print(smry)
                                    ## print(df.residual(mdl))
                                    ## print(coefficients(smry))
                                    tstats=coefficients(smry)[, "t value"]
                                    ## print(
                                    p.values=2 * pt(abs(tstats), df = df.residual(mdl), lower.tail = FALSE)
                                    
                                    csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%.2f,%.2f,%.5f,%.2f,%s",
                                                    groups, analysis, seedName, statLabel,
                                                    "all", makeGraphTitle(clusterWhereAmI[mean.col.itr]),
                                                    paste(c(mean.col.itr, clusters[mean.col.itr, c("CM RL", "CM AP", "CM IS")]), collapse=","), regression.name,
                                                    mgd.dim.before.complete.cases[1], dim(mgd.complete)[1], coefficients(smry)[mean.column, "t value"],
                                                    df.residual(mdl),
                                                    p.values[mean.column], coefficients(smry)[mean.column, "Estimate"], make.significance.indications(p.values[mean.column]))
                                    
                                    ## print(csvLine)
                                    
                                    ## print(class(smry))
                                    ## print(typeof(smry))                                    
                                    ## print(smry)
                                    ## print(names(smry))
                                    ## for (nn in names(smry)) {
                                    ##     cat(nn, ": ", as.character(smry[[nn]]), "\n")
                                    ## }
                                        
                                    ## print(mgd)
                                    print(addmargins(table(mgd.complete[ , c("Group", "timepoint")])))
                                    graph.regression(mgd, seedName, statLabel, mean.col.itr, clusterWhereAmI[mean.col.itr],
                                                     mean.column,         ## x axis 
                                                     regression.variable, ## y axis
                                                     regression.name,
                                                     "RSFC (Z-score)",
                                                     regression.name)
                                    ## stop()
                                } else {
                                    mgd=mgd[order(mgd$Subj, mgd$timepoint), ]
                                    ## print(mgd)
                                    mgd.diff = cbind (mgd[mgd$timepoint=="A", c("Subj", "Group")],
                                                      mgd[mgd$timepoint=="A", mean.column]         - mgd[ mgd$timepoint=="C", mean.column],
                                                      mgd[mgd$timepoint=="A", regression.variable] - mgd[ mgd$timepoint=="C", regression.variable])
                                    colnames(mgd.diff)=c("Subj", "Group", paste(c(mean.column, regression.variable), "diff", sep="."))
                                    ## print(mgd.diff)
                                    ## stop()
                                    mgd.diff.dim.before.complete.cases=dim(mgd.diff)
                                    mgd.diff=mgd.diff[complete.cases(mgd.diff), ]

                                    ## print(mgd.diff)
                                    print(addmargins(table(mgd.diff[ , c("Group")])))

                                    ct.all=cor.test(mgd.diff[, paste(mean.column,         "diff", sep=".")],
                                                    mgd.diff[, paste(regression.variable, "diff", sep=".")],
                                                    method=in.method)
                                    ## print(ct.all)
                                    csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%.2f,%.2f,%.5f,%.2f,%s",
                                                    groups, analysis, seedName, statLabel,
                                                    "all", makeGraphTitle(clusterWhereAmI[mean.col.itr]),
                                                    paste(c(mean.col.itr, clusters[mean.col.itr, c("CM RL", "CM AP", "CM IS")]), collapse=","), regression.name,
                                                    mgd.diff.dim.before.complete.cases[1], dim(mgd.diff)[1], ct.all$statistic,
                                                    ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
                                                    ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
                                    ## print(csvLine)
                                    ## print(ct.all$p.value)
                                    ## print(make.significance.indications(ct.all$p.value))
                                    push(in.results.stack, csvLine)
                                    
                                    if (ct.all$p.value < 0.1) {
                                        graph.correlation(mgd.diff, seedName, statLabel, "all", mean.col.itr, clusterWhereAmI[mean.col.itr],
                                                          paste(mean.column,         "diff", sep="."), ## x axis 
                                                          paste(regression.variable, "diff", sep="."), ## y axis
                                                          regression.name,
                                                          expression(paste(Delta, "RSFC (Z-score)", sep=" ")),
                                                          bquote(paste(Delta, .(regression.name), sep=" ")))

                                    }
                                    
                                    ct.all=NULL
                                    for (level in levels(mgd.diff$Group)) {
                                        ##print(mgd.diff)
                                        ##print(dim(mgd.diff))
                                        mgd.diff.ss=subset(mgd.diff, Group == level)
                                        ##print(mgd.diff.ss)
                                        ##print(dim(mgd.diff.ss))
                                        ##stop()
                                        ct.all=cor.test(mgd.diff.ss[, paste(mean.column,         "diff", sep=".")],
                                                        mgd.diff.ss[, paste(regression.variable, "diff", sep=".")],
                                                        method=in.method)
                                        ## print(ct.all)

                                        ## in the csv creation line
                                        ## below the
                                        ## dim(mgd.diff.ss)[1] is
                                        ## listed twice because the
                                        ## mgd.diff.dim.before.complete.cases
                                        ## variable only referrs to
                                        ## all subjects, irrespective
                                        ## of group, going intot he
                                        ## across group correlation
                                        ## performed just above. At
                                        ## this point all incomplete
                                        ## cases have already been
                                        ## filtered out and using ot
                                        ## for per group correlations
                                        ## would be non-sensical
                                        
                                        csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%.2f,%.2f,%.5f,%.2f,%s",
                                                        groups, analysis, seedName, statLabel,
                                                        level, makeGraphTitle(clusterWhereAmI[mean.col.itr]),
                                                        paste(c(mean.col.itr, clusters[mean.col.itr, c("CM RL", "CM AP", "CM IS")]), collapse=","),
                                                        regression.name,
                                                        dim(mgd.diff.ss)[1], dim(mgd.diff.ss)[1], ct.all$statistic,
                                                        ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
                                                        ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
                                        ## print(ct.all$p.value)
                                        ## print(make.significance.indications(ct.all$p.value))
                                        ## print(csvLine)
                                        push(in.results.stack, csvLine)
                                        
                                        if (ct.all$p.value < 0.1) {
                                            graph.correlation(mgd.diff.ss, seedName, statLabel, level, mean.col.itr, clusterWhereAmI[mean.col.itr],
                                                              paste(mean.column,         "diff", sep="."), ## x axis 
                                                              paste(regression.variable, "diff", sep="."), ## y axis
                                                              regression.name,
                                                              expression(paste(Delta, "RSFC (Z-score)", sep=" ")),
                                                              bquote(paste(Delta, .(regression.name), sep=" ")))
                                            
                                        }
                                    } ## end of for (level in levels(mgd.diff$Group))
                                    ## stop()
                                } ## end of if (n.timepoints > 2)
                            } ## end of for (mean.col.itr in seq_len(clusterCount))
                        } ## end of if (is.f.stat(statLabel) )
                    } ## end of for (correlation.variable in correlation.variables ) {
                } ## end of if (clusterCount > 0 ) {
           } ## end of if(file.exists(roistats.filename)) {
        } ## for (statLabel in in.stats.labels) {
    } ## end of for (seed in seeds) {

    return(in.results.stack)
} ## end of run.correlations


run.delta.correlations <- function (in.stats.labels, in.results.stack, in.method="pearson") {

    for (seed in seeds) {
        seedName=getSeedName(seed)

        for (statLabel in in.stats.labels) {
            
            cat("####################################################################################################\n")
            cat(sprintf("*** Graphing ROIs from the %s seed for the %s groups\n", seedName, groups))
            
            infix=sprintf("fwhm%s.%s.%s.%s.%s.%s", usedFwhm, task, groups, analysis, seedName, statLabel)

            roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))            
            if (is.f.stat(statLabel) ) {
                roistats.averageFvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageFvalue.txt", infix))
            } else if (is.z.stat(statLabel) ) {
                roistats.averageZvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageZValue.txt", infix))
                roistats.averageContrastValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageContrastValue.txt", infix))
            }
            
            if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the contrast in each ROI,
            ## you do not what to graph this

                roistats=readStatsTable(roistats.filename)
                ## print(roistats)
                ## print(class(roistats))
                ## stop()
                roistats$Sub.brick=NULL
                if (is.f.stat(statLabel) ) {
                    roistats.averageFvalue=readStatsTable(roistats.averageFvalue.filename)
                    roistats.averageFvalue$Sub.brick=NULL                
                } else if (is.z.stat(statLabel) ) {
                    roistats.averageZvalue=readStatsTable(roistats.averageZvalue.filename)
                    roistats.averageContrastValue=readStatsTable(roistats.averageContrastValue.filename)
                    roistats.averageZvalue$Sub.brick=NULL
                    roistats.averageContrastValue$Sub.brick=NULL                
                }
            
                clusterCount=length(grep("Mean", colnames(roistats)))
                if (clusterCount > 0 ) {
                    cat(sprintf("*** %d clusters found in %s\n", clusterCount, roistats.filename))
                
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
                
                    clustersFilename=file.path(group.results.dir, sprintf("clust.%s.txt", infix))
                    clusters=readClustersTable(clustersFilename)
                
                    ## this table contains the locations, as text, of the clusters and is the output of a perl script
                    clusterLocationsFilename=file.path(group.results.dir, sprintf("clusterLocations.%s.csv", infix))
                    clusterWhereAmI=readClusterLocationsTable(clusterLocationsFilename)
                
                    ## this file stores the order of the subjects in each of the following BRIK files
                    ##                     dataTable.MDD.and.NCL.L_whole_amygdala.3mm.any.timepoint.txt
                    dataTableFilename=file.path(group.data.dir, paste("dataTable.rsfc", groups, seedName, analysis, "txt", sep="."))
                    dataTable=fixDataTable(readDataTable(dataTableFilename))
                    ## print(dataTable)
                    ## print(head(demographics))
                    
                    ## stop()
                    cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(dataTable$Subj))))

                    for (cor.var.itr in seq_along(regressionVariables))  {

                        regression.variable=regressionVariables[[cor.var.itr]]$variable
                        regression.name=regressionVariables[[cor.var.itr]]$name                        
                        cat("*** Running correlations for on", regression.variable, "\n")

                        ## The use of the match and interaction comes
                        ## from:
                        ## http://stackoverflow.com/questions/6880450/match-two-columns-with-two-other-columns
                        ## this could probably be done more
                        ## efficiently using the data tables instead
                        ## of data frames
                        
                        mgd=cbind(dataTable, roistats, demographics[ match(
                                                           interaction( dataTable$Subj, dataTable$timepoint),
                                                           interaction( demographics$ID, demographics$timepoint)
                                                       ),
                                                      c("Gender", regression.variable)]) 

                        
                        ## mgd=cbind(dataTable, roistats, demographics[match(dataTable$Subj, demographics$ID), c("Gender")])
                        ## colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], "Gender")
                        ## print(mgd)
                        ## stop()
                        if (any(mgd$Subj=="378") ) {
                            mgd[mgd$Subj=="378", "Gender"]="F"
                        }
                        rownames(mgd)=NULL
                        ## ensure that subject is s factor
                        mgd$Subj=as.factor(mgd$Subj)
                        mgd=droplevels(mgd)
                        
                        ## print(clusterWhereAmI)
                        ## print(clusters)
                        ## print(roistats)
                        ## print(mgd)
                        ## cat("*** Some of the mgd data frame\n")
                        ## print(mgd)
                        ## stop("Check the mgd data frame\n")
                        
                        if (is.f.stat(statLabel) ) {
                            for (mean.col.itr in seq_len(clusterCount)) {
                                
                                mean.column = paste("Mean", mean.col.itr, sep="_")
                                
                                mgd=mgd[order(mgd$Subj, mgd$timepoint), ]
                                ## print(mgd)
                                mgd.diff = cbind (mgd[mgd$timepoint=="A", c("Subj", "Group")],
                                                  mgd[mgd$timepoint=="C", mean.column]         - mgd[ mgd$timepoint=="A", mean.column],
                                                  mgd[mgd$timepoint=="C", regression.variable] - mgd[ mgd$timepoint=="A", regression.variable])
                                colnames(mgd.diff)=c("Subj", "Group", paste(c(mean.column, regression.variable), "diff", sep="."))
                                ## print(mgd.diff)
                                ## stop()
                                mgd.diff.dim.before.complete.cases=dim(mgd.diff)
                                mgd.diff=mgd.diff[complete.cases(mgd.diff), ]
                                
                                ## print(mgd.diff)
                                print(addmargins(table(mgd.diff[ , c("Group")])))
                                
                                ct.all <- cor.test(mgd.diff[, paste(mean.column,         "diff", sep=".")],
                                                   mgd.diff[, paste(regression.variable, "diff", sep=".")],
                                                   method=in.method)
                                ## print(ct.all)
                                    
                                csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%.2f,%.2f,%.5f,%.2f,%s",
                                                groups, analysis, seedName, statLabel,
                                                "all", makeGraphTitle(clusterWhereAmI[mean.col.itr]),
                                                paste(c(mean.col.itr, clusters[mean.col.itr, c("CM RL", "CM AP", "CM IS")]), collapse=","), regression.name,
                                                mgd.diff.dim.before.complete.cases[1], dim(mgd.diff)[1], ct.all$statistic,
                                                ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
                                                ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
                                ## print(csvLine)
                                ## print(ct.all$p.value)
                                ## print(make.significance.indications(ct.all$p.value))
                                push(in.results.stack, csvLine)
                                
                                ## if (ct.all$p.value < 0.1) {
                                ##     graph.correlation(mgd.diff, seedName, statLabel, "all", mean.col.itr, clusterWhereAmI[mean.col.itr],
                                ##                       paste(mean.column,         "diff", sep="."), ## x axis 
                                ##                       paste(regression.variable, "diff", sep="."), ## y axis
                                ##                       regression.name,
                                ##                       expression(paste(Delta, "RSFC (Z-score)", sep=" ")),
                                ##                       bquote(paste(Delta, .(regression.name), sep=" ")))
                                    
                                ## }
                                
                                ct.all=NULL
                                for (level in levels(mgd.diff$Group)) {
                                    ##print(mgd.diff)
                                    ##print(dim(mgd.diff))
                                    mgd.diff.ss=subset(mgd.diff, Group == level)
                                    ##print(mgd.diff.ss)
                                    ##print(dim(mgd.diff.ss))
                                    ##stop()
                                    ct.all=cor.test(mgd.diff.ss[, paste(mean.column,         "diff", sep=".")],
                                                    mgd.diff.ss[, paste(regression.variable, "diff", sep=".")],
                                                    method=in.method)
                                    ## print(ct.all)
                                    
                                    ## in the csv creation line below
                                    ## the dim(mgd.diff.ss)[1] is
                                    ## listed twice because the
                                    ## mgd.diff.dim.before.complete.cases
                                    ## variable only referrs to all
                                    ## subjects, irrespective of
                                    ## group, going intot he across
                                    ## group correlation performed
                                    ## just above. At this point all
                                    ## incomplete cases have already
                                    ## been filtered out and using ot
                                    ## for per group correlations
                                    ## would be non-sensical
                                    
                                    csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%.2f,%.2f,%.5f,%.2f,%s",
                                                    groups, analysis, seedName, statLabel,
                                                    level, makeGraphTitle(clusterWhereAmI[mean.col.itr]),
                                                    paste(c(mean.col.itr, clusters[mean.col.itr, c("CM RL", "CM AP", "CM IS")]), collapse=","),
                                                    regression.name,
                                                    dim(mgd.diff.ss)[1], dim(mgd.diff.ss)[1], ct.all$statistic,
                                                    ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
                                                    ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
                                    ## print(csvLine)
                                    ## print(ct.all$p.value)
                                    ## print(make.significance.indications(ct.all$p.value))
                                    ## print(csvLine)
                                    push(in.results.stack, csvLine)
                                    
                                    ## if (ct.all$p.value < 0.1) {
                                    ##     graph.correlation(mgd.diff.ss, seedName, statLabel, level, mean.col.itr, clusterWhereAmI[mean.col.itr],
                                    ##                       paste(mean.column,         "diff", sep="."), ## x axis 
                                    ##                       paste(regression.variable, "diff", sep="."), ## y axis
                                    ##                       regression.name,
                                    ##                       expression(paste(Delta, "RSFC (Z-score)", sep=" ")),
                                    ##                       bquote(paste(Delta, .(regression.name), sep=" ")))
                                        
                                    ## }
                                } ## end of for (level in levels(mgd.diff$Group))
                                ## stop()
                                
                            } ## end of for (mean.col.itr in seq_len(clusterCount))
                        } ## end of if (is.f.stat(statLabel) )
                    } ## end of for (correlation.variable in correlation.variables ) {
                } ## end of if (clusterCount > 0 ) {
            } ## end of if(file.exists(roistats.filename)) {
        } ## for (statLabel in in.stats.labels) {
    } ## end of for (seed in seeds) {
    
    return(in.results.stack)
} ## end of run.delta.correlations



graph.regression <-function(in.mgd, in.seed.name, in.stat.label, in.roi.no, in.cluster.name,
                            in.x.axis, in.y.axis, in.rv.name, in.x.label, in.y.label) {

    my.dodge=position_dodge(.2)
    graph=ggplot(in.mgd, aes_string(x=in.x.axis, y=in.y.axis, color="timepoint", shape="timepoint"))
    graph = graph + geom_point(position=my.dodge, size=3)
    graph = graph + geom_smooth(method="lm", se=FALSE)
    graph = graph + scale_color_brewer(name="Timepoint:", palette="Set1")
    graph = graph + scale_shape_discrete(name="Timepoint:")
    graph = graph + labs(title = substituteShortLabels(in.cluster.name), x=in.x.label, y=in.y.label) 
    graph = graph + my.theme + theme(legend.position="bottom")
    
    image.filename=file.path(group.results.dir, paste(groups, analysis, in.seed.name, in.stat.label, in.roi.no, gsub("[[:space:]]+", ".", in.cluster.name, fixed=FALSE), "pdf", sep="."))
    ## zprint(graph)

    cat("*** Saving graph to" , image.filename, "\n")
    ggsave(image.filename, graph, width=1.45, height=1.6, units="in")
    stop()
}


## rv= regression variable
graph.correlation <- function (in.mgd.diff, in.seed.name, in.stat.label, in.group, in.roi.no, in.cluster.name,
                               in.x.axis, in.y.axis, in.rv.name, in.x.label, in.y.label) {

    my.dodge=position_dodge(.2)
    graph=ggplot(in.mgd.diff, aes_string(x=in.x.axis, y=in.y.axis))
    graph = graph + geom_point(position=my.dodge)
    graph = graph + geom_smooth(method="lm", se=FALSE, color="black")
    graph = graph + labs(title = substituteShortLabels(in.cluster.name), x=in.x.label, y=in.y.label) 
    graph = graph + my.theme
    
    image.filename=file.path(group.results.dir, paste(groups, analysis, in.seed.name, in.stat.label, in.group, in.y.axis, in.roi.no, gsub("[[:space:]]+", ".", in.cluster.name, fixed=FALSE), "pdf", sep="."))
    ## print(graph)

    cat("*** Saving graph to" , image.filename, "\n")
    ggsave(image.filename, graph)
}


        
generateGraphs <- function () {

    for (seed in seeds) {
        seedName=getSeedName(seed)

        for (statLabel in statsLabels) {
            publicationTableFilename=file.path(group.results.dir, paste("publicationTable", usedFwhm, groups, analysis, seedName, "csv", sep="."))
            if (file.exists(publicationTableFilename)) {
                file.remove(publicationTableFilename)
            }
            ## cat("*** Writing publication table to", publicationTableFilename, "\n")
            
            cat("####################################################################################################\n")
            cat(sprintf("*** Graphing ROIs from the %s seed for the %s groups\n", seedName, groups))
            cat(sprintf("Analysis=%s,StatLabel=%s\n", analysis, statLabel), file=publicationTableFilename, append=TRUE)
            
            infix=sprintf("fwhm%s.%s.%s.%s.%s.%s", usedFwhm, task, groups, analysis, seedName, statLabel)

            roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))            
            if (is.f.stat(statLabel) ) {
                roistats.averageFvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageFvalue.txt", infix))
            } else if (is.z.stat(statLabel) ) {
                roistats.averageZvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageZValue.txt", infix))
                roistats.averageContrastValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageContrastValue.txt", infix))
            }
            
            degreesOfFreedom.filename=file.path(group.results.dir, sprintf("text.%s.degreesOfFreedom.txt", infix))
            if ( file.exists(degreesOfFreedom.filename) ) {
                degrees.of.freedom=readDegreesOfFreedom(degreesOfFreedom.filename)
            } else {
                cat ("*** No such file", degreesOfFreedom.filename, "\n")
                degrees.of.freedom=NULL
            }
            cat("*** Degrees of freedom set to:", paste(degrees.of.freedom, collapse=", "), "\n")
            
            if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the contrast in each ROI,
            ## you do not what to graph this

                roistats=readStatsTable(roistats.filename)
                ## print(roistats)
                ## print(class(roistats))
                ## stop()
                roistats$Sub.brick=NULL
                if (is.f.stat(statLabel) ) {
                    roistats.averageFvalue=readStatsTable(roistats.averageFvalue.filename)
                    roistats.averageFvalue$Sub.brick=NULL                
                } else if (is.z.stat(statLabel) ) {
                    roistats.averageZvalue=readStatsTable(roistats.averageZvalue.filename)
                    roistats.averageContrastValue=readStatsTable(roistats.averageContrastValue.filename)
                    roistats.averageZvalue$Sub.brick=NULL
                    roistats.averageContrastValue$Sub.brick=NULL                
                }
            
                clusterCount=length(grep("Mean", colnames(roistats)))
                if (clusterCount > 0 ) {
                    cat(sprintf("*** %d clusters found in %s\n", clusterCount, roistats.filename))
                
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
                
                    clustersFilename=file.path(group.results.dir, sprintf("clust.%s.txt", infix))
                    clusters=readClustersTable(clustersFilename)
                
                    ## this table contains the locations, as text, of the clusters and is the output of a perl script
                    clusterLocationsFilename=file.path(group.results.dir, sprintf("clusterLocations.%s.csv", infix))
                    clusterWhereAmI=readClusterLocationsTable(clusterLocationsFilename)
                
                    ## this file stores the order of the subjects in each of the following BRIK files
                    ##                     dataTable.MDD.and.NCL.L_whole_amygdala.3mm.any.timepoint.txt
                    dataTableFilename=file.path(group.data.dir, paste("dataTable.rsfc", groups, seedName, analysis, "txt", sep="."))
                    dataTable=fixDataTable(readDataTable(dataTableFilename))
                    ## print(dataTable)
                    ## stop()
                    cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(dataTable$Subj))))

                    mgd=cbind(dataTable, roistats, demographics[match(dataTable$Subj, demographics$ID), c("Gender")])
                    colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], "Gender")
                    ## print(mgd)
                    ## stop()
                    if (any(mgd$Subj=="378") ) {
                        mgd[mgd$Subj=="378", "Gender"]="F"
                    }
                    rownames(mgd)=NULL
                    ## ensure that subject is s factor
                    mgd$Subj=as.factor(mgd$Subj)
                    mgd=droplevels(mgd)
                
                    ## print(clusterWhereAmI)
                    ## print(clusters)
                    ## print(roistats)
                    ## print(mgd)

                    if (setequal(levels(mgd$timepoint), c("A", "C"))) {
                        cat ("*** Changing A and C timepoint labels to Baseline and Follow-up, respectively\n")
                        levels(mgd$timepoint) =c("Baseline", "Follow-up")
                    } else if (setequal(levels(mgd$timepoint), c("A", "C", "D"))) {
                        cat ("*** Changing A, C, and D timepoint labels to Baseline and 3 Months, and 6 Months, respectively\n")
                        levels(mgd$timepoint) =c("Baseline", "3 Months", "6 Months")
                    } else {
                        stop("Don't know how to deal with levels:",  paste(levels(mgd$timepoint), collapse=","))
                    }
                    
                    cat("*** Some of the mgd data frame\n")
                    print(some(mgd))

                    ## stop("Check the mgd data frame\n")
                    print(addmargins(table(mgd[ , c("Group", "timepoint")])))
                    print(table(mgd[ , c("Group", "Gender", "timepoint")]))                    

                    ## stop()
                    if (is.f.stat(statLabel) ) {
                        summaryColumns=unlist(strsplit(gsub(".F", "", statLabel, fixed=TRUE), "X", fixed=TRUE))
                        cat("*** The following variables will be used for the summary columns in the publication table:", paste(summaryColumns, collapse=", "), "\n")
                        publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd,
                            roistats.averageFvalue,
                            inSummaryColumns=summaryColumns,
                            inStatColumnName="Average F value",
                            inCom=TRUE)

                        cat("*** Publication table\n")
                        print(publicationTable)
                        ## savePublicationTable(publicationTable, publicationTableFilename, TRUE)
                        ## stop("Check the publication data frame\n")
                        
                        melted.mgd=melt(mgd,  id.vars=c("Subj", summaryColumns),
                            measure.vars=paste("Mean_", seq.int(1, clusterCount), sep=""),
                            variable_name="cluster")
                        
                        melted.mgd$cluster=factor(melted.mgd$cluster,
                            levels=c(paste("Mean_", seq.int(1, clusterCount), sep="")),
                            ## ensure that the sequence number is 2 digits
                            ## and padded with 0 if necessary. Not as
                            ## elegent as zip in python but achieves the
                            ## same result
                            apply(cbind(seq.int(1, length(clusterWhereAmI)), clusterWhereAmI), 1, function(xx) { sprintf("%02d %s", as.integer(xx[1]), xx[2]) } ))
                        
                        ## cat("*** Some of the melted mgd data frame\n")
                        ## print (some(melted.mgd))
                        ## stop("Check the melted mgd data frame\n")

                        graph.f.stats(melted.mgd, groups, seedName, statLabel)
                        
                    } else if (is.z.stat(statLabel) ) {
                        if (is.main.effect.contrast.z(statLabel) ) {
                            summaryColumn=get.z.contrast.terms(statLabel)
                            cat("*** The following variables will be used for the summary columns in the publication table:", summaryColumn, "\n")
                            
                            publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd,
                                roistats.averageZvalue, roistats.averageContrastValue,
                                inSummaryColumns=summaryColumn,
                                inStatColumnName="Average Z value",
                                inContrastColumnName="Average Contrast Value", inCom=TRUE)
                            cat("*** Publication table\n")
                            print(publicationTable)
                            ## savePublicationTable(publicationTable, publicationTableFilename, TRUE)
                            ## stop("Check the publication data frame\n")
                            
                            melted.mgd=melt(mgd,  id.vars=c("subject", summaryColumn),
                                measure.vars=paste("Mean_", seq.int(1, clusterCount), sep=""),
                                variable_name="cluster")
                            
                            melted.mgd$cluster=factor(melted.mgd$cluster,
                                levels=c(paste("Mean_", seq.int(1, clusterCount), sep="")),
                                ## ensure that the sequence number is 2 digits
                                ## and padded with 0 if necessary. Not as
                                ## elegent as zip in python but achieves the
                                ## same result
                                apply(cbind(seq.int(1, length(clusterWhereAmI)), clusterWhereAmI), 1, function(xx) { sprintf("%02d %s", as.integer(xx[1]), xx[2]) } ))
                            
                            cat("*** Some of the melted mgd data frame\n")
                            print (some(melted.mgd))
                            ## stop("Check the melted mgd data frame\n")
                            
                            graph.z.stats(melted.mgd, groups, seedName, statLabel, summaryColumn)
                        } else if (is.interaction.contrast.z(statLabel) ) {
                            terms=get.z.contrast.terms(statLabel)
                            summaryColumns=convert.contrast.terms.to.factor.names(terms)
                            ## terms[1] will always be the group and
                            ## the remainder will be the other levels
                            ## of the second factor
                            cat("*** The following variables will be used for the summary columns in the publication table:", paste(summaryColumns, collapse=", "), "\n")
                            
                            mgd=droplevels(subset(mgd, mgd[ , summaryColumns[1]]==terms[1]))
                            rownames(mgd)=NULL

                            publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd,
                                roistats.averageZvalue, roistats.averageContrastValue,
                                inSummaryColumns=summaryColumns,
                                inStatColumnName="Average Z value",
                                inContrastColumnName="Average Contrast Value", inCom=TRUE)
                            cat("*** Publication table\n")
                            print(publicationTable)
                            ## savePublicationTable(publicationTable, publicationTableFilename, TRUE)
                            ## stop("Check the publication data frame\n")

                            melted.mgd=melt(mgd,  id.vars=c("subject", summaryColumns[2]),
                                measure.vars=paste("Mean_", seq.int(1, clusterCount), sep=""),
                                variable_name="cluster")
                            
                            melted.mgd$cluster=factor(melted.mgd$cluster,
                                levels=c(paste("Mean_", seq.int(1, clusterCount), sep="")),
                                ## ensure that the sequence number is 2 digits
                                ## and padded with 0 if necessary. Not as
                                ## elegent as zip in python but achieves the
                                ## same result
                                apply(cbind(seq.int(1, length(clusterWhereAmI)), clusterWhereAmI), 1, function(xx) { sprintf("%02d %s", as.integer(xx[1]), xx[2]) } ))
                            
                            cat("*** Some of the melted mgd data frame\n")
                            print (some(melted.mgd))
                            ## stop("Check the melted mgd data frame\n")

                            graph.z.stats(melted.mgd, groups, seedName, statLabel, summaryColumns[2])
                        }
                    } ## end of else if (is.z.stat(statLabel) ) {
                } else { ## end of if (clusterCount > 0 ) {
                    cat ("*** No clusters found\n")
                }
            } else {
                cat("No Clusters,\n\n", file=publicationTableFilename, append=TRUE)
                cat("*** No such file", roistats.filename, "\n")
            } ## end of if(file.exists(roistats.filename)) {
        } ## end of for (seed in seeds) {
    } ## for (statLabel in statsLabels) {
} ## end of generateGraphs definition


graph.f.stats <-function(inMeltedRoistats, inGroups, inSeed, inStatLabel) {

    imageDirectory=file.path(group.results.dir, inSeed, analysis, inStatLabel, inGroups)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory, recursive=TRUE)
    }

    for ( level in levels(inMeltedRoistats$cluster) ) {

        ss=subset(inMeltedRoistats, cluster==level)
        ## print (ss)

        ## the first line is used to change the level labels to match
        ## that used in the amygdala rsfc paper (MDD and HCL rather
        ## than MDD and NCL)
        levels(ss$Group)=c("MDD", "HCL")
        ## this line is used to ensure that HCL will be listed before
        ## MDD in the graphs
        ss$Group = relevel(ss$Group, ref="HCL")
        
        imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%s.%s.%s.%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, inGroups, inSeed, inStatLabel))
        cat(paste("*** Creating", imageFilename, "\n"))

        if (inStatLabel == "Group.F") {
            roistats.summary=summarySE(ss, measurevar="value", groupvars=c("Group", "cluster"))
            x.axis="Group"
            y.axis="value"
            shape="Group"
            color="Group"
            xlabel="Group"
            group=1
            plot.breaks=levels(inMeltedRoistats$Group)
            plot.labels=levels(inMeltedRoistats$Group)
        } else if (inStatLabel == "timepoint.F" ) {
            roistats.summary=summarySE(ss, measurevar="value", groupvars=c("timepoint", "cluster"))                              
            x.axis="timepoint"
            y.axis="value"
            shape="timepoint"
            color="timepoint"
            xlabel="Timepoint"
            group=1
        } else if (inStatLabel == "GroupXtimepoint.F" ) {
            roistats.summary=summarySE(ss, measurevar="value", groupvars=c("Group", "timepoint"))
            ## x.axis="Group"
            ## y.axis="value"
            ## shape="timepoint"
            ## color="timepoint"
            ## group="timepoint"
            ## xlabel="Group"

            x.axis="timepoint"
            y.axis="value"
            shape="Group"
            color="Group"
            group="Group"
            xlabel="Timepoint"

            ss.diff=cbind(ss[ss$timepoint=="Baseline", c("Subj", "Group")],
                          "timepoint"="Difference",
                          "cluster"=ss[ss$timepoint=="Baseline", "cluster"],
                          "value"=ss[ss$timepoint=="Follow-up", "value"] - ss[ss$timepoint=="Baseline", "value"])
            rownames(ss.diff)=NULL
            ## print(ss.diff)
            cat("*** Summary of Follow-up - Baseline values\n")
            print(summarySE(ss.diff, measurevar="value", groupvars=c("Group", "timepoint")))
        } else {
            stop(sprintf("Can't summarize the %s label. Stopping.\n", inStatLabel))
        }
        cat("*** Summary of data used to create scatter plot\n")
        print(roistats.summary)
        ## stop()
        my.dodge=position_dodge(.2)
        ## print(sprintf("*** Got to graph.f.stats %s\n", level))
        ## this works for time point or group
        interactionGraph = ggplot(data=roistats.summary, aes_string(x=x.axis, y=y.axis, color=color, shape=shape, group=group) )
        interactionGraph = interactionGraph + geom_point(size=2)
        interactionGraph = interactionGraph + geom_point(data=ss, size=2)         
        interactionGraph = interactionGraph + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2)
        interactionGraph = interactionGraph + scale_color_manual(values = c("#94bd51", "#f18440")) 
        ## interactionGraph = interactionGraph + scale_color_brewer(name="Timepoint:", palette="Set1")
        ## interactionGraph = interactionGraph + labs(title = substituteShortLabels(level), x=xlabel, y="RSFC (Z-score)")
        interactionGraph = interactionGraph + labs(title = substituteShortLabels(level), x=xlabel, y="RSFC (Z-score)")         
        interactionGraph = interactionGraph + my.theme
        
        if (inStatLabel == "GroupXtimepoint.F") {
            ## interactionGraph = interactionGraph + geom_line(position=my.dodge)
            interactionGraph = interactionGraph + geom_line()            
            interactionGraph = interactionGraph + theme(legend.position="right")
        }

        ## if (setequal(ss$timepoint, c("Baseline", "3 Months", "6 Months"))) {
        ##     interactionGraph = interactionGraph + geom_line()            
        ##     interactionGraph = interactionGraph + theme(legend.position="bottom")
        ## }

        
        
        if (FALSE && setequal(ss$Group, c("MDD", "HCL"))) {
            baseline.test=t.test(ss[ss$Group=="MDD" & ss$timepoint=="Baseline", "value"],
                                 ss[ss$Group=="HCL" & ss$timepoint=="Baseline", "value"])
            print(baseline.test)
            print(tes(baseline.test$statistic,
                      n.1=length(ss[ss$Group=="MDD" & ss$timepoint=="Baseline", "value"]),
                      n.2=length(ss[ss$Group=="HCL" & ss$timepoint=="Baseline", "value"])))
            
            follow.up.test=t.test(ss[ss$Group=="MDD" & ss$timepoint=="Follow-up", "value"],
                                  ss[ss$Group=="HCL" & ss$timepoint=="Follow-up", "value"])
            print(follow.up.test)
            print(tes(follow.up.test$statistic,
                      n.1=length(ss[ss$Group=="MDD" & ss$timepoint=="Baseline", "value"]),
                      n.2=length(ss[ss$Group=="HCL" & ss$timepoint=="Baseline", "value"])))

            mdd.change.ttest=t.test(ss[ss$Group=="MDD" & ss$timepoint=="Follow-up", "value"],
                                    ss[ss$Group=="MDD" & ss$timepoint=="Baseline",  "value"], paired=TRUE)
            print(mdd.change.ttest)
            print(tes(mdd.change.ttest$statistic,
                      n.1=length(ss[ss$Group=="MDD" & ss$timepoint=="Follow-up", "value"]),
                      n.2=length(ss[ss$Group=="MDD" & ss$timepoint=="Baseline",  "value"])))


            hcl.change.ttest=t.test(ss[ss$Group=="HCL" & ss$timepoint=="Follow-up", "value"],
                                    ss[ss$Group=="HCL" & ss$timepoint=="Baseline",  "value"], paired=TRUE)
            print(hcl.change.ttest)
            print(tes(hcl.change.ttest$statistic,
                      n.1=length(ss[ss$Group=="HCL" & ss$timepoint=="Follow-up", "value"]),
                      n.2=length(ss[ss$Group=="HCL" & ss$timepoint=="Baseline",  "value"])))}
        ## print(interactionGraph)
        ##stop()
        ## ggsave(imageFilename, interactionGraph, width=1.45, height=1.6, units="in")
        ggsave(imageFilename, interactionGraph, width=4, height=3, units="in")        
        ## stop()
        ## ggsave(imageFilename, interactionGraph, units="in") 
    } ## end of for ( level in levels(roistats.summary$cluster) )
} ## end of graphFStats


is.main.effect.contrast.z <- function (in.stat.label) {
    main.effects.pattern="[A-Za-z0-9]+-[A-Za-z0-9]+\\.Z"
    return(grepl(main.effects.pattern, in.stat.label, fixed=FALSE))
}


is.interaction.contrast.z <- function (in.stat.label) {
    interaction.pattern="([A-Za-z0-9]+)\\.[A-Za-z0-9]+-\\1\\.[A-Za-z0-9]+\\.Z"
    return(grepl(interaction.pattern, in.stat.label, fixed=FALSE))
}

get.z.contrast.terms <- function (in.stat.label) {
    if (is.main.effect.contrast.z(in.stat.label) ) {
        if (grepl ("MDD", in.stat.label, fixed=TRUE)) {
            return("Group")
        } else if (grepl ("F", in.stat.label, fixed=TRUE)) {
            return("Gender")
        }
    } else {
        interaction.pattern="([A-Za-z0-9]+)\\.([A-Za-z0-9]+)-\\1\\.([A-Za-z0-9]+)\\.Z"
        return(regmatches(in.stat.label, regexec(interaction.pattern, in.stat.label))[[1]][-1])
    }
}

convert.contrast.terms.to.factor.names <- function (in.terms) {

    retVec=vector(mode="character", length=2)

    if ( in.terms[1] %in% c("MDD", "NCL"))
        retVec[1] = "Group"
    if ( all(in.terms[-1] %in% c("F", "M")))
        retVec[2] = "Gender"

    return (retVec)
}
        
graph.z.stats <-function(inMeltedRoistats, inGroups, inSeed, inStatLabel, in.effect) {

    imageDirectory=file.path(group.results.dir, inSeed, analysis, inStatLabel)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory, recursive=TRUE)
    }

    for ( level in levels(inMeltedRoistats$cluster) ) {
        
        ss=subset(inMeltedRoistats, cluster==level)
        
        imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%s.%s.%s.%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, inGroups, inSeed, inStatLabel))
        cat(paste("*** Creating", imageFilename, "\n"))
    
        if (in.effect == "Group") {
            roistats.summary=summarySE(ss, measurevar="value", groupvars=c("Group", "cluster"))
            x.axis="Group"
            y.axis="value"
            shape="Group"
            color="Group"
            xlabel="Group"
            group=1
            plot.breaks=levels(inMeltedRoistats$Group)
            plot.labels=levels(inMeltedRoistats$Group)
        } else if (in.effect == "Gender" ) {
            roistats.summary=summarySE(ss, measurevar="value", groupvars=c("Gender", "cluster"))                              
            x.axis="Gender"
            y.axis="value"
            shape="Gender"
            color="Gender"
            xlabel="Gender"
            group=1
        } else {
            stop(sprintf("Can't summarize the %s label. Stopping.\n", in.effect))
        }

        my.dodge=position_dodge(.2)

        ## this works for time point or group
        graph=ggplot(data=roistats.summary, aes_string(x=x.axis, y=y.axis, color=color, shape=shape, group=group) ) 
        graph=graph + geom_point(position=my.dodge) 
        ## graph=graph + geom_jitter(data=ss)
        graph=graph + geom_point(data=ss)         
        graph=graph + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=my.dodge) 
        graph=graph + scale_shape_discrete(name=paste(shape, ":", sep="")) 
        graph=graph + scale_color_brewer(name=paste(shape, ":", sep=""), palette="Set1")
        graph=graph + labs(title = substituteShortLabels(level), x=xlabel, y="RSFC (Z-score)") 
        graph=graph + my.theme
        
        if (inStatLabel == "GroupXGender.F") {
            graph=graph + geom_line(position=my.dodge)
        }
        
        print(graph)
        stop()
        ## ggsave(imageFilename, interactionGraph, width=3, height=3, units="in")
        ggsave(imageFilename, graph, units="in") 
    } ## end of for ( level in levels(roistats.summary$cluster) )
} ## end of graph.main.interaction.z.stats

####################################################################################################
### End of functions
####################################################################################################


if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=normalizePath(file.path(data.dir, "admin"))
config.data.dir=normalizePath(file.path(data.dir, "config"))
seeds.data.dir=normalizePath(file.path(data.dir, "seeds"))

## used for the group-by-timepoint analysis
group.data.dir=normalizePath(file.path(data.dir, "Group.data.followup", "group.by.timepoint.lmes"))
group.results.dir=normalizePath(file.path(data.dir, "Group.results.followup", "group.by.timepoint.lmes"))

## used for the within MDD only timepoint analysis
## group.data.dir=normalizePath(file.path(data.dir, "Group.data.followup"))
## group.results.dir=normalizePath(file.path(data.dir, "Group.results.followup", "A.C.D"))

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
## demographicsFilename=file.path(admin.data.dir, "merged.demographics.and.neuropsych.csv")
demographicsFilename=file.path(admin.data.dir, "merged.demographics.and.neuropsych.2016.11.07.csv")
demographics=readCsvFile(demographicsFilename)
demographics=rename(demographics, c("Grp"="Group"))

clust.header = c("Volume", "CM RL", "CM AP", "CM IS", "minRL",
    "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
    "MI RL", "MI AP", "MI IS")

task="restingstate"
usedFwhm="7.98x8.05x7.32"


##################################################

do.graph.creation=TRUE
do.lmes=FALSE
do.delta.correlations=FALSE

##################################################


##################################################
## analysis="all.timepoints"
## analysis="any.timepoint"
## analysis.methods=c("all.timepoints", "any.timepoint")
## analysis.methods=c("all.timepoints")
analysis.methods=c("timepoint.A.and.C")
## analysis.methods=c("timepoint.A.and.C.and.D")

##################################################
## seedFiles=
##     sapply(c(
##         "juelich_whole_amygdala_seeds.txt",
##         "short_ACC_seed_list.txt",
##         "hippocampus_ventral_striatum_seeds.txt",
##     	"followup-dlpfc-ins-IP-MPFC-seeds.txt",
##         "Fox-Goldapple-seeds.txt",
##         "miller-dmn.txt",
##         "jacobs-seeds.txt",
##         "goldapple-ofc-seeds.txt",
##         "gabbay-striatum-seeds.txt",
##         "tremblay-seeds.txt"
##     ),
##     function(xx) {
##         file.path(config.data.dir, xx)
##     })
## seeds=unlist(sapply(seedFiles, readSeedsFile))

seedFiles=file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt")
seeds=readSeedsFile(file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt"))

numberOfSeeds=length(seeds)

my.base.size=18
my.theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        legend.position="none",
        ## legend.position="bottom",        
        ## panel.grid.major = element_blank(),
        ## panel.grid.minor = element_blank(),

        ##remove the panel border
        ## panel.border = element_blank(),

        ## add back the axis lines
        axis.line=element_line(colour = "grey50"),

        axis.text.x=element_text(size=my.base.size*0.9),
        axis.text.y=element_text(size=my.base.size*0.9),       
        
        axis.title.x=element_blank(),
        ## axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size, vjust=1))

## output.filename=file.path(group.results.dir, paste("roi-graphing-results-output", format(Sys.time(), "%Y%m%d-%H%M%Z"), "txt", sep="."))
## cat("*** Output table is in ", output.filename, "\n")
## ff=file(output.filename, open="w", encoding="utf-8")
## sink(ff, append=FALSE)

####################################################################################################
####################################################################################################
####################################################################################################

if (do.graph.creation) {
    
    groups.to.stats.labels=
        list(
            "MDD.and.NCL"=
                c(
                    ## Main and interaction effect F statistics follow
                    ## "Group.F",
                    ## "timepoint.F",
                    "GroupXtimepoint.F"#,
                    
                    ## contrast related Z scores follow
                    ##"MDD-NCL.Z",
                    ##"M-F.Z",
                    ##"MDD.F-MDD.M.Z",
                    ##"NCL.F-NCL.M.Z"
                )
            ## ,
            ## "MDD"=c(
            ##     "timepoint.F"
            ## )
        )
    
    for (analysis in analysis.methods) {
        for (groups in names(groups.to.stats.labels)) {
            statsLabels=groups.to.stats.labels[[groups]]
            
            cat("****************************************************************************************************\n")
            cat(sprintf("*** Processing stats labels for the %s group(s) for the %s combination\n", groups, analysis))
            print(statsLabels)
            
            generateGraphs()    
        }
    }
}
####################################################################################################
####################################################################################################
####################################################################################################

if (do.delta.correlations) {

    groups.to.stats.labels=
        list(
            "MDD.and.NCL"=
                c(
                    ## Main and interaction effect F statistics follow
                    "GroupXtimepoint.F"#,
                    
                    ## contrast related Z scores follow
                    ##"MDD-NCL.Z",
                    ##"M-F.Z",
                    ##"MDD.F-MDD.M.Z",
                    ##"NCL.F-NCL.M.Z"
                )
           ## ,
           ##  "MDD"=c(
           ##      "timepoint.F"
           ##  )
        )

    regressionVariables=list(
        ## list(variable="MASC.tscore",        name="MASC"),
        ## list(variable="BDI.II.Total",       name="BDI Inventory II"),
        ## list(variable="CDI.Total",          name="CDI"),
        ## list(variable="CDRS.t.score",       name="CDRS-R")
        ## list(variable="CGAS",               name="CGAS"),
        ## list(variable="RADS.Total.tscore",  name="RADS")

        ## list(variable="RSQ.Total",             name="Ruminative Responses Styles Questionnaire"),
        list(variable="CDRS.t.score",         name="CDRS-R"),
        list(variable="RADS.DM.tscore",        name="RADS Dysphoric Mood"),        
        list(variable="RADS.AN.tscore",        name="RADS Anhedonia/Negative Affect"),
        list(variable="RADS.NS.tscore",        name="RADS Negative Self-evaluation"),
        list(variable="RADS.SC.tscore",        name="RADS Somatic Complaints"),
        list(variable="RADS.Total.tscore",     name="RADS Total"),
        list(variable="C_Irritability.Total",  name="Caprara Irritability Scale"),
        list(variable="C_EmoSuscept.Total",    name="Caprara Emotional Susceptibility Scale")
    )

    results.stack=stack()
    for (analysis in analysis.methods) {
        for (groups in names(groups.to.stats.labels)) {
            statsLabels=groups.to.stats.labels[[groups]]
            
            cat("****************************************************************************************************\n")
            cat(sprintf("*** Processing stats labels for the %s group(s) for the %s combination\n", groups, analysis))
            print(statsLabels)

            results.stack=run.delta.correlations(statsLabels, results.stack)

        }
    }

    header="groups,analysis,seedName,statLabel,group,cluster,ROI.No,RL,AP,IS,RegressionName,N.with.incomplete.cases,N,Stat,DoF,pValue,r,Significance\n"
    results.table.filename=file.path(group.results.dir, paste("correlations", format(Sys.time(), "%Y%m%d-%H%M%Z"), "csv", sep="."))
    cat("*** Correlation table is in ", results.table.filename, "\n")

    cat (header, file=results.table.filename, append=FALSE)
    l=results.stack$value()
    for (i in 1:length(l)) {
        cat (l[[i]], "\n", file=results.table.filename, append=TRUE)
    }
}

## sink()
