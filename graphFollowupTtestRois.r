rm(list=ls())
graphics.off()

## car for the some function
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
    require(plyr)

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
                                 inStatColumnName="Default Stat Name", inContrastColumnName="Default Contrast Name", inCom=TRUE) {
  hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
  ##print(hemisphere)
  if ( inCom ) {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
      colnames(locations) = c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS")
  } else {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
      colnames(locations) = c("Structure", "Hemisphere", "Volume", "MI RL", "MI AP", "MI IS")      
  }

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
  ddply.agg=ddply(inRoistats, .(Group, timepoint),
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

  ##print(which (! grepl("Mean", colnames(ddply.agg))))
  ##ddply.agg[, which (! grepl("Mean", colnames(ddply.agg)))]
  ## now make the column names for the timepoint x group 
  cnames=apply(ddply.agg[, which (! grepl("Mean", colnames(ddply.agg)))], 1,
      function(xx) {
          ## xx[1] is group, xx[2] is timepoint
          return(sprintf("%s (%s)", xx[1], xx[2]))
      })
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
      colnames(pubTable)=c(colnames(pubTable)[-length(colnames(pubTable))], inStatColumnName)
      ## print(pubTable)
  } else {
      pubTable=cbind(locations, mns)
  }

  ## Formula for Cohen's d_z from http://journal.frontiersin.org/article/10.3389/fpsyg.2013.00863/full
  dz=sapply(which(grepl("Mean", colnames(inRoistats))), function (xx) {
      difference=inRoistats[inRoistats$timepoint=="C", xx] - inRoistats[inRoistats$timepoint=="A", xx]
      ## cat("*** difference\n")
      ## print(difference)
      mean.difference=mean(difference)
      ## cat("*** mean.difference\n")      
      ## print(mean.difference)
      dz=mean.difference / sqrt(sum((((difference - mean.difference)^2)/(dim(inRoistats)[1]-1))))
      ## cat("*** dz\n")            
      ## print(dz)
  })

  ## cat("*** inRoistats\n")
  ## print(inRoistats)
  ## cat("*** dz\n")  
  ## print(dz)

  pubTable=cbind(pubTable, round(dz, 4))
  colnames(pubTable)=c(colnames(pubTable)[-length(colnames(pubTable))], "Cohen's d_z")
  ## stop()
  
  if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageContrastValue) ) {
      ## cat("Adding average contrast values\n")      
      pubTable=cbind(pubTable, round(t(inRoistats.averageContrastValue), 2))
      colnames(pubTable)=c(colnames(pubTable)[-length(colnames(pubTable))], inContrastColumnName)
      ## print(pubTable)      
  }

  ## pubTable=cbind(pubTable, mns)
  ## print(pubTable)
  ## stop()
  ## if (! is.null(inRoistats.averageStatValue) ) {
  ##     colnames(pubTable)=
  ##         c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", inStatColumnName, "Cohen's d_z", colnames(mns))
  ## } else {
  ##     colnames(pubTable)=c("Structure", colnames(pubTable)[-1], "Cohen's d_z")
  ## }

  ## if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageContrastValue) ) {
  ##     colnames(pubTable)=
  ##         c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", inStatColumnName, inContrastColumnName, colnames(mns))
  ## }
  
  rownames(pubTable)=NULL
  ## print(pubTable)
  ## stop()
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

readSubjectOrderTable <- function (inFilename) {
    cat("*** Reading", inFilename, "\n")
    subjectOrder=read.csv(inFilename, header=T)

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

readDataTable <- function (inFilename) {
    cat("*** Reading", inFilename, "\n")
    dataTable=read.table(inFilename, header=T, allowEscapes=TRUE)

    return(dataTable)
}

fixDataTable <- function (inDataTable) {
    inDataTable$Subj=gsub("S\\.([0-9]{3})", "\\1", as.character(inDataTable$Subj), fixed=FALSE)
    if (any(grepl("300", inDataTable$Subj, fixed=TRUE))) {
        inDataTable$Subj=gsub("300", "169/300", as.character(inDataTable$Subj), fixed=TRUE)
        inDataTable$Subj=as.factor(inDataTable$Subj)
    }
    return(inDataTable)
}

generateGraphs <- function (seed.list) {

    for (seed in seeds) {
        seedName=getSeedName(seed)

        publicationTableFilename=file.path(group.results.dir, paste("publicationTable.ttest", usedFwhm, groups, seedName, "all.timepoints", "csv", sep="."))
        if (file.exists(publicationTableFilename)) {
            file.remove(publicationTableFilename)
        }
        ## cat("*** Writing publication table to", publicationTableFilename, "\n")
        
        cat("####################################################################################################\n")
        cat(sprintf("*** Graphing ROIs from the %s seed for the %s groups\n", seedName, groups))
        
        infix=sprintf("ttest.fwhm%s.%s.%s.timepoint.A.and.C", usedFwhm, groups, seedName)
        ## roiStats.ttest.fwhm7.98x8.05x7.32.MDD.R_whole_amygdala.3mm.timepoint.A.and.C.txt
        roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))

        roistats.averageTvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageTValue.txt", infix))
        roistats.averageContrastValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageContrastValue.txt", infix))
            
        if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the contrast in each ROI,
            ## you do not what to graph this
            
            roistats=readStatsTable(roistats.filename)
            roistats$Sub.brick=NULL

            roistats.averageTvalue=readStatsTable(roistats.averageTvalue.filename)
            roistats.averageContrastValue=readStatsTable(roistats.averageContrastValue.filename)
            roistats.averageTvalue$Sub.brick=NULL
            roistats.averageContrastValue$Sub.brick=NULL                
            
            clusterCount=length(grep("Mean", colnames(roistats)))
            if (clusterCount > 0 ) {
                cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
                
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
                
                clustersFilename=file.path(group.results.dir, sprintf("clust.%s.txt", infix))
                clusters=readClustersTable(clustersFilename)
                
                ## this table contains the locations, as text, of the clusters and is the output of a perl script
                clusterLocationsFilename=file.path(group.results.dir, sprintf("clusterLocations.%s.csv", infix))
                clusterWhereAmI=readClusterLocationsTable(clusterLocationsFilename)
                
                ## this file stores the order of the subjects in each of the following BRIK files
                dataTableFilename=file.path(group.data.dir, paste("dataTable.ttest.rsfc", groups, seedName, "all.timepoints", "txt", sep="."))
                dataTable=fixDataTable(readDataTable(dataTableFilename))
                ## print(dataTable)

                mgd=cbind(dataTable, roistats, demographics[match(dataTable$Subj, demographics$ID), c("Gender")])
                colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], "Gender")
                ## print(mgd)
                ## stop()

                rownames(mgd)=NULL
                ## ensure that subject is s factor
                mgd$Subj=as.factor(mgd$Subj)
                mgd=droplevels(mgd)
                
                ## print(clusterWhereAmI)
                ## print(clusters)
                ## print(roistats)
                ## print(roistats.averageTvalue)
                ## print(roistats.averageContrastValue)                
                ## print(mgd)
                cat("*** Some of the mgd data frame\n")
                print(some(mgd))
                print(addmargins(table(mgd$timepoint)))

                ## stop("Check the mgd data frame\n")

                publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd, inCom=TRUE)
                print(publicationTable)
                ## stop()
                savePublicationTable(publicationTable, publicationTableFilename, TRUE)

                ## print(publicationTable)
                ## stop("Check the publication data frame\n")
                melted.mgd=melt(mgd,  id.vars=c("Subj", "Group", "timepoint"),
                    measure.vars=paste("Mean_", seq(1, clusterCount), sep=""),
                    variable_name="cluster")
                
                melted.mgd$cluster=factor(melted.mgd$cluster,
                    levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                    labels=paste(sprintf("%02d", seq(1, clusterCount)), clusterWhereAmI))
                
                print (head((melted.mgd)))
                ## stop("Check the melted mgd data frame\n")
                
                graphTStats(melted.mgd, groups, seedName)
                ## stop("Do the graphs look ok?\n")
                
            } ## end of if (clusterCount > 0 ) {
        } else {
            cat("No Clusters,\n\n", file=publicationTableFilename, append=TRUE)
        } ## end of if(file.exists(roistats.filename)) {
    } ## end of for (seed in seeds) {
    
} ## end of generateGraphs definition


graphTStats <-function(inMeltedRoistats, inGroups, inSeed) {

    imageDirectory=file.path(group.results.dir, inSeed, analysis, "ttest", inGroups)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory, recursive=TRUE)
    }

    for ( level in levels(inMeltedRoistats$cluster) ) {

        ss=subset(inMeltedRoistats, cluster==level)

        imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%s.%s.%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, inGroups, inSeed))
        cat(paste("*** Creating", imageFilename, "\n"))
        
        roistats.summary=summarySE(ss, measurevar="value", groupvars=c("timepoint", "cluster"))
        print(roistats.summary)
        x.axis="timepoint"
        y.axis="value"
        xlabel="Timepoint"

        my.dodge=position_dodge(.2)

        ## this works for time point or group
        graph=ggplot(data=roistats.summary, aes_string(x=x.axis, y=y.axis) )
        ##geom_point(position=my.dodge, size=0.5)
        graph=graph + geom_bar(stat="identity")
        ## geom_jitter(data=ss)
        graph=graph + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.5, size=1, color="black", position=my.dodge)
        graph=graph + labs(title = substituteShortLabels(level), x=xlabel, y="RSFC (Z-score)")
        graph=graph + my.theme 
        ## print(graph)
        ## stop()
        ## ggsave(imageFilename, graph, width=4, height=3.5, units="in")
        ggsave(imageFilename, graph)
    } ## end of for ( level in levels(roistats.summary$cluster) )
}

run.correlations <- function (in.results.stack, in.method="pearson") {

    for (seed in seeds) {
        seedName=getSeedName(seed)
        
        cat("####################################################################################################\n")
        cat(sprintf("*** Graphing ROIs from the %s seed for the %s groups\n", seedName, groups))
        
        infix=sprintf("ttest.fwhm%s.%s.%s.timepoint.A.and.C", usedFwhm, groups, seedName)
        ## roiStats.ttest.fwhm7.98x8.05x7.32.MDD.R_whole_amygdala.3mm.timepoint.A.and.C.txt
        roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))

        roistats.averageTvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageTValue.txt", infix))
        roistats.averageContrastValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageContrastValue.txt", infix))
            
        if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the contrast in each ROI,
            ## you do not what to graph this
            
            roistats=readStatsTable(roistats.filename)
            roistats$Sub.brick=NULL

            roistats.averageTvalue=readStatsTable(roistats.averageTvalue.filename)
            roistats.averageContrastValue=readStatsTable(roistats.averageContrastValue.filename)
            roistats.averageTvalue$Sub.brick=NULL
            roistats.averageContrastValue$Sub.brick=NULL                
            
            clusterCount=length(grep("Mean", colnames(roistats)))
            if (clusterCount > 0 ) {
                cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
                
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
                
                clustersFilename=file.path(group.results.dir, sprintf("clust.%s.txt", infix))
                clusters=readClustersTable(clustersFilename)
                
                ## this table contains the locations, as text, of the clusters and is the output of a perl script
                clusterLocationsFilename=file.path(group.results.dir, sprintf("clusterLocations.%s.csv", infix))
                clusterWhereAmI=readClusterLocationsTable(clusterLocationsFilename)
                
                ## this file stores the order of the subjects in each of the following BRIK files
                dataTableFilename=file.path(group.data.dir, paste("dataTable.ttest.rsfc", groups, seedName, "all.timepoints", "txt", sep="."))
                dataTable=fixDataTable(readDataTable(dataTableFilename))
                ## print(dataTable)

                ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
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
                                                       interaction(dataTable$Subj,  dataTable$timepoint),
                                                       interaction(demographics$ID, demographics$timepoint)
                                                   ),
                                                   c("Gender", regression.variable)])
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
                    ## print(mgd)
                    ## stop("Check the mgd data frame\n")
                    
                    for (mean.col.itr in seq_len(clusterCount)) {
                        
                        mean.column = paste("Mean", mean.col.itr, sep="_")
                        
                        mgd=mgd[order(mgd$Subj, mgd$timepoint), ]
                        ## print(mgd)
                        mgd.diff = cbind (mgd[mgd$timepoint=="A", c("Subj", "Group")],
                                          mgd[mgd$timepoint=="C", mean.column]         - mgd[ mgd$timepoint=="A", mean.column],
                                          mgd[mgd$timepoint=="C", regression.variable] - mgd[ mgd$timepoint=="A", regression.variable])
                        colnames(mgd.diff)=c("Subj", "Group", paste(c(mean.column, regression.variable), "diff", sep="."))
                        rownames(mgd.diff)=NULL
                        ## print(mgd.diff)
                        ## stop()
                        mgd.diff.dim.before.complete.cases=dim(mgd.diff)
                        mgd.diff=mgd.diff[complete.cases(mgd.diff), ]
                        print(mgd.diff)
                        ## stop()
                        print(addmargins(table(mgd.diff[ , c("Group")])))
                        
                        ct.all=cor.test(mgd.diff[, paste(mean.column,         "diff", sep=".")],
                                        mgd.diff[, paste(regression.variable, "diff", sep=".")],
                                        method=in.method)
                        es=res(ct.all$estimate, n=dim(mgd.diff)[1], verbose=FALSE)
                        ## print(ct.all)
                        
                        csvLine=
                            sprintf("%s,%s,%s,%s,%s,%s,%s,%d,%d,t(%0.2f)=%0.3f,%0.4f,%0.4f,%0.4f,%s",
                                    groups, analysis, seedName, "all",
                                    substituteShortLabels(clusterWhereAmI[mean.col.itr]),
                                    paste(c(mean.col.itr, clusters[mean.col.itr, c("CM RL", "CM AP", "CM IS")]), collapse=","),
                                    regression.name,
                                    mgd.diff.dim.before.complete.cases[1],
                                    dim(mgd.diff)[1],
                                    ifelse(is.null(ct.all$parameter), NA, ct.all$parameter), ct.all$statistic,
                                    es$d,
                                    ct.all$p.value,
                                    ct.all$estimate,
                                    make.significance.indications(ct.all$p.value))
                        
                        ## print(csvLine)
                        ## print(ct.all$p.value)
                        ## print(make.significance.indications(ct.all$p.value))
                        push(in.results.stack, csvLine)
                        
                        if (ct.all$p.value < 0.1) {
                            graph.correlation(mgd.diff, seedName, "ttest", "all", mean.col.itr, clusterWhereAmI[mean.col.itr],
                                              paste(mean.column,         "diff", sep="."), ## x axis 
                                              paste(regression.variable, "diff", sep="."), ## y axis
                                              regression.name,
                                              expression(paste(Delta, "RSFC (Z-score)", sep=" ")),
                                              bquote(paste(Delta, .(regression.name), sep=" ")))
                        } ## end of if (ct.all$p.value < 0.1)
                        ## stop()
                    } ## end of for (mean.col.itr in seq_len(clusterCount))
                } ## end of for (cor.var.itr in seq_along(regressionVariables))
                ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##                
            } ## end of if (clusterCount > 0 ) {
        } else {
            cat("*** No Clusters\n")
        } ## end of if(file.exists(roistats.filename)) {
    } ## end of for (seed in seeds) {
    return(in.results.stack)
} ## end of run.correlations definition


## rv= regression variable
graph.correlation <- function (in.mgd.diff, in.seed.name, in.stat.label, in.group, in.roi.no, in.cluster.name,
                               in.x.axis, in.y.axis, in.rv.name, in.x.label, in.y.label) {

    my.dodge=position_dodge(.2)
    graph=ggplot(in.mgd.diff, aes_string(x=in.x.axis, y=in.y.axis))
    graph = graph + geom_point(position=my.dodge)
    graph = graph + geom_smooth(method="lm", se=FALSE, color="black")
    graph = graph + labs(title = substituteShortLabels(in.cluster.name), x=in.x.label, y=in.y.label) 
    graph = graph + my.theme
    
    image.filename=file.path(group.results.dir, paste(groups, analysis, in.seed.name, in.stat.label, in.group,
                                                      in.y.axis, in.roi.no,
                                                      gsub("[[:space:]]+", ".", in.cluster.name, fixed=FALSE), "pdf", sep="."))
    ## print(graph)

    cat("*** Saving graph to" , image.filename, "\n")
    ggsave(image.filename, graph)
}


run.clinical.ttests <- function(in.results.stack) {


    seedName="L_whole_amygdala.3mm"
    ## this file stores the order of the subjects in each of the following BRIK files
    dataTableFilename=file.path(group.data.dir, paste("dataTable.ttest.rsfc", groups, seedName, "all.timepoints", "txt", sep="."))
    dataTable=fixDataTable(readDataTable(dataTableFilename))
    ## print(dataTable)
    
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
    for (cor.var.itr in seq_along(regressionVariables))  {
        
        regression.variable=regressionVariables[[cor.var.itr]]$variable
        regression.name=regressionVariables[[cor.var.itr]]$name                        
        cat("*** Running t tests for on", regression.variable, "\n")
        
        ## The use of the match and interaction comes
        ## from:
        ## http://stackoverflow.com/questions/6880450/match-two-columns-with-two-other-columns
        ## this could probably be done more
        ## efficiently using the data tables instead
        ## of data frames
        
        mgd=cbind(dataTable, demographics[ match(
                                 interaction(dataTable$Subj,  dataTable$timepoint),
                                 interaction(demographics$ID, demographics$timepoint)
                             ),
                             c("Gender", regression.variable)])
        ## stop()
        if (any(mgd$Subj=="378") ) {
            mgd[mgd$Subj=="378", "Gender"]="F"
        }

        ## print(mgd)
        mgd=mgd[complete.cases(mgd), ]

        ## cat("*** Before filtering to ensure that", regression.variable, "is available at both timepoints\n")
        rv.by.timepoint=table(mgd[, c("Subj", "timepoint")])
        ## print(addmargins(rv.by.timepoint))
        rv.by.timepoint.df=as.data.frame.matrix(rv.by.timepoint)
        rv.by.timepoint.df$Subj = rownames(rv.by.timepoint.df)
        rownames(rv.by.timepoint.df)=NULL
        rv.by.timepoint.df$at.both.timepoints=(rv.by.timepoint.df$A + rv.by.timepoint.df$C) == 2
        rv.by.timepoint.df=rv.by.timepoint.df[, c("Subj", "A", "C", "at.both.timepoints")]
        rv.by.timepoint.df=subset(rv.by.timepoint.df, at.both.timepoints == TRUE)
        rownames(rv.by.timepoint.df)=NULL
        
        ## cat("*** AFTER *** Filtering to include only subjects with", regression.variable, "data at both time points\n")
        mgd=subset(mgd, mgd$Subj %in% rv.by.timepoint.df$Subj)
        mgd = droplevels(mgd [ order(mgd$Group, mgd$Subj, mgd$timepoint), ])
        
        mgd$Subj = as.factor(paste("S", mgd$Subj, sep="."))
        rownames(mgd)=NULL
        
        n.table=table(mgd[, c("Group", "timepoint")])
        n.table=addmargins(n.table)
        ## print(n.table)

        ## print(mgd)
        t.test.formula=as.formula(paste(regression.variable, "~", "timepoint"))
        rv.ttest=t.test(t.test.formula, data=mgd, paired=TRUE)
        ## print(rv.ttest)

        rv.summary=summarySE(mgd, measurevar=regression.variable, groupvars=c("timepoint"))
        ## print(rv.summary)
        csvLine=sprintf("%s,%d,%d,t(%0.2f)=%0.2f,%0.2f,%s,%02f,%0.2f",
                        regression.name, n.table["MDD", "A"], n.table["MDD", "C"], rv.ttest$parameter, rv.ttest$statistic, rv.ttest$p.value,
                        make.significance.indications(rv.ttest$p.value),
                        rv.summary[rv.summary$timepoint=="A", regression.variable],
                        rv.summary[rv.summary$timepoint=="C", regression.variable]                        
                        )
        ## cat(csvLine, "\n")
        push(in.results.stack, csvLine)
    }
    return(in.results.stack)
}


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

data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
seeds.data.dir=file.path(data.dir, "seeds")

group.data.dir=file.path(data.dir, "Group.data.followup")
group.results.dir=file.path(data.dir, "Group.results.followup", "new_ttests")

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographicsFilename=file.path(admin.data.dir, "merged.demographics.and.neuropsych.2016.11.07.csv")
demographics=readCsvFile(demographicsFilename)
demographics=rename(demographics, c("Grp"="Group"))

clust.header = c("Volume", "CM RL", "CM AP", "CM IS", "minRL",
    "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
    "MI RL", "MI AP", "MI IS")

##################################################

do.graph.creation=TRUE
do.correlations=FALSE

## this controls the running of ttests for clinical variables defined
## in the regressionVariables list below are run. This involves ttests
## on the values of each of the clinical variables for the subjects
## included in the MD only baseline vs. follow-up voxelwise ttests
do.clinical.ttests=FALSE

task="restingstate"
usedFwhm="7.98x8.05x7.32"
analysis="all.timepoints"
groups="MDD"

##################################################

seedFiles=
    sapply(c("short_ACC_seed_list.txt",
             "followup-dlpfc-ins-IP-MPFC-seeds.txt",
             "Fox-Goldapple-seeds.txt"
             ),
           function(xx) {
               file.path(config.data.dir, xx)
           })
seeds=unlist(sapply(seedFiles, readSeedsFile))

## seeds=readSeedsFile(file.path(config.data.dir, "short_ACC_seed_list.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "miller-dmn.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "jacobs-seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "goldapple-ofc-seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "gabbay-striatum-seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "tremblay-seeds.txt"))

numberOfSeeds=length(seeds)
cat(sprintf("*** Found %02d seeds in the seed file\n", length(seeds)))

my.base.size=14
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
        
        ##axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))

regressionVariables=list(

    list(variable="RADS.DM.tscore",       name="RADS Dysphoric Mood"),
    list(variable="RADS.AN.tscore",       name="RADS Anhedonia/Negative Affect"),
    list(variable="RADS.NS.tscore",       name="RADS Negative Self-evaluation"),
    list(variable="RADS.SC.tscore",       name="RADS Somatic Complaints"),
    list(variable="RADS.Total.tscore",    name="RADS Total"),
    list(variable="C_Irritability.Total", name="Caprara Irritability Scale"),
    list(variable="CGAS",                 name="CGAS"),
    list(variable="MASC.tscore",          name="MASC"),
    list(variable="CDRS.t.score",         name="CDRS-R")
)

output.filename=file.path(group.results.dir, paste("effectsizes.ttests.results.output", format(Sys.time(), "%Y%m%d-%H%M%Z"), "txt", sep="."))
cat("*** Output table is in ", output.filename, "\n")
ff=file(output.filename, open="w", encoding="utf-8")
sink(ff, append=FALSE)

if (do.graph.creation) {
    generateGraphs(seed.list)
}

if (do.clinical.ttests) {
    
    results.stack=stack()
    
    results.stack=run.clinical.ttests(results.stack)
    header="Variable,Baseline N,Follow-up N,Statistic,p value,Significance,Baseline Mean,Follow-up Mean\n"
    results.table.filename=file.path(group.results.dir, paste("clinical.ttest", format(Sys.time(), "%Y%m%d-%H%M%Z"), "csv", sep="."))
    cat("*** T tests table is in ", results.table.filename, "\n")
    
    cat (header, file=results.table.filename, append=FALSE)
    l=results.stack$value()
    for (i in 1:length(l)) {
        cat (l[[i]], "\n", file=results.table.filename, append=TRUE)
    }

}


if (do.correlations) {

    results.stack=stack()

    results.stack=run.correlations(results.stack)
    
    header="groups,analysis,seedName,statLabel,cluster,ROI.No,RL,AP,IS,RegressionName,N.with.incomplete.cases,N,Stat,Cohen's d,pValue,r,Significance\n"
    results.table.filename=file.path(group.results.dir, paste("sj.correlations.ttest", format(Sys.time(), "%Y%m%d-%H%M%Z"), "csv", sep="."))
    cat("*** Correlation table is in ", results.table.filename, "\n")

    cat (header, file=results.table.filename, append=FALSE)
    l=results.stack$value()
    for (i in 1:length(l)) {
        cat (l[[i]], "\n", file=results.table.filename, append=TRUE)
    }
}
sink()
