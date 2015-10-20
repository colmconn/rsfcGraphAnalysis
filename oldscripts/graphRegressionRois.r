rm(list=ls())
graphics.off()

##library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
library(grid)
library(plyr)
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

makePublicationTable <- function(inClusterWhereAmI, inClusters, inRoistats, inRoistats.averageStatValue=NULL, inRoistats.averageCoefficientValue=NULL,
                                 inStatColumnName="Default Stat Name", inCoefficientColumnName="Default Coefficient Name", inCom=TRUE) {
  hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
  ## print(hemisphere)
  if ( inCom ) {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
  } else {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
  }

  ## cat("Locations: Volume and coordinates\n")
  ## print(locations)

  ##cat ("Columns matching Mean: ", grep("Mean", colnames(inRoistats)), "\n")
  ##cat ("Data from the above columns:", inRoistats[, grep("Mean", colnames(inRoistats))], "\n")
  ##print(list(inRoistats$Group))
  agg=aggregate(inRoistats[, grep("Mean", colnames(inRoistats))], list(inRoistats$Group), mean)
  ## cat("agg: mean for each group in each ROI\n")  
  ## print(agg)
  mns=round(t(agg[, -1]), 2)
  ## cat("mns: transposed mean for each group in each ROI\n")    
  colnames(mns)=levels(agg[,1])
  ## print(mns)

  ## cat("inRoistats.averageStatValue\n")
  ## print (inRoistats.averageStatValue)
  if (! is.null(inRoistats.averageStatValue) ) {
      ## cat("Adding average t stats\n")
      pubTable=cbind(locations, round(t(inRoistats.averageStatValue), 2))
      ## print(pubTable)
  } else {
      pubTable=cbind(locations, mns)
  }

  if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageCoefficientValue) ) {
      ## cat("Adding average coefficient values\n")      
      pubTable=cbind(pubTable, round(t(inRoistats.averageCoefficientValue), 2))
      ## print(pubTable)      
  }

  pubTable=cbind(pubTable, mns)
  if (! is.null(inRoistats.averageStatValue) ) {
      colnames(pubTable)=
          c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", inStatColumnName, colnames(mns))
  } else {
      colnames(pubTable)=c("Structure", colnames(pubTable)[-1])
  }

  if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageCoefficientValue) ) {
      colnames(pubTable)=
          c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", inStatColumnName, inCoefficientColumnName, colnames(mns))
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

readSubjectOrderTable <- function (inFilename) {
    cat("*** Reading", inFilename, "\n")
    subjectOrder=read.csv(inFilename, header=T)

    return(subjectOrder)
}

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
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

replaceNa <- function (inData, inColumns, inSubjectColumnName="ID", inGroupColumnName="Grp", inReplaceNaWith=NULL) {

    n=length(inColumns)
    if ( ! is.null (inReplaceNaWith) && length(inReplaceNaWith) != n) {
        stop("*** ERROR: Length of inReplaceNaWith does not match that of inColumns. Cannot continue. Stopping.\n")
    }
    for (ii in seq.int(1, n)) {
        column=inColumns[ii]
        if (any(is.na(inData[, column]))) {
            cat ("****************************************************************************************************\n")
            cat (sprintf("*** The following subjects have NA data for %s\n", column))

            cat (paste (as.vector ( inData[is.na(inData[, column]), inSubjectColumnName]), collapse=" "), "\n")
            cat (paste (as.vector ( inData[is.na(inData[, column]), inGroupColumnName]), collapse=" "), "\n")            
            
            ##cat(paste(as.vector(is.na(inData[inData[, column]) & ! is.na(inData[, column]), inSubjectColumnName]), collapse=" "), "\n")
            ##cat(paste(as.vector(is.na(inData[inData[, column]) & ! is.na(inData[, column]), inGroupColumnName]), collapse=" "), "\n")
            ##cat(paste(as.vector(is.na(inData[inData[, column]) & ! is.na(inData[, column]), column]), collapse=" "), "\n")
            if (! is.null(inReplaceNaWith)) {

                replacement.value=inReplaceNaWith[ii]
                if ( ! replacement.value %in% levels(inData[, column]) ) {
                    cat("*** Warning: The replacement value", replacement.value, "is not among the levels of the", column, "column.\n")
                    cat("*** Warning: Converting", column, "to character and back to a factor to accommodate this.\n")
                    inData[, column] = as.character(inData[, column])
                }

                cat ("*** Now setting these values to NA\n")
                inData[ which(is.na(inData[, column])), column]=replacement.value

                if ( ! inReplaceNaWith %in% levels(inData[, column]) ) {
                    inData[, column] = as.factor(inData[, column])
                }
            }
            cat ("****************************************************************************************************\n")      
        }
    } ## end of for (column in inColumns) {
    
    return(inData)
} ## end of checkIsNa


generateGraphs <- function (group.data.dir, group.results.dir, rvVariable, rvName, publicationTableFilename, seed.list) {


    if (groups == "mddAndCtrl" & grepl("diff", rvVariable)) {
        print("*** You are trying to regress/graph  the control and MDD subjects aagainst a time A to C difference. you cannot do that!\n")
        return()
    }
    
    for (seed in seeds) {
        seedName=getSeedName(seed)
        
        for (polarity in c("negative", "positive")) {

            cat(sprintf("Seed=%s,Variable=%s (%s),Polarity=%s,\n", seedName, gsub("\n", " ", rvName), rvVariable, polarity), file=publicationTableFilename, append=TRUE)
            
            cat("####################################################################################################\n")
            cat(sprintf("*** Graphing ROIs for the %s seed that was regressed against %s (%s) in the %s group\n", seedName, rvName, rvVariable, groups))

            infix=sprintf("regression.fwhm%s.%s.%s.%s.and.%s.%s", usedFwhm, task, groups, seedName, rvVariable, polarity)

            roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))
            roistats.averageTvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageTValue.txt", infix))
            roistats.averageCoefficientValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageCoefficientValue.txt", infix))

            if(file.exists(roistats.filename)) {
                
                ## roistats contains the avergae from the contrast in each ROI,
                ## you do not what to graph this
            
                roistats=readStatsTable(roistats.filename)
                roistats.averageTvalue=readStatsTable(roistats.averageTvalue.filename)
                roistats.averageCoefficientValue=readStatsTable(roistats.averageCoefficientValue.filename)

                roistats$Sub.brick=NULL
                roistats.averageTvalue$Sub.brick=NULL
                roistats.averageCoefficientValue$Sub.brick=NULL                
                
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
                    subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", groups, seedName, "csv", sep="."))
                    subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))
                    
                    cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

                    ## we branch on the regreeeion variable name because the MASC
                    ## tscores used are not in the demographics df but are
                    ## computed a few lines below from the MASC.total. By the time
                    ## the model df is created the MASC.tscore variable will be in
                    ## the mgd df as a result of calling computeMascScore
                    if ( rvVariable == "MASC.tscore" ) {
                        mgd=cbind(subjectOrder, roistats,
                            demographics[match(subjectOrder$subject, demographics$ID), c("Group", "Gender", "DOB", "MRI", "MASC.total")],
                            medsAndTreatment[match(subjectOrder$subject, medsAndTreatment$SubjID), c("Medication", "Treatment")])
                    } else if (grepl("diff", rvVariable) ) {

                        mgd=cbind(subjectOrder, roistats,
                            demographics[match(subjectOrder$subject, demographics$ID), c("Group", "Gender", "DOB", "MRI")],
                            medsAndTreatment[match(subjectOrder$subject, medsAndTreatment$SubjID), c("Medication", "Treatment")])

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
                        mgd=cbind(subjectOrder, roistats,
                            demographics[match(subjectOrder$subject, demographics$ID), c("Group", "Gender", "DOB", "MRI", rvVariable)],
                            medsAndTreatment[match(subjectOrder$subject, medsAndTreatment$SubjID), c("Medication", "Treatment")])
                    }
                    ## ensure that subject is s factor
                    mgd$subject=as.factor(mgd$subject)
                    mgd=droplevels(mgd)
                    
                    if ( rvVariable == "MASC.tscore" ) {
                        mgd=fixDates(mgd)
                        mgd=computeAge(mgd)
                        
                        mgd=computeMascScore(mgd)
                    }

                    ## print(mgd)
                    mgd=replaceNa(mgd, c("Medication", "Treatment"), inSubjectColumnName="subject", inGroupColumnName="Group", inReplaceNaWith=c("No Rx Info", "No Tx Info"))
                    
                    ## print(clusterWhereAmI)
                    ## print(clusters)
                    ##print(mgd)
                    
                    ##stop("Check the mgd data frame\n")
                    publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd, roistats.averageTvalue, roistats.averageCoefficientValue,
                        inStatColumnName="Average t value", inCoefficientColumnName="Average Coefficient Value", inCom=TRUE)
                    savePublicationTable(publicationTable, publicationTableFilename, TRUE)

                    print(publicationTable)
                    ##stop("Check the publication data frame\n")
                    ##melted.mgd=melt(mgd,  id.vars=c("subject", "Group", "Gender", "DOB", "MRI", "Medication", "Treatment", rvVariable),
                    melted.mgd=melt(mgd,  id.vars=c("subject", "Group", "Gender", "DOB", "MRI", "Treatment", rvVariable),                        
                        measure.vars=paste("Mean_", seq(1, clusterCount), sep=""),
                        variable_name="cluster")
                    
                    melted.mgd$cluster=factor(melted.mgd$cluster,
                        levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                        labels=paste(seq(1, clusterCount), clusterWhereAmI))

                    ##print (melted.mgd)
                    ##print(levels(mgd$Treatment))
                    
                    ##stop("Check the melted mgd data frame\n")
                           
                    graphRegressions(melted.mgd, group.results.dir, rvVariable, rvName, polarity, seedName)
                    
                } ## end of if (clusterCount > 0 ) {
            } else {
                cat("No Clusters,\n\n", file=publicationTableFilename, append=TRUE)
            } ## end of if(file.exists(roistats.filename)) {
        } ## end of for (polarity in c("negative", "positive")) {
    } ## end of for (seed in seeds) {

} ## end of generateGraphs definition


graphRegressions <- function(melted.mgd, group.results.dir, rvVariable, rvName, polarity, seedName) {

    imageDirectory=file.path(group.results.dir, seedName)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }

    for ( level in levels(melted.mgd$cluster) ) {

        ss=droplevels(subset(melted.mgd, cluster==level))
        
        imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%s.%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, rvVariable, polarity))
        cat(paste("*** Creating", imageFilename, "\n"))

        roistats.summary=summarySE(ss, measurevar="value", groupvars=c("cluster"))

        y.axis="value"

        ## cat("*** Medication levels and labels\n")
        ## cat("*** Length of medication levels", length(levels(melted.mgd$Medication)), "\n")
        ## cat("*** Length of medication labels", length(sub("Px", "Medication", sub("Info", "Information", levels(melted.mgd$Medication), fixed=TRUE), fixed=TRUE)), "\n")
        ## print(levels(melted.mgd$Medication))
        ## print(sub("Px", "Medication", sub("Info", "Information", levels(melted.mgd$Medication), fixed=TRUE), fixed=TRUE))
        
        ## cat("*** Treatment levels and labels\n")
        ## cat("*** Length of treatment levels", length(levels(melted.mgd$Treatment)), "\n")
        ## cat("*** Length of treatment labels", length(sub("Tx", "Treatment", sub("Info", "Information", levels(melted.mgd$Treatment), fixed=TRUE), fixed=TRUE)), "\n")
        ## print(levels(melted.mgd$Treatment))
        ## print(sub("Tx", "Treatment", sub("Info", "Information", levels(melted.mgd$Treatment), fixed=TRUE), fixed=TRUE))

        graph=ggplot(ss, aes_string(x=rvVariable, y=y.axis)) +
            stat_smooth(method="rlm", se=FALSE, color="black") +
                geom_point(mapping=aes_string(color="Treatment", shape="Treatment")) +
                    scale_color_brewer(name="Treatment:", palette="Set1",
                                       breaks=levels(ss$Treatment),
                                       labels=sub("Tx", "Treatment", sub("Info", "Information", levels(ss$Treatment), fixed=TRUE), fixed=TRUE)) +
                                           scale_shape_discrete(name="Treatment:",
                                                                breaks=levels(ss$Treatment),
                                                                labels=sub("Tx", "Treatment", sub("Info", "Information", levels(ss$Treatment), fixed=TRUE), fixed=TRUE)) +
                                                                    labs(title = substituteShortLabels(level), x=rvName, y="RSFC (Baseline Z-score)") +
                                                                        my_theme

        ## scale_shape_discrete(name="Medication:",
        ##                      breaks=levels(melted.mgd$Treatment),
        ##                      labels=sub("Rx", "Medication", sub("Info", "Information", levels(melted.mgd$Medication), fixed=TRUE), fixed=TRUE)) +

        
                                 
        
        ## print(graph)
        
        ## ggsave(imageFilename, graph, width=3, height=3, units="in")
        ggsave(imageFilename, graph, units="in")
        ## stop("Check graph\n")
    } ## end of for ( level in levels(roistats.summary$cluster) )

} ## end of graphRegressions




####################################################################################################
### End of functions
####################################################################################################

clust.header = c("Volume", "CM RL", "CM AP", "CM IS", "minRL",
    "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
    "MI RL", "MI AP", "MI IS")

task="restingstate"
usedFwhm="4.2"

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


## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)
demographics=rename(demographics, c("Grp"="Group"))

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

medsAndTreatment.filename=file.path(admin.data.dir, "medsAndTreatment.csv")
medsAndTreatment=readCsvFile(medsAndTreatment.filename, "SubjID")

regressionVariables=list(
    list(variable="CDRSR.diff",            name="Children's Depression Rating Scale\n(Baseline to 3 Months Change)"),
    list(variable="MASC.tscore.diff",      name="Multidimensional Anxiety Scale for Children\n(Standardized; Baseline to 3 Months Change)"),
    list(variable="CGAS.diff",             name="Children's Global Assessment Scale\n(Baseline to 3 Months Change)"),
    list(variable="RADS.Total.Tscore.diff", name="Reynolds Adolescent Depression Scale Total\n(Standardized; Baseline to 3 Months Change)")

    ## list(variable="BDI.diff",             name="Beck Depression Inventory II (A to C Change)"),
    ## list(variable="CDI.diff",             name="Children's Depression Inventory (A to C Change)"),
    
    ##list(variable="MASC.tscore",    name="Multidimensional Anxiety Scale for Children (Standardized)")
    #list(variable="CDRS.tscore",    name="Children's Depression Rating Scale (Standardized)")
    ##list(variable="BDI.II",         name="Beck Depression Inventory II")
    ## list(variable="RADS.Total.Tscore",    name="Reynolds Adolescent Depression Scale Total (Standardized)")
    )

groups="mddOnly"

my.base.size=14
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
        legend.position="bottom",        
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

## seeds=readSeedsFile(file.path(config.data.dir, "juelich_amygdala_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "juelich_bla_amygdala_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt"))
## numberOfSeeds=length(seeds)

## seedFiles=
##     sapply(c("juelich_amygdala_seeds.txt",
##              "juelich_whole_amygdala_seeds.txt",
##              "Harvard-Oxford_amygdala_seeds.txt"),
##            function(xx) {
##                file.path(config.data.dir, xx)
##            })

seedFiles=file.path(config.data.dir, "juelich_amygdala_seeds.txt")

for (seedFile in seedFiles) {

    seeds=readSeedsFile(seedFile)

    regressionCount=1
    for ( regressionVariableCount in 1:length(regressionVariables ) ) {

        ## this is the version of the regressionVariable id that will be
        ## used in the model matrix and its accompaning formula and the
        ## output brik filename
        ##rvName=cleanRegressionVariable( regressionVariables[[regressionVariableCount]]$variable)
        rvVariable=regressionVariables[[regressionVariableCount]]$variable
        rvName=regressionVariables[[regressionVariableCount]]$name

        group.data.dir=normalizePath(file.path(data.dir, paste("Group.data", rvVariable, "withAandC", sep=".")))
        group.results.dir=normalizePath(file.path(data.dir, paste("Group.results", rvVariable, "withAandC", sep=".")))

        publicationTableFilename=file.path(group.results.dir, paste("publicationTable", usedFwhm, groups, "csv", sep="."))
        ## if (file.exists(publicationTableFilename)) {
        ##     file.remove(publicationTableFilename)
        ## }

        generateGraphs(group.data.dir, group.results.dir, rvVariable, rvName, publicationTableFilename, seed.list)

        regressionCount=regressionCount+1
    } ## end of for ( regressionVariableCount in 1:length(regressionVariables ) ) {
}
