rm(list=ls())
graphics.off()

## for the some function
library(car)
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
                                 inStatColumnName="Default Stat Name", inContrastColumnName="Default Contrast Name", inCom=TRUE) {
  hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
  ##print(hemisphere)
  if ( inCom ) {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
  } else {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
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
        
generateGraphs <- function () {

    for (seed in seeds) {
        seedName=getSeedName(seed)

        for (cor.var.itr in seq_along(regressionVariables))  {
            regression.variable=regressionVariables[[cor.var.itr]]$variable
            regression.name=regressionVariables[[cor.var.itr]]$name                        
            cat("*** Running correlations for on", regression.variable, "\n")
            
            publicationTableFilename=file.path(group.results.dir, paste("publicationTable", usedFwhm, groups, regression.variable, seedName, "csv", sep="."))
            if (file.exists(publicationTableFilename)) {
                file.remove(publicationTableFilename)
            }

            statLabel=sprintf("timepointX%s.F", regression.variable)
            ## cat("*** Writing publication table to", publicationTableFilename, "\n")
            
            cat("####################################################################################################\n")
            cat(sprintf("*** Graphing ROIs from the %s seed for the %s groups with the %s variable \n", seedName, groups, regression.variable))
            cat(sprintf("StatLabel=%s, Regression Variable=%s\n", statLabel, regression.variable), file=publicationTableFilename, append=TRUE)
            
            infix=sprintf("fwhm%s.%s.%s.%s.%s.%s", usedFwhm, task, groups, seedName, regression.variable, statLabel)
            ## roiStats.fwhm7.98x8.05x7.32.restingstate.MDD.BDI.II.Total.R_whole_amygdala.3mm.timepointXBDI.II.Total.F.txt
            ## roiStats.fwhm7.98x8.05x7.32.restingstate.MDD.RADS.Total.tscore.choi_ventral_striatum_right.3mm.timepointXRADS.Total.tscore.txt


            ## no : roiStats.fwhm7.98x8.05x7.32.restingstate.MDD.BDI.II.Total.L_whole_amygdala.3mm.timepointXBDI.II.Total.txt: No such file or directory
            ## yes: roiStats.fwhm7.98x8.05x7.32.restingstate.MDD.BDI.II.Total.L_whole_amygdala.3mm.timepointXBDI.II.Total.F.txt
            roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))            
            if (is.f.stat(statLabel) ) {
                roistats.averageFvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageFvalue.txt", infix))
            } else if (is.z.stat(statLabel) ) {
                roistats.averageZvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageZValue.txt", infix))
                roistats.averageContrastValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageContrastValue.txt", infix))
            }
            
            if(file.exists(roistats.filename)) {
            
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
                    dataTableFilename=file.path(group.data.dir, paste("dataTable.rsfc", groups, seedName, regression.variable, "timepoint.A.and.C.txt", sep="."))

                    dataTable=fixDataTable(readDataTable(dataTableFilename))
                    ## print(dataTable)
                    ## stop()
                    cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(dataTable$Subj))))

                    mgd=cbind(dataTable, roistats) ## , demographics[match(dataTable$Subj, demographics$ID), c("Gender")])
                    ## colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], "Gender")
                    ## print(mgd)
                    ## stop()
                    ## if (any(mgd$Subj=="378") ) {
                    ##     mgd[mgd$Subj=="378", "Gender"]="F"
                    ## }
                    rownames(mgd)=NULL
                    ## ensure that subject is s factor
                    mgd$Subj=as.factor(mgd$Subj)
                    mgd=droplevels(mgd)
                
                    ## print(clusterWhereAmI)
                    ## print(clusters)
                    ## print(roistats)
                    ## print(mgd)
                    cat("*** Some of the mgd data frame\n")
                    print(mgd)
                    ## stop("Check the mgd data frame\n")
                    print(addmargins(table(mgd[ , c("Group", "timepoint")])))
                    ## stop()
                    if (is.f.stat(statLabel) ) {
                        summaryColumns=unlist(strsplit(gsub(".F", "", statLabel, fixed=TRUE), "X", fixed=TRUE))
                        cat("*** The following variables will be used for the summary columns in the publication table:", paste(summaryColumns, collapse=", "), "\n")
                        publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd,
                            roistats.averageFvalue,
                            inStatColumnName="Average F value", inCom=TRUE)

                        cat("*** Publication table\n")
                        print(publicationTable)
                        ## savePublicationTable(publicationTable, publicationTableFilename, TRUE)
                        ## stop("Check the publication data frame\n")
                        ## stop()
                        ## melted.mgd=melt(mgd,  id.vars=c("Subj", summaryColumns),
                        ##     measure.vars=paste("Mean_", seq.int(1, clusterCount), sep=""),
                        ##     variable_name="cluster")
                        
                        ## melted.mgd$cluster=factor(melted.mgd$cluster,
                        ##     levels=c(paste("Mean_", seq.int(1, clusterCount), sep="")),
                        ##     ## ensure that the sequence number is 2 digits
                        ##     ## and padded with 0 if necessary. Not as
                        ##     ## elegent as zip in python but achieves the
                        ##     ## same result
                        ##     apply(cbind(seq.int(1, length(clusterWhereAmI)), clusterWhereAmI), 1, function(xx) { sprintf("%02d %s", as.integer(xx[1]), xx[2]) } ))
                        
                        ## cat("*** Some of the melted mgd data frame\n")
                        ## print (some(mgd))
                        ## stop("Check the melted mgd data frame\n")

                        graph.f.stats(mgd, groups, seedName, statLabel, regression.variable, regression.name, clusterWhereAmI)
                        
                    }
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


graph.f.stats <-function(inRoistats, inGroups, inSeed, inStatLabel, in.rv, in.rname, clusterWhereAmI) {

    imageDirectory=file.path(group.results.dir, inSeed, inStatLabel, inGroups)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory, recursive=TRUE)
    }

    mean.columns=grep("Mean", colnames(inRoistats))
    n.mean.column=length(mean.columns)
    mean.column.count=1
    for ( roi.column in mean.columns) {

        roi.number=as.numeric(sub("Mean_", "", colnames(inRoistats)[roi.column]))
        imageFilename=file.path(imageDirectory, sprintf("%02d.%s.fwhm%s.%s.%s.%s.%s.pdf", mean.column.count, gsub(" +", ".", clusterWhereAmI[roi.number]),  usedFwhm, task, inGroups, inSeed, inStatLabel))
        cat(paste("*** Creating", imageFilename, "\n"))

        x.axis=in.rv
        y.axis=colnames(inRoistats)[roi.column]
        shape="timepoint"
        color="timepoint"
        group="timepoint"
        xlabel=in.rname


        ## print(roistats.summary)
        my.dodge=position_dodge(.2)
        ## print(sprintf("*** Got to graph.f.stats %s\n", level))
        ## this works for time point or group
        graph = ggplot(data=inRoistats, aes_string(x=x.axis, y=y.axis, color=color, shape=shape, group=group) )
        graph = graph + geom_point(position=my.dodge, size=2.5)
        graph = graph + geom_smooth(method="lm", se=FALSE)
        graph = graph + scale_shape_discrete(name="Timepoint:") 
        graph = graph + scale_color_brewer(name="Timepoint:", palette="Set1")
        graph = graph + labs(title = substituteShortLabels(clusterWhereAmI[roi.number]), x=xlabel, y="RSFC (Z-score)") 
        graph = graph + my.theme
        graph = graph + theme(legend.position="bottom")
     
        ## print(graph)
        ## stop()
        ## ggsave(imageFilename, graph, width=3, height=3, units="in")
        ggsave(imageFilename, graph, units="in")

        mean.column.count = mean.column.count + 1
    } ## end of for ( level in levels(roistats.summary$cluster) )
} ## end of graphFStats


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
group.data.dir=normalizePath(file.path(data.dir, "Group.data.followup"))
group.results.dir=normalizePath(file.path(data.dir, "Group.results.followup", "new_lmes"))

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

##################################################


##################################################
## analysis="all.timepoints"
## analysis="any.timepoint"
## analysis.methods=c("all.timepoints", "any.timepoint")
analysis.methods=c("all.timepoints")

##################################################
## seedFiles=
##     sapply(c("short_ACC_seed_list.txt",
##              "juelich_whole_amygdala_seeds.txt",
##              "hippocampus_ventral_striatum_seeds.txt",
##              "followup-dlpfc-ins-IP-MPFC-seeds.txt"),
##            function(xx) {
##                file.path(config.data.dir, xx)
##            })
## seeds=unlist(sapply(seedFiles, readSeedsFile))

seeds=readSeedsFile(file.path(config.data.dir, "followup-MPFC-seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "hippocampus_ventral_striatum_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt"))

numberOfSeeds=length(seeds)

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

####################################################################################################
####################################################################################################
####################################################################################################

output.filename=file.path(group.results.dir, paste("clinical.lmes.results.output", format(Sys.time(), "%Y%m%d-%H%M%Z"), "txt", sep="."))
cat("*** Output table is in ", output.filename, "\n")
ff=file(output.filename, open="w", encoding="utf-8")
sink(ff, append=FALSE)

if (do.graph.creation) {

    groups="MDD"
    regressionVariables=list(
        list(variable="MASC.tscore",        name="MASC"),
        list(variable="BDI.II.Total",       name="BDI Inventory II"),
        list(variable="CDI.Total",          name="CDI"),
        list(variable="CDRS.t.score",       name="CDRS-R"),
        list(variable="CGAS",               name="CGAS"),
        list(variable="RADS.Total.tscore",  name="RADS")
    )
  
    cat("****************************************************************************************************\n")
    generateGraphs()    
}

## {MASC.tscore,BDI.II.Total,CDI.Total,CDRS.t.score,CGAS,RADS.Total.tscore}
sink()
