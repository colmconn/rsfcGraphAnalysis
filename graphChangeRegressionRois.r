rm(list=ls())
##graphics.off()

##library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
library(grid)
library(plyr)
library(scales)

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

makePublicationTable <- function(inClusterWhereAmI, inClusters, inRoistats,
                                 inRoistats.averageStatValue=NULL, inRoistats.averageCoefficientValue=NULL, inRoistats.averageBiasValue=NULL,
                                 inStatColumnName="Default Stat Name", inCoefficientColumnName="Default Coefficient Name", inBiasColumnName="Default Bias Name",
                                 inCom=TRUE) {

    hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
    ## print(hemisphere)
    if ( inCom ) {
        locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
    } else {
        locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
    }

    ## cat("Locations: Volume and coordinates\n")
    ## print(locations)

    if (length(grep("Mean", colnames(inRoistats))) == 0 ) {
        stop("*** Cat there are no columns in the inRoistats data fram e that begin with Mean_\n")
    }
    
    ## cat ("Columns matching Mean: ", grep("Mean", colnames(inRoistats)), "\n")
    ## cat ("Data from the above columns:\n")
    ## print(inRoistats[, grep("Mean", colnames(inRoistats))])
    ## print(list(inRoistats$Group))
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
        colnames(pubTable) = c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", inStatColumnName)
        ## print(pubTable)
    } else {
        pubTable=locations
        colnames(pubTable) = c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS")        
    }

    if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageCoefficientValue) ) {
        ## cat("Adding average coefficient values\n")      
        pubTable=cbind(pubTable, round(t(inRoistats.averageCoefficientValue), 2))
        colnames(pubTable)[dim(pubTable)[2]] = inCoefficientColumnName
        ## print(pubTable)      
    }

    ## only add the bias column if the stats and coefficient arguments are also not null
    if (! is.null(inRoistats.averageStatValue) & ! is.null(inRoistats.averageCoefficientValue) & ! is.null (inRoistats.averageBiasValue) ) {
        pubTable=cbind(pubTable, round(t(inRoistats.averageBiasValue), 4))
        colnames(pubTable)[dim(pubTable)[2]] = inBiasColumnName
    }

    pubTable=cbind(pubTable, mns)
    
    rownames(pubTable)=NULL
    ## print(pubTable)
    return(pubTable)
}


savePublicationTable <- function (inPublicationTable, inPublicationTableFilename, append=TRUE) {
    cat("*** Writing publication table to", inPublicationTableFilename, "\n")
    if ( append ) {
        write.table(inPublicationTable, file=inPublicationTableFilename, quote=F, col.names=TRUE, row.names=FALSE, sep=",", append=TRUE)
    } else {
        write.table(inPublicationTable, file=inPublicationTableFilename, quote=F, col.names=FALSE, row.names=FALSE, sep=",", append=FALSE)
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
    clusterWhereAmI=gsub(" +", " ", scan(file=inFilename, what='character', sep=',', quiet=TRUE))

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


read.change.score.file <- function (in.variable) {

    if (grepl("both|rstandard", rvVariable, fixed=FALSE))
        change.score.filename=file.path(admin.data.dir, sprintf("new.mdd.%s.scores.csv", in.variable))
    else if (baselineOnly)
        change.score.filename=file.path(admin.data.dir, sprintf("new.mdd.%s.timepoint.a.score.csv", in.variable))
    else if (grepl("diff", rvVariable, fixed=TRUE))
        change.score.filename=file.path(admin.data.dir, sprintf("new.mdd.%s.change.score.csv", in.variable))
    else
        change.score.filename=file.path(admin.data.dir, sprintf("new.mdd.%s.scores.csv", in.variable))        
        ##change.score.filename=file.path(admin.data.dir, sprintf("new.mdd.%s.rstandard.score.csv", in.variable))        
    
    cat("*** Change score file is", change.score.filename, "\n")
    if (file.exists(change.score.filename))
        change.scores=readCsvFile(change.score.filename)
    else
        stop("Couldn't find the change score file to go with a variable named ", in.variable, "\n")

    return(change.scores)
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

                if ( ! replacement.value %in% levels(inData[, column]) ) {
                    inData[, column] = as.factor(inData[, column])
                }
            }
            cat ("****************************************************************************************************\n")      
        }
    } ## end of for (column in inColumns) {
    
    return(inData)
} ## end of checkIsNa


generateGraphs <- function (group.data.dir, group.results.dir, rvVariable, rvName, publicationTableFilename, seed.list, bootstrapped=FALSE) {


    if (groups == "mddAndCtrl" && grepl("diff", rvVariable)) {
        print("*** You are trying to regress/graph  the control and MDD subjects aagainst a time A to C difference. you cannot do that!\n")
        return()
    }
    
    for (seed in seed.list) {
        seedName=getSeedName(seed)

        cat(sprintf("Seed=%s,Variable=%s (%s),\n", seedName, gsub("\n", " ", rvName), rvVariable), file=publicationTableFilename, append=TRUE)
        
        cat("####################################################################################################\n")
        cat(sprintf("*** Graphing ROIs for the %s seed that was regressed against %s (%s) in the %s group\n", seedName, rvName, rvVariable, groups))

        infix=sprintf("regression.fwhm%s.%s.%s.%s.and.%s", usedFwhm, task, groups, seedName, rvVariable)

        roistats.filename=file.path(group.results.dir, sprintf("roiStats.%s.txt", infix))
        roistats.averageTvalue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageTValue.txt", infix))
        roistats.averageCoefficientValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageCoefficientValue.txt", infix))
        if (bootstrapped) {
            roistats.averageBiasValue.filename=file.path(group.results.dir, sprintf("roiStats.%s.averageBiasValue.txt", infix))
        }

            
        if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the contrast in each ROI,
            ## you do not what to graph this
            
            roistats=readStatsTable(roistats.filename)
            roistats.averageTvalue=readStatsTable(roistats.averageTvalue.filename)
            roistats.averageCoefficientValue=readStatsTable(roistats.averageCoefficientValue.filename)
            if (bootstrapped)
                roistats.averageBiasValue=readStatsTable(roistats.averageBiasValue.filename)
            
            roistats$Sub.brick=NULL
            roistats.averageTvalue$Sub.brick=NULL
            roistats.averageCoefficientValue$Sub.brick=NULL
            if (bootstrapped)
                roistats.averageBiasValue$Sub.brick=NULL
            
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

                change.scores=read.change.score.file(rvVariable)

                mgd=cbind(subjectOrder, roistats)
                mgd=merge(mgd, change.scores, by.x="subject", by.y="ID", sort=FALSE)

                if (! all ( mgd$subject == subjectOrder$subject) ) {
                    stop("*** The order of the subjects in the mgd data frame is not the same as the order of subjects in subjecOrder.\n",
                         "*** Cannot continue.\n")
                }
                
                if (dim(mgd)[1] != dim(subjectOrder)[1] ) {
                    cat("*** The number of subjects in the merged data frame is not the same as the number of subjects in the subjectOrder file.\n")
                    cat("*** Cannot continue\n")
                    stop(status=1)
                }

                ## print(head(subjectOrder))
                ## print(head(mgd))

                ss.medsAndTreatment=subset(medsAndTreatment, Timepoint=="C", select=c("SubjectNumber", "SubjID", "Timepoint", "Group", "Rx.Summary", "Tx.Summary"))
                ## print(ss.medsAndTreatment)
                
                ## now try to add the meds and treatment info for each subject
                ## mgd=cbind(mgd,
                ##     medsAndTreatment[match(subjectOrder$subject, medsAndTreatment$SubjID), c("Medication", "Treatment")])
                mgd=cbind(mgd,
                    ss.medsAndTreatment[match(subjectOrder$subject, ss.medsAndTreatment$SubjID), c("Rx.Summary", "Tx.Summary")])
                rownames(mgd)=NULL

                mgd=rename(mgd, c("Grp"="Group"))
                
                if (any(grepl("Tx.Summary", colnames(mgd), fixed=TRUE))) {
                    mgd=rename(mgd, c("Tx.Summary"="Treatment"))
                }
                if (any(grepl("Rx.Summary", colnames(mgd), fixed=TRUE))) {                
                    mgd=rename(mgd, c("Rx.Summary"="Medication"))
                }
                
                ## print(head(mgd))
                ## print(table(mgd[, c("Medication")]))
                ## print(table(mgd[, c("Treatment")]))                
                      
                mgd=replaceNa(mgd, c("Medication", "Treatment"), inSubjectColumnName="subject", inGroupColumnName="Group", inReplaceNaWith=c("No Rx Info", "No Tx Info"))
                ## stop()
                
                ## cat("*** subjectOrder:\n")
                ## print(subjectOrder)
                ## cat("*** clusterWhereAmI:\n")
                ## print(clusterWhereAmI)
                ## cat("*** clusters:\n")
                ## print(clusters)
                ## cat("*** mgd:\n")                
                ## print(mgd)
                ## cat("*** roistats.averageTvalue:\n")                
                ## print(roistats.averageTvalue)
                ## cat("*** roistats.averageCoefficientValue:\n")                
                ## print(roistats.averageCoefficientValue)

                ## print(head(mgd))
                ## stop("Check the mgd data frame\n")
                if (bootstrapped) {
                    publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd,
                        inRoistats.averageStatValue=roistats.averageTvalue,
                        inRoistats.averageCoefficientValue=roistats.averageCoefficientValue,
                        inRoistats.averageBiasValue=roistats.averageBiasValue,
                        inStatColumnName="Average t value",
                        inCoefficientColumnName="Average Coefficient Value",
                        inBiasColumnName="Average Bias Value",
                        inCom=TRUE)
                } else {
                    publicationTable=makePublicationTable(clusterWhereAmI, clusters, mgd,
                        inRoistats.averageStatValue=roistats.averageTvalue, inRoistats.averageCoefficientValue=roistats.averageCoefficientValue,
                        inStatColumnName="Average t value", inCoefficientColumnName="Average Coefficient Value",
                        inCom=TRUE)
                }
                savePublicationTable(publicationTable, publicationTableFilename, TRUE)

                print(publicationTable)
                ## stop("Check the publication data frame\n")
                if ( grepl ("both.scaled", rvVariable, fixed=TRUE) ) {
                    melted.mgd=melt(mgd,  id.vars=c("subject", "Group", "Treatment", "Medication", sub("both.scaled", "C.scaled", rvVariable, fixed=TRUE)),
                        measure.vars=paste("Mean_", seq(1, clusterCount), sep=""),
                        variable_name="cluster")
                    graph.variable=sub("both.scaled", "C.scaled", rvVariable, fixed=TRUE)
                } else if ( grepl ("both", rvVariable, fixed=TRUE) ) {
                    melted.mgd=melt(mgd,  id.vars=c("subject", "Group", "Treatment", "Medication", sub("both", "C", rvVariable, fixed=TRUE)),
                        measure.vars=paste("Mean_", seq(1, clusterCount), sep=""),
                        variable_name="cluster")
                    graph.variable=sub("both", "C", rvVariable, fixed=TRUE)
                } else {
                    ## melted.mgd=melt(mgd,  id.vars=c("subject", "Group", "Gender", "DOB", "MRI", "Medication", "Treatment", rvVariable),
                    melted.mgd=melt(mgd,  id.vars=c("subject", "Group", rvVariable),                        
                        measure.vars=paste("Mean_", seq(1, clusterCount), sep=""),
                        variable_name="cluster")
                    graph.variable=rvVariable
                }
                
                melted.mgd$cluster=factor(melted.mgd$cluster,
                    levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                    labels=paste(sprintf("%02d", seq(1, clusterCount)), clusterWhereAmI))

                ## print (head((melted.mgd)))
                ## print(levels(mgd$Treatment))

                ## print(mgd)
                ## stop("Check the melted mgd data frame\n")
                
                graphRegressions(melted.mgd, group.results.dir, graph.variable, rvName, seedName, indicate.treatments=FALSE)
                
            } ## end of if (clusterCount > 0 ) {
        } else {
            cat ("*** No such file", roistats.filename, "\n")
            cat("No Clusters,\n\n", file=publicationTableFilename, append=TRUE)
        } ## end of if(file.exists(roistats.filename)) {

    } ## end of for (seed in seeds) {

} ## end of generateGraphs definition


graphRegressions <- function(melted.mgd, group.results.dir, rvVariable, rvName, seedName, indicate.treatments=FALSE) {

    ## print(head(melted.mgd))
    imageDirectory=file.path(group.results.dir, seedName)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }

    for ( level in levels(melted.mgd$cluster) ) {

        ss=droplevels(subset(melted.mgd, cluster==level))
        ## print(ss)
        ## stop()
        
        imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, rvVariable))
        cat(paste("*** Creating", imageFilename, "\n"))

        roistats.summary=summarySE(ss, measurevar="value", groupvars=c("cluster"))
        
        if (grepl("\\.reversed$", group.results.dir)) {
            cat("*** Using reversed axis variables and variable names\n")
            x.axis="value"
            y.axis=rvVariable
            x.axis.label="RSFC (Z-score)"
            y.axis.label=rvName
        } else if (grepl("both", group.results.dir)) {
            x.axis="value"
            y.axis=rvVariable
            x.axis.label="RSFC (Z-score)"
            y.axis.label=rvName
        } else {
            x.axis="value"
            y.axis=rvVariable
            x.axis.label="RSFC (Z-score)"
            y.axis.label=rvName
        }

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

        graph=ggplot(ss, aes_string(x=x.axis, y=y.axis, label="subject")) +
            ## stat_smooth(method="rlm", se=FALSE, color="black") +
                labs(title = substituteShortLabels(level), x=x.axis.label, y=y.axis.label) +
                    my_theme
        if (grepl("scaled", group.results.dir, fixed=TRUE)) {
            graph=graph + scale_y_continuous(labels = percent)
        } ## else { 
        ##     graph=graph + scale_x_continuous(labels = percent)
        ## }

        if ( indicate.treatments ) {
            cat ("*** Adding treatment indicators to graphs\n")
            ## print(head(ss))
            ## geom_point(mapping=aes_string(color="Treatment", shape="Treatment")) +
            ## geom_point(mapping=aes_string(color="Treatment", shape="Medication")) +
            ## colors generated here http://phrogz.net/css/distinct-colors.html
            ## my.colors = c("#8c0000", "#ff8080", "#a67c7c", "#b7ff00", "#5f664d", "#33cc8a", "#0024f2", "#131c4d", "#a3abd9", "#a62995")

            ## from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
            ## The palette with black:
            cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

            graph = graph +
                geom_point(mapping=aes_string(color="Treatment", shape="Medication")) +
                    scale_color_manual(name="Treatment:",
                                       guide = guide_legend(nrow=3),
                                       breaks=levels(ss$Treatment),
                                       labels=sub("Tx", "Treatment", sub("Info", "Information", levels(ss$Treatment), fixed=TRUE), fixed=TRUE),
                                       values=cbbPalette
                                       )

                    
                ## scale_color_brewer(name="Treatment:", palette="Dark2",
                ##                    guide = guide_legend(nrow=3),
                ##                    breaks=levels(ss$Treatment),
                ##                    labels=sub("Tx", "Treatment", sub("Info", "Information", levels(ss$Treatment), fixed=TRUE), fixed=TRUE)
                ##                    )
            graph = graph +
                scale_shape_discrete(name="Medication:",
                                     guide = guide_legend(nrow=3),
                                     breaks=levels(melted.mgd$Medication),
                                     labels=capwords(sub("Rx", "Medication", sub("Info", "Information", levels(melted.mgd$Medication), fixed=TRUE), fixed=TRUE)))
                ## scale_shape_discrete(name="Treatment:",
                ##                      guide = guide_legend(nrow=2),                                                                    
                ##                      breaks=levels(ss$Treatment) ##,
                ##                      ## labels=sub("Tx", "Treatment", sub("Info", "Information", levels(ss$Treatment), fixed=TRUE), fixed=TRUE)
                ##                      )
        } else {
            graph= graph + geom_point()
        }
                                                                
        ## scale_shape_discrete(name="Medication:",
        ##                      breaks=levels(melted.mgd$Treatment),
        ##                      labels=sub("Rx", "Medication", sub("Info", "Information", levels(melted.mgd$Medication), fixed=TRUE), fixed=TRUE)) +

        
                                 
        
        ## print(graph)
        ## stop()
        ggsave(imageFilename, graph, width=4, height=3, units="in")
        ## ggsave(imageFilename, graph, width=7, height=7, units="in")
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

data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
seeds.data.dir=file.path(data.dir, "seeds")


## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
## demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
## demographics=readCsvFile(demographicsFilename)
## demographics=rename(demographics, c("Grp"="Group"))

## medsAndTreatment.filename=file.path(admin.data.dir, "medsAndTreatment.csv")
medsAndTreatment.filename=file.path(admin.data.dir, "new_medsAndTreatment.csv")
medsAndTreatment=readCsvFile(medsAndTreatment.filename, "SubjID")

regressionVariables=list(
    ## list(variable="CDRS.t.score.diff",             name="Children's Depression Rating Scale\n(Baseline to 3 Months Change)"),
    ## list(variable="MASC.tscore.diff",              name="Multidimensional Anxiety Scale for Children\n(Standardized; Baseline to 3 Months Change)"),
    ## list(variable="CGAS.diff",                     name="Children's Global Assessment Scale\n(Baseline to 3 Months Change)"),
    ## list(variable="RADS.Total.tscore.diff",        name="Reynolds Adolescent Depression Scale Total\n(Standardized; Baseline to 3 Months Change)"),
    
    ## list(variable="CDRS.t.score.scaled.diff",      name="Children's Depression\nRating Scale\n(Baseline to 3 Months Change)")#,

    ## list(variable="CDRS.t.score.scaled",           name="Children's Depression\nRating Scale (Baseline)"),

    ## baseline 
    ## list(variable="CDRS.t.score",                     name="Children's Depression\nRating Scale (Baseline)"),    

    ## predictive regressions
    list(variable="CDRS.t.score.scaled.diff",      name=expression(paste(Delta, " CDRS-R")))#,

    ## list(variable="CDRS.t.score.rstandard",        name="CDRS-R Residual")#,

    ## list(variable="CDRS.t.score.both.scaled",         name="Follow-up CDRS-R")#,    
    ## list(variable="CDRS.t.score.both",                name="Follow-up CDRS-R")#,
    
    ## list(variable="CGAS.scaled.diff",              name="Children's Global Assessment Scale\n(Baseline to 3 Months Change)")
    
    ## list(variable="MASC.tscore.scaled.diff",       name="Multidimensional Anxiety Scale for Children\n(Standardized; Baseline to 3 Months Change)"),
    ## list(variable="RADS.Total.tscore.scaled.diff", name="Reynolds Adolescent Depression Scale Total\n(Standardized; Baseline to 3 Months Change)")
    
    ## list(variable="BDI.II.Total.diff",             name="Beck Depression Inventory II (A to C Change)"),
    ## list(variable="CDI.Total.diff",                name="Children's Depression Inventory (A to C Change)")

    ## list(variable="BDI.II.Total",             name="Beck Depression Inventory II")
    ## list(variable="MASC.tscore",              name="Multidimensional Anxiety Scale\nfor Children (Standardized)")    
)

groups="mddOnly"

my.base.size=14
my_theme=
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

## seeds=readSeedsFile(file.path(config.data.dir, "juelich_amygdala_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "juelich_bla_amygdala_seeds.txt"))
## seeds=readSeedsFile(file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt"))
## numberOfSeeds=length(seeds)

## seedFiles=
##     sapply(c("juelich_amygdala_seeds_weights.txt",
##              "juelich_whole_amygdala_seeds.txt"),
##            function(xx) {
##                file.path(config.data.dir, xx)
##            })

##seedFiles=file.path(config.data.dir, "juelich_amygdala_seeds_weights.txt")
seedFiles=file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt")

baselineOnly=FALSE

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

        group.data.dir=file.path(data.dir, paste("Group.data", rvVariable, sep="."))
        group.results.dir=file.path(data.dir, paste("Group.results", rvVariable, sep="."))##, "wholeBrain")
        
        publicationTableFilename=file.path(group.results.dir, paste("publicationTable", usedFwhm, groups, "csv", sep="."))
        ## if (file.exists(publicationTableFilename)) {
        ##     file.remove(publicationTableFilename)
        ## }

        generateGraphs(group.data.dir, group.results.dir, rvVariable, rvName, publicationTableFilename, seeds, bootstrapped=TRUE)

        regressionCount=regressionCount+1
    } ## end of for ( regressionVariableCount in 1:length(regressionVariables ) ) {
}
