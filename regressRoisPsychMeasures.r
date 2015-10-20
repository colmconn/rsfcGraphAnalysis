rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
library(boot)

source("scoreMasc.r")

##########################################################################################################################################################################
### START OF FUNCTIONS ###################################################################################################################################################
##########################################################################################################################################################################

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
    returnSubstitutedLabel = gsub("[0-9]+ ", "", gsub("Inf", "Inferior", gsub("Sup", "Superior", gsub("Gy", "Gyrus", gsub("R", "Right", gsub("L", "Left", inLevel, fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE))
    
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

makeGraphTitle <- function(inLevel) {

    ## graphTitle= gsub("^[0-9]+[ ]*", "", gsub("Inf", "Inferior", gsub("Gy", "Gyrus", gsub("Sup", "Superior", inLevel))))
    graphTitle= gsub("Inf", "Inferior", gsub("Gy", "Gyrus", gsub("Sup", "Superior", inLevel)))    

    return(graphTitle)
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


readStatsTable <- function (inFilename) {

    cat("*** Reading" , inFilename, "\n")
    statsTable=read.table(inFilename, header=T, sep="")
    ## dump the first column as it's only the file name
    statsTable=statsTable[, -c(1, 2)]
    return(statsTable)
}

readClustersTable <- function (inFilename){
    cat("*** Reading", file.path(group.results.dir, inFilename), "\n")
    clusters=read.table(file.path(group.results.dir, inFilename))
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
    cat("*** Reading", subjectOrderFilename, "\n")
    subjectOrder=read.csv(subjectOrderFilename, header=T)

    return(subjectOrder)
}

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
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

computeMascScore <- function (inData) {
    inData.dim=dim(inData)
    for (r in seq(1, inData.dim[1]) ) {
        ## cat("##################################################\n")
        subjectNumber=inData[r, "subject"]
        gender=inData[r, "Gender"]
        age=round(inData[r, "age.in.years"], 0)
        old.masc.tscore=inData[r, "MASC.tscore"]

        new.masc.tscore=scoreMasc(gender, age, inData[r, "MASC.total"])
        if (is.na(new.masc.tscore) ) {
            warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, inData[r, "MASC.total"]))
        }
        
        inData[r, "MASC.tscore"]=new.masc.tscore

        ## cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f new MASC tscore=%0.0f\n",
        ##             r, subjectNumber, gender, age, inData[r, "MASC.total"], old.masc.tscore, new.masc.tscore))

    }
    return (inData)
}

runRegressions <- function (inSeed, inDataFrame, inRegressionVariables) {
    for ( j in 1:length(inRegressionVariables ) ) {
        ## clCounter=clusterCounter
        clCounter=1
        for ( level in levels(inDataFrame$cluster) ) {
            cat("####################################################################################################\n")
            cat("*** Seed: ", inSeed, "\n")
            cat("*** Running regression model for ", level, "on", inRegressionVariables[[j]]$name, "\n")
            
            regressionFormula=as.formula(paste("value",  "~", inRegressionVariables[[j]]$variable, sep=" "))
            cat("*** Regression formula: ")
            print(regressionFormula)


            df.all=   inDataFrame[inDataFrame$cluster %in% level, ]

            df.model=df.all[, c("value", inRegressionVariables[[j]]$variable, "Gender")]
            regression.formula=as.formula(paste("value", "~", inRegressionVariables[[j]]$variable, sep=" "))
            lm.model.nogender=lm(regression.formula, data=df.model)
            print(summary(lm.model.nogender))


            regression.formula=as.formula(paste("value", "~", inRegressionVariables[[j]]$variable, "+", "Gender", sep=" "))
            lm.model.withgender=lm(regression.formula, data=df.model)
            print(summary(lm.model.withgender))
            
            clCounter=clCounter+1
            cat("clCounter is now", clCounter, "\n")
        } ## end of for ( level in levels(inDataFrame$cluster) )
        
        ##stop("Stopping")
    } ## end of for ( level in levels(roistats.summary$cluster) )
    
}


runRegressionsForKajaWithGrpInModel <- function (inSeed, inDataFrame, inRegressionVariables, csv.file.name=NULL, append=TRUE) {

    regression.stack=stack()

    inDataFrame$Grp=relevel(inDataFrame$Grp, ref="NCL")
    
    for ( j in 1:length(inRegressionVariables ) ) {
        ## clCounter=clusterCounter
        clCounter=1
        for ( level in levels(inDataFrame$cluster) ) {
            cat("####################################################################################################\n")
            cat("*** Seed: ", inSeed, "\n")
            cat("*** Running regression model for ", level, "on", inRegressionVariables[[j]]$name, "(", inRegressionVariables[[j]]$variable, ")\n")

            df.all=   inDataFrame[inDataFrame$cluster %in% level, ]
            ## print(head(df.all))
            df.model=df.all[, c("value", inRegressionVariables[[j]]$variable, "Grp")]

            regression.formula=as.formula(paste("value",  "~", inRegressionVariables[[j]]$variable, "+ Grp", sep=" "))
            ## cat("*** Regression formula: ")
            ## print(regression.formula)

            lm.model=lm(regression.formula, data=df.model)
            sm=summary(lm.model)
            ## print(sm)

            coeff=coefficients(sm)
            
            for ( term in rownames(coeff) ) {
                csv.line=
                    sprintf("%s,%s,%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%s",
                            getSeedName(inSeed), level, term, inRegressionVariables[[j]]$variable,
                            coeff[term, "Estimate"],
                            coeff[term, "Std. Error"],
                            coeff[term, "t value"],
                            coeff[term, "Pr(>|t|)"],
                            make.significance.indications(coeff[term, "Pr(>|t|)"])
                            )
               push(regression.stack, csv.line)
            }

            clCounter=clCounter+1
            ## cat("clCounter is now", clCounter, "\n")
        } ## end of for ( level in levels(inDataFrame$cluster) )

        ##stop("Stopping")
    } ## end of for ( level in levels(roistats.summary$cluster) )
    
    if ( ! file.exists(csv.file.name)  ) {
        header="Seed,cluster,term,Regression Variable,estimate,std. error,t value,p value,significance"
        cat(header, "\n", file=csv.file.name)
    }
    if (is.null(csv.file.name) ) {
        header="Seed,cluster,term,Regression Variable,estimate,std. error,t value,p value,significance"
        cat(header, "\n")
    }
    
    l=regression.stack$value()
    if (! is.null(csv.file.name)) {
        sink(csv.file.name, append=append)
    }
   
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
    
    if ( ! is.null(csv.file.name)) {
        sink()
    }
}

runRegressionsForKajaWithGrpSeparated <- function (inSeed, inDataFrame, inRegressionVariables, csv.file.name=NULL, append=TRUE) {

    regression.stack=stack()

    inDataFrame$Grp=relevel(inDataFrame$Grp, ref="NCL")
    
    for ( j in 1:length(inRegressionVariables ) ) {
        ## clCounter=clusterCounter
        clCounter=1
        for ( level in levels(inDataFrame$cluster) ) {
            for ( gg in levels(inDataFrame$Grp) ) {
                cat("####################################################################################################\n")
                cat("*** Seed: ", inSeed, "\n")
                cat("*** Running regression model for ", level, "on", inRegressionVariables[[j]]$name, "(", inRegressionVariables[[j]]$variable, ") in the", gg, "group\n")

                df.all=   inDataFrame[inDataFrame$cluster %in% level & inDataFrame$Grp == gg, ]
                ## print(head(df.all))
                df.model=df.all[, c("value", inRegressionVariables[[j]]$variable, "Grp")]

                regression.formula=as.formula(paste("value",  "~", inRegressionVariables[[j]]$variable, sep=" "))
                ## cat("*** Regression formula: ")
                ## print(regression.formula)

                lm.model=lm(regression.formula, data=df.model)
                sm=summary(lm.model)
                ## print(sm)

                coeff=coefficients(sm)
            
                for ( term in rownames(coeff) ) {
                    csv.line=
                        sprintf("%s,%s,%s,%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%s",
                                getSeedName(inSeed), level, term, inRegressionVariables[[j]]$variable, gg,
                                coeff[term, "Estimate"],
                                coeff[term, "Std. Error"],
                                coeff[term, "t value"],
                                coeff[term, "Pr(>|t|)"],
                                make.significance.indications(coeff[term, "Pr(>|t|)"])
                                )
                    push(regression.stack, csv.line)
                }
                
                clCounter=clCounter+1
                ## cat("clCounter is now", clCounter, "\n")
            } ## end of for ( level in levels(inDataFrame$cluster) )
            
            ##stop("Stopping")
        } ## end of for ( gg in levels(inDataFrame$Grp) ) {
    } ## end of for ( level in levels(roistats.summary$cluster) )
        
    if ( ! file.exists(csv.file.name)  ) {
        header="Seed,cluster,term,Regression Variable,group,estimate,std. error,t value,p value,significance"
        cat(header, "\n", file=csv.file.name)
    }
    if (is.null(csv.file.name) ) {
        header="Seed,cluster,term,Regression Variable,group,estimate,std. error,t value,p value,significance"
        cat(header, "\n")
    }
    
    l=regression.stack$value()
    if (! is.null(csv.file.name)) {
        sink(csv.file.name, append=append)
    }
   
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
    
    if ( ! is.null(csv.file.name)) {
        sink()
    }
}


runCorrelations <- function (inSeed, inDataFrame, inRegressionVariables, inClusters, inMethod="pearson", results.table.filename=NULL) {

    results.stack=stack()

    for ( j in 1:length(inRegressionVariables ) ) {
        ## clCounter=clusterCounter
        clCounter=1
        for ( level in levels(inDataFrame$cluster) ) {
            cat("####################################################################################################\n")
            cat("*** Seed: ", inSeed, "\n")
            cat("*** Running regression model for ", level, "on", inRegressionVariables[[j]]$name, "\n")
            
            ## regressionFormula=as.formula(paste("value",  "~", inRegressionVariables[[j]]$variable, sep=" "))
            ## cat("*** Regression formula: ")
            ## print(regressionFormula)

            df.all=   inDataFrame[inDataFrame$cluster %in% level, ]

            ct.all=cor.test(df.all[, regressionVariables[[j]]$variable], df.all$value, method=inMethod)
            ##print(ct.all)

            csvLine=sprintf("%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
                getSeedName(seed), "all", makeGraphTitle(level), paste(inClusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.all$statistic,
                ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
                ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
            ## cat(csvLine, "\n")
            push(results.stack, csvLine)
            
            clCounter=clCounter+1
            ## cat("clCounter is now", clCounter, "\n")
        } ## end of for ( level in levels(inDataFrame$cluster) )
        
        ##stop("Stopping")
    } ## end of for ( level in levels(roistats.summary$cluster) )

    if ( ! is.null(results.table.filename)) {
        cat("*** Correlation table is in ", results.table.filename, "\n")
        sink(results.table.filename, append=FALSE)
    } else {
        cat("################################################################################\n");
        cat("Correlation statistics table\n")
    }

    l=results.stack$value()
    header="Seed,gender,cluster,RL,AP,IS,RegressionVariable,S,DoF,pValue,Rho,Significance"
    cat(header, "\n")
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
    if ( ! is.null(results.table.filename)) {
        sink()
    }
}


graphRegressions <- function (inSeed, inDataFrame, inRegressionVariables, inClusters, inMethod="pearson") {

    seedName=getSeedName(inSeed)
    
    for ( j in 1:length(inRegressionVariables ) ) {
        ## clCounter=clusterCounter
        clCounter=1
        for ( level in levels(inDataFrame$cluster) ) {
            cat("####################################################################################################\n")
            cat("*** Seed: ", inSeed, "\n")
            cat("*** Running regression model for ", level, "on", inRegressionVariables[[j]]$name, "\n")

            df.all=   inDataFrame[inDataFrame$cluster %in% level, ]

            ct.all=cor.test(df.all[, regressionVariables[[j]]$variable], df.all$value, method=inMethod)
            if (ct.all$p.value < 0.05 ) {
            
                ## make scatter plots for the ROIs that are significant
                imageDirectory=file.path(group.results.dir, seedName)
                if ( ! file.exists(imageDirectory) ) {
                    dir.create(imageDirectory)
                }
                imageFilename=file.path(imageDirectory,
                    sprintf("%s.fwhm%s.%s.and.%s.pdf",
                            gsub(" +", ".", level),  usedFwhm, seedName, regressionVariables[[j]]$variable))
            
                message(paste("*** Creating", imageFilename, "\n"))
            
                ylabel="RSFC (z-score)"
                graphTitle=makeGraphTitle(level)
                ##sprintf("%s", toupper(seed))
                my.base.size=18

                x.axis=regressionVariables[[j]]$variable
                y.axis="value"
                
                graph=ggplot(df.all, aes_string(x=x.axis, y=y.axis)) +
                    theme_bw(base_size =  my.base.size) +
                        geom_point() +
                            geom_smooth(method="lm", color="black") +
                                scale_fill_brewer(palette="Set1") +
                                    labs(title = graphTitle, y=ylabel, x=regressionVariables[[j]]$name) +
                                        my_theme + theme(legend.position="none")
                
                ggsave(imageFilename, graph)
            } ## end of if (ct.all$p.value < 0.05 ) {
           
            clCounter=clCounter+1
            cat("clCounter is now", clCounter, "\n")
        } ## end of for ( level in levels(inDataFrame$cluster) )
        
        ##stop("Stopping")
    } ## end of for ( level in levels(roistats.summary$cluster) )
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

scripts.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))
data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/"))
admin.data.dir=normalizePath(file.path(data.dir, "admin"))
config.data.dir=normalizePath(file.path(data.dir, "config"))
group.data.dir=normalizePath(file.path(data.dir, "Group.data"))
seeds.data.dir=normalizePath(file.path(data.dir, "seeds"))

group.results.dir=normalizePath(file.path(data.dir, "Group.results"))

####################################################################################################
## These variables control what barcharts are created

createIndividualRegressionImageFiles=FALSE
createSingleRegressionPdfFile=TRUE

##runRegressions=TRUE

####################################################################################################

clust.header = c("Volume", "CM RL", "CM AP", "CM IS", "minRL",
    "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
    "MI RL", "MI AP", "MI IS")

groups="mddAndCtrl"
usedFwhm="4.2"

my.base.size=14
my_theme=
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
        
        ##axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))


##demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_10152013.csv")
demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
slesFilename=file.path(admin.data.dir, "SLES_20140820.csv")

cat("*** Reading", slesFilename, "\n")
sles=read.csv(slesFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
cat(sprintf("*** Read SLES data for %s unique subjects\n",  length(unique(sles$subject))))
sles=sles[, c(1:92)]

cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))


selectedColumns=c(
    "Grp", "Gender", "DOB", "MRI", "Race", "Hand", "SES", "Tanner1", "Tanner2", "TannerAvg",
    "CGI.CGAS",                 "CDRS.raw",                  "CDRS.tscore",               "SFMQ.total",              
    "MADRS.raw",                "WASI.PERF",                 "WASI.Full.4",               "PSWQ", 		"CoRum",
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

regressionVariables=list(
    ## list(variable="CDRSR.diff",             name="Children's Depression Rating Scale\n(Baseline to 3 Months Change)"),
    ## list(variable="MASC.tscore.diff",       name="Multidimensional Anxiety Scale for Children\n(Standardized; Baseline to 3 Months Change)"),
    ## list(variable="CGAS.diff",              name="Children's Global Assessment Scale\n(Baseline to 3 Months Change)"),
    ## list(variable="RADS.Total.Tscore.diff", name="Reynolds Adolescent Depression Scale Total\n(Standardized; Baseline to 3 Months Change)")

    ## list(variable="num.stressful.events", name="SLES: Number of stressful events"),
    ## list(variable="num.severe.events",    name="SLES: Number of severe events"),
    ## list(variable="total.sum.stress",     name="SLES: Total Sum Stress")
    
    ## list(variable="BDI.diff",             name="Beck Depression Inventory II (A to C Change)"),
    ## list(variable="CDI.diff",             name="Children's Depression Inventory (A to C Change)"),
    
    list(variable="MASC.tscore",        name="Multidimensional Anxiety Scale for Children (Standardized)"),
    list(variable="CDRS.tscore",        name="Children's Depression Rating Scale (Standardized)"),
    list(variable="CGI.CGAS",           name="Children's Global Assessment Scale"),
    ## list(variable="RADS.Total.Tscore",  name="Reynolds Adolescent Depression Scale Total (Standardized)")

    ## list(variable="C_Irritability", name="Caprara Irritability Scale"),
    ## list(variable="C_EmoSuscept",   name="Caprara Emotional Susceptibility Scale")

    list(variable="RSQ.fixed", name="Rumination Response Styles Questionnaire")
    )

## ## extract the list of variable names from the regressionVariables list rvs=regression variables
rvs=unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))]
## ## select only the columns we want to perform regressions on from the hnrcData data frame
##m=match(unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))], colnames(demographics))
## ## remove NAs caused by the BDI and POMS not coming from the hnrcData frame
##m=m[!is.na(m)]
## m now contains only the column numbers of the neuropsych variables we'll regress the %cs against

## this should be set to the number of subjects you have in the bucket
## files and hence the subjectOrder file.
expectedNumberOfSubjects=101

## seeds=readSeedsFile("../data/config/juelich_amygdala_seeds.txt")
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

## seedFiles=file.path(config.data.dir, "juelich_amygdala_seeds.txt")

seedFiles=file.path(config.data.dir, "juelich_amygdala_seeds_weights.txt")

for (seedFile in seedFiles) {

    seeds=readSeedsFile(seedFile)
    for (seed in seeds ) {

        seedName=getSeedName(seed)
        
        cat("####################################################################################################\n")
        cat(sprintf("*** Graphing ROIs for the %s seedName\n", toupper(seedName)))

        ## roiStats.fwhm4.2.mddAndCtrl.L_BLA.3mm.averageContrastValue.txt
        ## roiStats.fwhm4.2.mddAndCtrl.L_BLA.3mm.averageTValue.txt
        ## roiStats.fwhm4.2.mddAndCtrl.L_BLA.3mm.txt
        
        roistats.filename=file.path(group.results.dir, sprintf("roiStats.fwhm%s.%s.%s.txt", usedFwhm, groups, seedName))
        roistats.averageTvalue.filename=file.path(group.results.dir, sprintf("roiStats.fwhm%s.%s.%s.averageTValue.txt", usedFwhm, groups, seedName))
        roistats.averageContrastValue.filename=file.path(group.results.dir, sprintf("roiStats.fwhm%s.%s.%s.averageContrastValue.txt", usedFwhm, groups, seedName))
        
        if(file.exists(roistats.filename)) {
            
            roistats=readStatsTable(roistats.filename)
            roistats.averageTvalue=readStatsTable(roistats.averageTvalue.filename)
            roistats.averageContrastValue=readStatsTable(roistats.averageContrastValue.filename) 
            if (dim(roistats)[1] != expectedNumberOfSubjects) {
                stop(sprintf("Check the roistats variable. It has dim[%d, %d] but am expecting %d subjects\n",
                             dim(roistats)[1], dim(roistats)[2], expectedNumberOfSubjects))
            }
            
            clusterCount=length(grep("Mean", colnames(roistats)))
            if (clusterCount > 0 ) {
                cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
                
                clustersFilename=sprintf("clust.fwhm%s.%s.%s.txt", usedFwhm, groups, seedName)
                clusters=readClustersTable(clustersFilename)
                
                ## this table contains the locations, as text, of the clusters and is the output of a perl script
                clusterLocationsFilename=file.path(group.results.dir, sprintf("clusterLocations.fwhm%s.%s.%s.csv", usedFwhm, groups, seedName))
                clusterWhereAmI=readClusterLocationsTable(clusterLocationsFilename)
                
                ## this file stores the order of the subjects in each of the following BRIK files
                subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", groups, seedName, "csv", sep="."))
                subjectOrder=fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename))
                
                cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

                mgd=cbind(subjectOrder,
                    demographics[match(subjectOrder$subject, demographics$ID), selectedColumns]
                          )

                mgd$subjectWithTimePoint=as.factor(paste(mgd$subject, "A", sep="_"))
                mgd=cbind(
                    mgd,
                    sles [match(mgd$subjectWithTimePoint, sles$subject),
                          c("num.stressful.events", "num.severe.events", "total.sum.stress", "num.1.little.to.no.effect", "num.2.some.effect", "num.3.moderate.effect", "num.4.great.effect")])
                ## remove the now unnecessary subject with timepoint column
                mgd$subjectWithTimePoint=NULL

                mgd$MASC.tscore=as.numeric(mgd$MASC.tscore)
                
                mgd=droplevels(mgd) 

                mgd=fixDates(mgd)
                mgd=computeAge(mgd)

                mgd=computeMascScore(mgd)

                mgd=cbind(roistats, mgd)
                
                melted.roistats=melt(mgd, id.vars=c("subject", "Grp", "Gender", rvs),
                    measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                
                melted.roistats$cluster=factor(melted.roistats$cluster,
                    levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                    labels=paste(1:length(clusterWhereAmI), clusterWhereAmI))


                ## now drop the normal controls
                mgd=mgd[mgd$Grp=="MDD", ]

                ## runRegressionsForKajaWithGrpInModel(seed, melted.roistats, regressionVariables, file.path(group.results.dir, "regressionsForKaja.csv"), append=TRUE)
                ## runRegressionsForKajaWithGrpSeparated(seed, melted.roistats, regressionVariables, file.path(group.results.dir, "regressionsForKajaSeparatedByGrp.csv"), append=TRUE)
                

                runCorrelations(seed, melted.roistats, regressionVariables, clusters, inMethod="spearman",
                                file.path(group.results.dir, paste("correlations", seedName, "csv", sep=".")))

                ## runCorrelations(seed, melted.roistats, regressionVariables, clusters, inMethod="spearman", NULL)
                
                ## graphRegressions(seed, melted.roistats, regressionVariables, clusters, inMethod="spearman")
                
            } ## end of if (clusterCount > 0 )
        } else {
            cat("### ", roistats.filename, "does not exist. Skipping\n")
        } ## end of if(file.exists(roistats.filename))
    } ## end of for (seedName in seedNames)
} ## end of for (seedFile in seedFiles) {
