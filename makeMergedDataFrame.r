rm(list=ls())
graphics.off()

library(reshape)
library(nlme)
## library(dplyr)
library(multcomp)

source("scoreMasc.r")

########################################################################################################################################################################################################
### START OF FUNCTIONS #################################################################################################################################################################################
########################################################################################################################################################################################################

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c(
                                            "NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "#NAME?", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", "",
                                            "Test_Subject_", "Test_Subject", "NO SUBJECT", " NO SUBJECT", "NO SUBJECT ", " NO SUBJECT ", "Cancelled",
                                            "Head too big", "IGNORE THIS SUBJ", "NOT ELIGIBLE", "Permanent retainer", "PNS, fibromyalgia", "Psychosis/MDD"
    ))

    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
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

        cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f,\n", r, subjectNumber, gender, age, inData[r, "MASC.total"], old.masc.tscore))
        
        new.masc.tscore=scoreMasc(gender, age, inData[r, "MASC.total"])
        if (is.na(new.masc.tscore) ) {
            warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, inData[r, "MASC.total"]))
        }
        
        inData[r, "MASC.tscore"]=new.masc.tscore
        
        ## cat (sprintf("Old MASC tscore=%0.0f, new MASC tscore=%0.0f\n", old.masc.tscore, new.masc.tscore))
    }
    return (inData)
}

fixNonNumericColumns <- function (in.data, in.indices) {

    ## print (in.indices)
    for (ii in in.indices) {

        if ( is.factor(mgd[, ii]) ) {
            
            ## cat("*** The column being processed is named:", colnames(in.data)[ii], "\n")
        
            cn=grep ("-?[0-9]+", names(summary( mgd[, ii] )))
            ## cat("*** cn is", cn, "\n")
            sn=setdiff(names(summary( mgd[, ii] )) , names(summary( mgd[, ii] ))[cn])

            ## cat("*** Summary of column\n")
            ## print(summary( mgd[, ii] ))

            ## cat("*** Summary of non numeric factor levels\n")            
            ## print(summary( mgd[, ii] )[sn])
            
            ## now check those columns not matching NA to see if they have
            ## more than 0 entries in the column in mgd, if so we can't
            ## fix up that particular column in the data frame

            ## the "Other" is in the regexp below solely to capture
            ## problems with the C_EmoSuscept.Total column
            nonNumericLevelsWithMoreThanZeroEntries=summary( mgd[, ii] )[sn[grep("NA|Other", sn, invert=TRUE)]] > 0
            
            if ( (any ( nonNumericLevelsWithMoreThanZeroEntries) ) ) {
                cat("*** The following non-numeric entries in",  colnames(in.data)[ii] ,"have more than 0 elements. Cannot convert this column to a numeric type\n")
                cat ("** Summary of column", colnames(in.data)[ii], ":\n")
                print(summary( mgd[, ii] ))
                cat ("** Summary of factor level causing problems in", colnames(in.data)[ii], ":\n")                
                print(summary( mgd[, ii] )[sn[grep("NA", sn, invert=TRUE)]] ) 
            } else {
                cat("*** Fixing column", colnames(in.data)[ii], "\n")
                in.data[, ii] = as.numeric(as.character(in.data[, ii]))
            }
            ## break
        } else {
            cat("*** Skipping non factor column", colnames(in.data)[ii], "\n")
        }
    }

    return(in.data)
}

fixGroupMembership <- function() {

    ss=demo[, c("ID", "Grp")]

    validSubjectToGroupMapping=ss[apply(ss, 1, function (xx) { ! any(is.na(xx)) } ), ]

    ## cat("*** fixGroupMembership before\n")
    ## print(mgd[mgd$ID %in% c("343", "371", "422"), c("ID", "SubjID", "Grp", "BDI.II.Total", "CDI.Total")])

    mgd$Grp = demo[match(mgd$ID, validSubjectToGroupMapping$ID), "Grp"]

    ## cat("*** fixGroupMembership after\n")    
    ## print(mgd[mgd$ID %in% c("343", "371", "422"), c("ID", "SubjID", "Grp", "BDI.II.Total", "CDI.Total")])

    return (mgd)
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
    ncpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
    ncpus=8
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))    
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

scripts.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts")
data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
group.data.dir=file.path(data.dir, "Group.data")
group.results.dir=file.path(data.dir, "Group.results")
seeds.data.dir=file.path(data.dir, "seeds")

## define the filenames used to read in the CSV files containing all of the data
demo.filename   =file.path(admin.data.dir, "0-data_entry_current_2014withProperTimepointColumns.demo.csv")
papers.filename =file.path(admin.data.dir, "0-data_entry_current_2014withProperTimepointColumns.papers.csv")
nps.filename    =file.path(admin.data.dir, "0-data_entry_current_2014withProperTimepointColumns.nps.csv")
survey1.filename=file.path(admin.data.dir, "0-data_entry_current_2014withProperTimepointColumns.survey1.csv")
survey2.filename=file.path(admin.data.dir, "0-data_entry_current_2014withProperTimepointColumns.survey2.csv")
survey3.filename=file.path(admin.data.dir, "0-data_entry_current_2014withProperTimepointColumns.survey3.csv")
survey4.filename=file.path(admin.data.dir, "0-data_entry_current_2014withProperTimepointColumns.survey4.csv")
sles.filename   =file.path(admin.data.dir, "SLES_20140820.csv")
ctq.filename    =file.path(admin.data.dir, "Exisiting CTQ Scored_with RC_3_10_2015.csv")
wasiFilename    =file.path(admin.data.dir, "WASI.csv")


## now read in the CSV files and perfrom preliminary removal of
## unwanted rows and columns
demo=readCsvFile(demo.filename)
demo=demo[, 1:15]

papers=readCsvFile(papers.filename)
papers=papers[-c(1, 2), c(1, 2, 3, 6, 7, 8, 9)]

nps=readCsvFile(nps.filename, inSubjectColumnName="SubID")
nps=nps[-c(1, 2), c(1, 8:18)]

survey1=readCsvFile(survey1.filename)
survey1=survey1[-c(1, 2), 1:16]

survey2=readCsvFile(survey2.filename)
survey2=survey2[-c(1, 2), 1:11]
## eliminate Age and Gender columns from survey2
survey2=survey2[, -grep("Age|Gender", colnames(survey2))]

survey3=readCsvFile(survey3.filename)
survey3=survey3[-c(1, 2), 1:14]

survey4=readCsvFile(survey4.filename)
survey4=survey4[-c(1, 2), 1:9]

sles=readCsvFile(sles.filename, "subject")
## keep only columns that do NOT start with an X
sles=sles[, grep("^X", colnames(sles), fixed=FALSE, invert=TRUE)]
## delet ethe Gndr column, it will be provided by another data frame
sles=sles[, -grep("Gndr", colnames(sles))]
## Ensure that the columns are tagged with the name SLES so we know
## for sure what's in those columns
colnames(sles)=c(colnames(sles)[1:3], paste("SLES", colnames(sles)[4:dim(sles)[2]], sep=".") )
## rename the first two columns to be consistent with other data
## frames
colnames(sles)[1:2]=c("SubjID", "ID")

ctq=readCsvFile(ctq.filename)
## keep only columns that do NOT start with an X
ctq=ctq[, grep("^X", colnames(ctq), fixed=FALSE, invert=TRUE)]
## dump rows of all NAs
ctq=ctq[apply(ctq, 1, function(xx) { ! all(is.na(xx)) }), ]

wasi=readCsvFile(wasiFilename, inSubjectColumnName="SubID")
## tag the WASI columns with "WASI." so that is obvious what the
## columns mean, since they could get swamped among the rest of the
## columns
colnames(wasi)=c(colnames(wasi)[1], paste("WASI", colnames(wasi)[-1], sep="."))

## now perform the merging of all the data frames
## mgd1=merge(demo, papers,  by.x="ID",                           by.y="ID",                           all=TRUE)
## mgd2=merge(mgd1, nps,     by.x="ID",                           by.y="SubID" ,                       all=TRUE)
## mgd3=merge(mgd2, survey1, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
## mgd4=merge(mgd3, survey2, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
## mgd5=merge(mgd4, survey3, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
## mgd6=merge(mgd5, survey4, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
## mgd7=merge(mgd6, sles,    by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
## mgd= merge(mgd7, ctq, by.x="ID", by.y="ID", all=TRUE)


##mgd1=merge(demo,   papers,  by.x="ID",                           by.y="ID",                           all=TRUE)
##mgd1=merge(papers, nps,     by.x="ID",                           by.y="SubID" ,                       all=TRUE)

## replace teh WASI from the huge excel file that contaisn everyting
## with the correctly scored WASI data from a supplemental xl file
mgd1=merge(papers, wasi,    by.x="ID",                           by.y="SubID" ,                       all=TRUE)
mgd2=merge(mgd1,   survey1, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
mgd3=merge(mgd2,   survey2, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
mgd4=merge(mgd3,   survey3, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
mgd5=merge(mgd4,   survey4, by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
mgd6=merge(mgd5,   sles,    by.x=c("SubjID", "ID", "timepoint"), by.y=c("SubjID", "ID", "timepoint"), all=TRUE)
mgd= merge(mgd6,   ctq,     by.x="ID",                           by.y="ID",                           all=TRUE)

selectedColumns=c("Grp", "Gender", "DOB", "MRI", "Ethnicity", "Race", "Hand", "SES", "Tanner1", "Tanner2", "TannerAvg")
mgd=cbind(
    mgd,
    demo    [match(mgd$ID, demo$ID), selectedColumns]
)

## this just reorders the columns in the mgd data frame
mgd=mgd[,
    c(grep("ID|SubjID|timepoint|Grp|Gender|DOB|MRI|Ethnicity|Race|Hand|SES|Tanner1|Tanner2|TannerAvg", colnames(mgd), fixed=FALSE),
      grep("ID|SubjID|timepoint|Grp|Gender|DOB|MRI|Ethnicity|Race|Hand|SES|Tanner1|Tanner2|TannerAvg", colnames(mgd), fixed=FALSE, invert=TRUE))]

## these are redundant but useful for debugging purposes, comment the
## line below to keep them
## rm(list=c(paste("mgd", 1:7, sep="")))

## remove unwanted columns (age, gender and MASC.tscore will be added
## back in later
##mgd=mgd[, -c(grep ("Custom\\.SubjID|SubjID\\.Custom\\.Data|Custom\\.Data\\.SubjID|Age|Gender", colnames(mgd), fixed=FALSE))]
mgd=mgd[, -c(grep ("Custom\\.SubjID|SubjID\\.Custom\\.Data|Custom\\.Data\\.SubjID|Age|eth\\.code|race\\.code", colnames(mgd), fixed=FALSE))]

## mgd=fixGroupMembership()

## delete subjects 375 and 397 as these are to be ignored according to
## the demographics spreadsheet
## delete 305 and 327 as these have no MRI data for some reason
mgd=mgd[  ! (mgd$ID %in% c("305", "327", "375", "408") ), ]

mgd=rename(mgd, c("MASC.t.score"="MASC.tscore", "MASC.TOTAL"="MASC.total", "CGI.CGAS"="CGAS"))
## remove any row (subject) with a NA ID
mgd=mgd[ ! is.na(mgd$ID), ]
mgd$ID=droplevels(mgd$ID)

## fix the columns that should be numer ic but are not
mgd=fixNonNumericColumns(mgd, seq(12, dim(mgd)[2]))

## fix the subject with the extra space in their gender
mgd[which(mgd$Gender=="M "), "Gender"]="M"
mgd$Gender=droplevels(mgd$Gender)

## 378 identified as transgender, was biologically female and not on any hormone therapy
mgd[mgd$ID=="378", "Gender"]="F"
mgd$Gender=droplevels(mgd$Gender)

mgd=fixDates(mgd)
mgd=computeAge(mgd)

mgd=computeMascScore(mgd)

## delete rows with NA for timepoint
mgd=mgd[ ! is.na(mgd$timepoint), ]
mgd$timepoint=droplevels(mgd$timepoint)
## drop rows where group is NA
mgd=mgd[ ! is.na(mgd$Grp), ]
## drop rows where Grp is not MDD or NCL
mgd=mgd[mgd$Grp %in% c("MDD", "NCL"), ]
mgd$Grp=droplevels(mgd$Grp)
##remove rownames
rownames(mgd)=NULL

## sort the data frame by ID and then by timepoint
mgd=mgd[order(mgd$ID, mgd$timepoint), ]

## save the mgd file
date.tag=format(Sys.Date(), "%Y.%m.%d")
merged.file.name=file.path(admin.data.dir, sprintf("merged.demographics.and.neuropsych.%s.csv", date.tag))
cat("*** Wrinting merged demographics and neuropsychiactric measures to", merged.file.name, "\n")
cat("*** WARNING: the MASC t-scores for the timepoints other than A will not be accurate since the MRI data refers to only timepoint A\n")
write.csv(mgd, merged.file.name, row.names=FALSE, col.names=TRUE, quote=TRUE)
