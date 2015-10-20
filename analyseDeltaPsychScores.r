rm(list=ls())
graphics.off()

library(ggplot2)
library(reshape2)
library(nlme)
library(effects)
library(plyr)
library(multcomp)

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

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}

readSubjectOrderTable <- function (inFilename) {
    cat("*** Reading", subjectOrderFilename, "\n")
    subjectOrder=read.csv(subjectOrderFilename, header=T)

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

checkColumnAgainstThreshold <- function (inData, inColumns, inThreshold, inDirection=c("greater", "less"), replacement=NULL) {
    inDirection=match.arg(inDirection)

    for (column in inColumns) {
        ## print(inData[, c("subject", "Grp", column)])
    
        if (inDirection == "greater") {
            overThresholdRows=inData[, column] > inThreshold
            if (any(overThresholdRows, na.rm = TRUE)) {
                ## print(inData[, column] > inThreshold)
                cat ("****************************************************************************************************\n")
                cat (sprintf("*** The following subjects have %s data greater than %f\n", column, inThreshold))
                
                cat (paste (as.vector ( inData[which(overThresholdRows), "subject"]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(overThresholdRows), "Grp"    ]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(overThresholdRows), column   ]), collapse=" "), "\n")
                if( ! is.null(replacement)) {
                    cat("Replacing with", replacement, "\n")
                    inData=inData[which(overThresholdRows), column] = replacement
                }
                cat ("****************************************************************************************************\n")
            }
        } else if (inDirection == "less") {
            underThresholdRows=inData[, column] < inThreshold
            if (any(underThresholdRows, na.rm = TRUE)) {
                cat ("****************************************************************************************************\n")
                cat (sprintf("*** The following subjects have %s data less than %f\n", column, inThreshold))
                
                cat (paste (as.vector ( inData[which(underThresholdRows), "subject"]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(underThresholdRows), "Grp"    ]), collapse=" "), "\n")
                cat (paste (as.vector ( inData[which(underThresholdRows), column   ]), collapse=" "), "\n")
                if( ! is.null(replacement)) {
                    cat("Replacing with", replacement, "\n")
                    inData=inData[which(underThresholdRows), column] = replacement
                }
                cat ("****************************************************************************************************\n")                
            }
        }
    }
    return(inData)
}

createMergedDataFrames <- function(subjectOrder, demographics, change.data.frame, subject.id.column, change.column.names, value.name) {

    mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender")])

    mgd=droplevels(
        cbind(mgd,
              change.data.frame[match(subjectOrder$subject, change.data.frame[, subject.id.column]), change.column.names]))
    
    ## Rename the NCL level to HCL
    levels(mgd$Grp)[levels(mgd$Grp)=="NCL"]="HCL"
    mgd$Grp=relevel(mgd$Grp, ref="HCL")
    
    mgd.melted=melt(mgd, id.vars=c("subject", "Grp", "id", "Gender"), measure.vars=change.column.names, value.name = value.name)
    mgd.melted <- rename(mgd.melted, c("variable"="timepoint"))
    mgd.melted <- rename(mgd.melted, c("value"=value.name))    
    mgd.melted$timepoint = factor(mgd.melted$timepoint,
        levels=levels(mgd.melted$timepoint), labels=c("Baseline", "3 months follow-up"))

    return(list("wide"=mgd, "long"=mgd.melted))

}

createMergedDiffDataFrame <- function(subjectOrder, demographics, change.data.frame, subject.id.column, change.column.name) {

    mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender")])

    ## print(head(change.data.frame))
    
    mgd=droplevels(
        cbind(mgd,
              change.data.frame[match(subjectOrder$subject, change.data.frame[, subject.id.column]), change.column.name]))
    ## Rename the NCL level to HCL
    levels(mgd$Grp)[levels(mgd$Grp)=="NCL"]="HCL"
    mgd$Grp=relevel(mgd$Grp, ref="HCL")

    ##print(colnames(mgd))
    ## print(c(colnames(mgd)[-length(colnames(mgd))], change.column.name))
    colnames(mgd)=c(colnames(mgd)[-length(colnames(mgd))], change.column.name)


    return(mgd)
}

createInteractionGraph <- function (data, y.axis, ylab, title=NULL) {
    
    data.sumamry.df=summarySE(data, measurevar=y.axis, groupvars=c("Grp", "timepoint"), na.rm=TRUE)
    my.dodge=position_dodge(.2)
    interaction.graph=ggplot(data, aes_string(x="Grp", y=y.axis, group="timepoint", color="timepoint", shape="timepoint")) +
        geom_jitter() +
            geom_errorbar(data=data.sumamry.df, aes_string(ymin=paste(y.axis, "se", sep="-"), ymax=paste(y.axis, "se", sep="+")), width=.2, position=my.dodge) +
                stat_summary(fun.y = mean, geom="point", position=my.dodge) +
                    stat_summary(fun.y=mean, geom="line", position=my.dodge) +
                        scale_color_brewer(name="Timepoint:", palette="Set1") +
                            scale_shape_discrete(name="Timepoint:") +                    
                                labs(x="Group", y=ylab) +
                                    my.theme
    if(! is.null(title)) {
        interaction.graph = interaction.graph + ggtitle(title)
    }

    return(interaction.graph)
}

new.createInteractionGraph <- function (data, y.axis, ylab, title=NULL) {
    
    data.sumamry.df=summarySE(data, measurevar=y.axis, groupvars=c("Grp", "timepoint"), na.rm=TRUE)
    my.dodge=position_dodge(.2)
    interaction.graph=ggplot(data, aes_string(x="timepoint", y=y.axis, group="Grp", color="Grp", shape="Grp")) +
        geom_jitter() +
            geom_errorbar(data=data.sumamry.df, aes_string(ymin=paste(y.axis, "se", sep="-"), ymax=paste(y.axis, "se", sep="+")), width=.2, position=my.dodge) +
                stat_summary(fun.y = mean, geom="point", position=my.dodge) +
                    stat_summary(fun.y=mean, geom="line", position=my.dodge) +
                        scale_color_brewer(name="Group:", palette="Set1") +
                            scale_shape_discrete(name="Group:") +                    
                                labs(x="Timepoint", y=ylab) +
                                    my.theme
    if(! is.null(title)) {
        interaction.graph = interaction.graph + ggtitle(title)
    }

    return(interaction.graph)
}

analyse <- function (in.data, y.axis, na.action="na.omit") {

    cat("####################################################################################################\n")
    cat("### Analysis of", y.axis, "\n")
    model.formula=as.formula(paste(y.axis, "~", "Grp * timepoint"))
    ## random.formula=as.formula("~ 1 + timepoint | subject")
    random.formula=as.formula("~ 1 | subject")    

    in.data=in.data[complete.cases(in.data), ]
    cat("### Table of n\n")
    print(with(in.data, addmargins(table(Grp, timepoint))))
    
    model=lme(data=in.data, fixed=model.formula, random=random.formula)
    model.anova=anova(model)

    print(summarySE(in.data, measure=y.axis, groupvars=c("Grp", "timepoint")))

    ss=subset(in.data, Grp=="MDD")
    ## print(ss)
    ## print(ss[ss$timepoint=="Baseline", ])
    ## print(ss[ss$timepoint=="3 months follow-up", ])
    
    mdd.only.ttest=t.test(ss[ss$timepoint=="Baseline", y.axis], ss[ss$timepoint=="3 months follow-up", y.axis], paired=TRUE)

    cat("\n\nMDD ONLY T TEST\n")
    print(mdd.only.ttest)
    cat("\n\n###LME:\n")
    print(model)
    cat("\n\n###ANOVA:\n")    
    print(model.anova)

    cat("\n\n###Effects:\n")    
    ## see
    ## http://cran.r-project.org/web/packages/car/vignettes/embedding.pdf
    ## for the source of this solution to the the scoping rules of R
    ## causing problems with not being able to find variables
    assign("model.formula", model.formula, envir=.GlobalEnv)
    assign("in.data", in.data, envir=.GlobalEnv)
    print(Effect("Grp", model))
    print(Effect("timepoint", model))
    print(allEffects(model, formula=model.formula))
    remove("model.formula", envir=.GlobalEnv)
    remove("in.data", envir=.GlobalEnv)            
    
}


new.analyse <- function (in.data, y.axis, na.action="na.omit") {

    cat("####################################################################################################\n")
    cat("### Analysis of", y.axis, "\n")

    cat("### Summary of", y.axis, "\n")
    print(summarySE(in.data, measure=y.axis, groupvars=c("Grp")))
    print(summarySE(in.data, measure=y.axis, groupvars=c("timepoint")))
    smry=summarySE(in.data, measure=y.axis, groupvars=c("Grp", "timepoint"))
    print(smry)
    print(make.table.strings(smry, y.axis))

    cat("### Table of n\n")
    print(with(in.data, addmargins(table(Grp, timepoint))))
    
    ## setup the froomulae used in the models
    model.formula.nointeraction=as.formula(paste(y.axis, "~", "Grp + timepoint"))
    model.formula.interaction=as.formula(paste(y.axis, "~", "Grp * timepoint"))    
    ## random.formula=as.formula("~ 1 + timepoint / subject")
    random.formula=as.formula("~ 1 | subject")    

    in.data=in.data[complete.cases(in.data), ]

    ## setup the LME models and ANOVAs derived therefrom
    model.lme.interaction=lme(data=in.data, fixed=model.formula.interaction, random=random.formula)
    ## print(model.lme.interaction)
    model.lme.nointeraction=lme(data=in.data, fixed=model.formula.nointeraction, random=random.formula)
    ## print(model.lme.nointeraction)
    model.anova.interaction=anova(model.lme.interaction)
    model.anova.nointeraction=anova(model.lme.nointeraction)

    cat("\n### Interaction LME summary\n")
    print(summary(model.lme.interaction))
    cat("\n### Nointeraction LME summary\n")
    print(summary(model.lme.nointeraction))    

    cat("\n### Interaction ANOVA summary\n")
    print(model.anova.interaction)
    cat("\n### Nointeraction ANOVA summary\n")
    print(model.anova.nointeraction)

    use.interaction.model=FALSE
    if (model.anova.interaction[4, "p-value"] < 0.05) {
        use.interaction.model=TRUE
        cat("\n")
        cat("**************************************************************\n")
        cat("*** The Grp * timepoint is significant: using interation model\n")
        cat("**************************************************************\n")
    }

    
    ## ## original from
    ## ## http://cran.r-project.org/web/packages/multcomp/vignettes/multcomp-examples.pdf
    ## ## mod <- lm(breaks ~ wool * tension, data = warpbreaks)
    ## ## tmp <- expand.grid(tension = unique(warpbreaks$tension), wool = unique(warpbreaks$wool))
    ## ## X <- model.matrix(~ wool * tension, data = tmp)
    ## ## glht(mod, linfct = X)
    ## ## Tukey <- contrMat(table(warpbreaks$tension), "Tukey")
    ## ## K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
    ## ## rownames(K1) <- paste(levels(warpbreaks$wool)[1], rownames(K1), sep = ":")
    ## ## K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
    ## ## rownames(K2) <- paste(levels(warpbreaks$wool)[2], rownames(K2), sep = ":")
    ## ## K <- rbind(K1, K2)
    ## ## colnames(K) <- c(colnames(Tukey), colnames(Tukey))
    ## ## ## and perform the tests via
    ## ## summary(glht(mod, linfct = K %*% X))
    
    if (use.interaction.model==FALSE) {
        cat("### Main effects in nointeraction model\n")
        assign("model.formula.nointeraction", model.formula.nointeraction, envir=.GlobalEnv)
        K1 <- glht(model.lme.nointeraction, mcp(Grp = "Tukey"))$linfct
        K2 <- glht(model.lme.nointeraction, mcp(timepoint = "Tukey"))$linfct
        print(summary(glht(model.lme.nointeraction, linfct = rbind(K1, K2))))
        remove("model.formula.nointeraction", envir=.GlobalEnv)        
    }

    if (use.interaction.model==TRUE) {
        cat("### Interaction effects\n")
        tmp <- expand.grid(timepoint = unique(in.data$timepoint), Grp = unique(in.data$Grp))
        X.a <- model.matrix(~ Grp * timepoint, data = tmp)
        ## print(glht(model.lme, linfct = X.a))
        Tukey.a <- contrMat(table(in.data$timepoint), "Tukey")
        Ka1 <- cbind(Tukey.a, matrix(0, nrow = nrow(Tukey.a), ncol = ncol(Tukey.a)))
        rownames(Ka1) <- paste(levels(in.data$Grp)[1], rownames(Ka1), sep = ":")
        Ka2 <- cbind(matrix(0, nrow = nrow(Tukey.a), ncol = ncol(Tukey.a)), Tukey.a)
        rownames(Ka2) <- paste(levels(in.data$Grp)[2], rownames(Ka2), sep = ":")
        K.a <- rbind(Ka1, Ka2)
        colnames(K.a) <- c(colnames(Tukey.a), colnames(Tukey.a))

        print(summary(glht(model.lme.interaction, linfct = K.a %*% X.a)))


        ##tmp <- expand.grid(Grp = unique(in.data$Grp), timepoint = unique(in.data$timepoint))
        X.b <- model.matrix(~ timepoint * Grp, data = tmp)
        ## print(glht(model.lme, linfct = X.b))
        Tukey.b <- contrMat(table(in.data$Grp), "Tukey")
        Kb1 <- cbind(Tukey.b, matrix(0, nrow = nrow(Tukey.b), ncol = ncol(Tukey.b)))
        rownames(Kb1) <- paste(levels(in.data$timepoint)[1], rownames(Kb1), sep = ":")
        Kb2 <- cbind(matrix(0, nrow = nrow(Tukey.b), ncol = ncol(Tukey.b)), Tukey.b)
        rownames(Kb2) <- paste(levels(in.data$timepoint)[2], rownames(Kb2), sep = ":")
        K.b <- rbind(Kb1, Kb2)
        colnames(K.b) <- c(colnames(Tukey.b), colnames(Tukey.b))
        ## and perform the tests via
        print(summary(glht(model.lme.interaction, linfct = K.b %*% X.b)))

        assign("in.data", in.data, envir=.GlobalEnv)
        assign("model.formula.interaction", model.formula.interaction, envir=.GlobalEnv)
        
        cat("*** Tables of means\n")
        print(allEffects(model.lme.interaction, formula=model.formula.interaction))

        remove("model.formula.interaction", envir=.GlobalEnv)
        remove("in.data", envir=.GlobalEnv)
    }
}

make.line.entry <- function (in.data, in.variable, in.group, in.timepoint) {

    nacount=in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "nacount"]
    nacount.string=""
    if(nacount > 0) {
        nacount.string=sprintf(" [%d]", nacount)
    }
    
    line.entry=sprintf("%s, %s: \"%0.1f ± %0.1f (%0.1f, %0.1f) %s\"",
        in.group, in.timepoint, 
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, in.variable],
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "se"],
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "min"],
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "max"],
        nacount.string)
    
    return(line.entry)
}

make.table.strings <- function (in.data, in.variable) {

    group.by.timepoint=expand.grid(unique(in.data$Grp), unique(in.data$timepoint))

    lines=apply(group.by.timepoint, 1,
        function(xx) {
            make.line.entry(in.data, in.variable, xx[1], xx[2])
        })

    return(lines)
}

########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

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
## group.data.dir=file.path(data.dir, "Group.data")
group.results.dir=file.path(data.dir, "Group.results")
seeds.data.dir=file.path(data.dir, "seeds")

demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)

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


my.base.size=14
my.theme=
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

## sink("delta.psych.scores.with.simple.random.effects.txt")
####################################################
## Change in CDRSR scores from time A to C by group
## group.data.dir=file.path(data.dir, "Group.data.CDRSR.diff.withAandC")
group.data.dir=file.path(data.dir, "Group.data.CDRS.t.score.scaled.diff.withAandC")
subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", "mddAndCtrl", "L_BLA.weight.3mm", "csv", sep="."))
subjectOrder=splitSubjectOrderIntoIdAndTimepoint(fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename)))

dfs=createMergedDataFrames(subjectOrder, demographics, cdrsr.change, "SubjID", c("CDRSR.tscore.A", "CDRSR.tscore.C"), "CDRSR")
mgd=dfs[["wide"]]
mgd.melted=dfs[["long"]]

interaction.graph=new.createInteractionGraph(mgd.melted, y.axis="CDRSR", ylab="CDRSR T score", title="Change in CDRSR Score from\nbaseline to 3 months follow-up")
print(interaction.graph)
new.analyse(mgd.melted, y.axis="CDRSR", na.action="na.omit")
stop()
####################################################
## Change in CGAS scores from time A to C by group
group.data.dir=file.path(data.dir, "Group.data.CGAS.diff.withAandC")
subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", "mddAndCtrl", "L_BLA.weight.3mm", "csv", sep="."))
subjectOrder=splitSubjectOrderIntoIdAndTimepoint(fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename)))

mgd.diff=createMergedDiffDataFrame(subjectOrder, demographics, cgas.change, "SubjNum", "CGAS.diff")

## drop subject 111 how has a 32 point change in CGAS from baseline to
## 3 months
## mgd.diff=mgd.diff[mgd.diff$id != "111_A", ]

cgas.diff.summary=summarySE(mgd.diff, measure="CGAS.diff", groupvars=c("Grp"))
cgas.diff.graph=ggplot(cgas.diff.summary, aes(x=Grp, y=CGAS.diff, color=Grp)) +
    geom_jitter(data=mgd.diff) +
        geom_errorbar(aes(ymin=CGAS.diff-se, ymax=CGAS.diff+se), width=0.5) +
            stat_summary(data=mgd.diff, fun.y="mean", geom="point") +
            scale_color_brewer(name="Group:", palette="Set1") +  
                labs(title="Change in CGAS from\nbaseline to 3 months followup", x="Group", y=expression(paste(Delta, " CGAS"))) +
                    my.theme
dev.new(); print(cgas.diff.graph)
ggsave(file.path(data.dir, "Group.results.CGAS.diff.withAandC", "delta.cgas.baseline.3months.pdf"), cgas.diff.graph)

dfs=createMergedDataFrames(subjectOrder, demographics, cgas.change, "SubjNum", c("CGAS.A", "CGAS.C"), "CGAS")

mgd=dfs[["wide"]]
## drop subject 111 how has a 32 point change in CGAS from baseline to
## 3 months
## mgd=mgd[mgd$id != "111_A", ]

mgd.melted=dfs[["long"]]
## drop subject 111 how has a 32 point change in CGAS from baseline to
## 3 months
## mgd.melted=mgd.melted[mgd.melted$id != "111_A", ]


interaction.graph=new.createInteractionGraph(mgd.melted, y.axis="CGAS", ylab="CGAS", title="Change in CGAS Score from\nbaseline to 3 months follow-up")
dev.new(); print(interaction.graph)
new.analyse(mgd.melted, y.axis="CGAS", na.action="na.omit")
cat("\n*c** Summary of CGAS.diff\n")
print(cgas.diff.summary)

####################################################
## Change in RADS scores from time A to C by group
group.data.dir=file.path(data.dir, "Group.data.RADS.Total.Tscore.diff.withAandC")
subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", "mddAndCtrl", "L_BLA.weight.3mm", "csv", sep="."))
subjectOrder=splitSubjectOrderIntoIdAndTimepoint(fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename)))

dfs=createMergedDataFrames(subjectOrder, demographics, rads.change, "SubjNum", c("RADS.Total.Tscore.A", "RADS.Total.Tscore.C"), "RADS")
mgd=dfs[["wide"]]
mgd.melted=dfs[["long"]]

interaction.graph=new.createInteractionGraph(mgd.melted, y.axis="RADS", ylab="RADS", title="Change in RADS Score from\nbaseline to 3 months follow-up")
dev.new(); print(interaction.graph)
new.analyse(mgd.melted, y.axis="RADS", na.action="na.omit")

## ####################################################
## ## Change in MASC scores from time A to C by group
## group.data.dir=file.path(data.dir, "Group.data.MASC.tscore.diff.withAandC")
## subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", "mddAndCtrl", "L_BLA.3mm", "csv", sep="."))
## subjectOrder=splitSubjectOrderIntoIdAndTimepoint(fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename)))

## dfs=createMergedDataFrames(subjectOrder, demographics, masc.change, "SubjNum", c("MASC.tscore.A", "MASC.tscore.C"), "MASC")
## mgd=dfs[["wide"]]
## mgd.melted=dfs[["long"]]

## interaction.graph=createInteractionGraph(mgd.melted, y.axis="MASC", ylab="MASC T score", title="Change in MASC Score from\nbaseline to 3 months follow-up")
## dev.new(); print(interaction.graph)
## analyse(mgd.melted, y.axis="MASC", na.action="na.omit")



####################################################
## Change in BDI scores from time A to C by group
group.data.dir=file.path(data.dir, "Group.data.BDI.diff.withAandC")
subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", "mddAndCtrl", "L_BLA.weight.3mm", "csv", sep="."))
subjectOrder=splitSubjectOrderIntoIdAndTimepoint(fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename)))

dfs=createMergedDataFrames(subjectOrder, demographics, bdi.change, "SubjNum", c("BDI.A", "BDI.C"), "BDI")
mgd=dfs[["wide"]]
mgd.melted=dfs[["long"]]

interaction.graph=new.createInteractionGraph(mgd.melted, y.axis="BDI", ylab="BDI-II score", title="Change in BDI-II Score from\nbaseline to 3 months follow-up")
dev.new(); print(interaction.graph)
new.analyse(mgd.melted, y.axis="BDI", na.action="na.omit")


## ####################################################
## ## Change in CDI scores from time A to C by group
## group.data.dir=file.path(data.dir, "Group.data.CDI.diff.withAandC")
## subjectOrderFilename=file.path(group.data.dir, paste("subjectOrder", "mddAndCtrl", "L_BLA.3mm", "csv", sep="."))
## subjectOrder=splitSubjectOrderIntoIdAndTimepoint(fixSubjectOrderTable(readSubjectOrderTable(subjectOrderFilename)))

## dfs=createMergedDataFrames(subjectOrder, demographics, cdi.change, "SubjNum", c("CDI.A", "CDI.C"), "CDI")
## mgd=dfs[["wide"]]
## mgd.melted=dfs[["long"]]

## interaction.graph=createInteractionGraph(mgd.melted, y.axis="CDI", ylab="CDI-II score", title="Change in CDI-II Score from\nbaseline to 3 months follow-up")
## dev.new(); print(interaction.graph)
## analyse(mgd.melted, y.axis="CDI", na.action="na.omit")


##sink()
