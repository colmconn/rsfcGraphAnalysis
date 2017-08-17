rm(list=ls())
graphics.off()



## cd /data/sanDiego

## subjects=$( ls -d [0-9][0-9][0-9]_[ACD]* ) 

subjects=dir(path="/data/sanDiego", pattern="^[0-9]{3,3}_[ACD][0-9]?.*", include.dirs=TRUE, no..=TRUE)

##subjects=head(subjects, 10)

## cleanedSubjects=""
## for ss in $subjects ; do
##     cc=$( echo ${ss} | sed -e "s/\([0-9]\{3\}_[ACD][0-9]\?\).*/\1/" )
##     cleanedSubjects="$cleanedSubjects $cc"

## done

cleanedSubjects=sapply(subjects, function(xx) {
    return(sub("([0-9]{3}_[ACD][0-9]?).*", "\\1", xx))
}, USE.NAMES = FALSE)

subjectNumbers=unique(sort(sapply(cleanedSubjects, function(xx) {
    return(sub("([0-9]{3})(_[ACD][0-9]?)", "\\1", xx))
}, USE.NAMES = FALSE)))

## echo $cleanedSubjects | tr " " "\n" | uniq |  tr "\n" " "
## echo

cleanedSubjects=unique(sort(cleanedSubjects))
cat(sprintf("*** List of %d subjects with folders (reconstructed or not) or tarballs or all of these\n", length(cleanedSubjects)))
cat(cleanedSubjects, "\n")

bothTimePoints  = vector(mode="character", length=length(cleanedSubjects))
threeTimePoints = vector(mode="character", length=length(cleanedSubjects))
bb=1

for (subject in subjectNumbers ) {
    cat(subject, " ")
    if ( paste(subject, "A", sep="_") %in% cleanedSubjects &&
         paste(subject, "C", sep="_") %in% cleanedSubjects) {
        bothTimePoints[bb]=subject
        bb=bb+1
    }
}
cat("\n")
tt=1
for (subject in  subjectNumbers ) {
    cat(subject, " ")
    if ( paste(subject, "A", sep="_") %in% cleanedSubjects &&
         paste(subject, "C", sep="_") %in% cleanedSubjects &&         
         paste(subject, "D", sep="_") %in% cleanedSubjects) {
        threeTimePoints[tt]=subject
        tt=tt+1
    }
}

bothTimePoints=unique(bothTimePoints)
cat(sprintf("*** List of %d subjects with both A and C time points\n", length(bothTimePoints)))
cat(bothTimePoints, "\n")
threeTimePoints=unique(threeTimePoints)
cat(sprintf("*** List of %d subjects with A, C, and D time points\n", length(threeTimePoints)))
cat(threeTimePoints, "\n")
