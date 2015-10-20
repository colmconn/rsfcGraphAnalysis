rm(list=ls())
graphics.off()



## cd /data/sanDiego

## subjects=$( ls -d [0-9][0-9][0-9]_[AC]* ) 

subjects=dir(path="/data/sanDiego", pattern="^[0-9]{3,3}_[AC][0-9]?.*", include.dirs=TRUE, no..=TRUE)

##subjects=head(subjects, 10)

## cleanedSubjects=""
## for ss in $subjects ; do
##     cc=$( echo ${ss} | sed -e "s/\([0-9]\{3\}_[AC][0-9]\?\).*/\1/" )
##     cleanedSubjects="$cleanedSubjects $cc"

## done

cleanedSubjects=sapply(subjects, function(xx) {
    return(sub("([0-9]{3}_[AC][0-9]?).*", "\\1", xx))
}, USE.NAMES = FALSE)

## echo $cleanedSubjects | tr " " "\n" | uniq |  tr "\n" " "
## echo

cleanedSubjects=unique(sort(cleanedSubjects))
cat(sprintf("*** List of %d subjects with folders (reconstructed or not) or tarballs or all of these\n", length(cleanedSubjects)))
cat(cleanedSubjects, "\n")

bothTimePoints=vector(mode="character", length=length(cleanedSubjects))
nn=1
for (subject in sub("_[AC]", "", cleanedSubjects, fixed=FALSE) ) {
##    cat(subject, "\n")
    if ( paste(subject, "A", sep="_") %in% cleanedSubjects &&
         paste(subject, "C", sep="_") %in% cleanedSubjects) {
        bothTimePoints[nn]=subject
        nn=nn+1
    }
}
bothTimePoints=unique(bothTimePoints)
cat(sprintf("*** List of %d subjects with both A and C time points\n", length(bothTimePoints)))
cat(bothTimePoints, "\n")
