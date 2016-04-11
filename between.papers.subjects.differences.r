rm(list=ls())

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}


emm.hcl.list=data.frame("subject" = c())
emm.mdd.list=data.frame("subject" = c())

sgacc.paper.group.data.dir=file.path(root.dir, "sanDiego/restingstate/data/Group.data")
amygdala.paper.group.data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.data")

sgacc.hcl.filename=file.path(sgacc.paper.group.data.dir, "subjectOrder.ctrlOnly.acci8L.REML.csv")
sgacc.mdd.filename=file.path(sgacc.paper.group.data.dir, "subjectOrder.mddOnly.acci8L.REML.csv")

amygdala.hcl.filename=file.path(amygdala.paper.group.data.dir, "subjectOrder.ctrlOnly.R_whole_amygdala.3mm.csv")
amygdala.mdd.filename=file.path(amygdala.paper.group.data.dir, "subjectOrder.mddOnly.R_whole_amygdala.3mm.csv")


sgacc.hcl.list=read.csv(sgacc.hcl.filename, header=TRUE)
sgacc.mdd.list=read.csv(sgacc.mdd.filename, header=TRUE)

amygdala.hcl.list=read.csv(amygdala.hcl.filename, header=TRUE)
amygdala.mdd.list=read.csv(amygdala.mdd.filename, header=TRUE)


emm.hcl.list=data.frame("subject" =
                            c("119_A", "121_A", "122_A", "123_A", "126_A", "127_A", "131_A", "135_A", "138_A", "139_A", "141_A", "142_A", "143_A",
                              "145_A", "146_A", "148_A", "151_A", "152_A", "155_A", "157_A", "159_A", "162_A", "163_A", "165_A", "306_A", "312_A",
                              "321_A", "326_A", "328_A", "337_A", "346_A", "348_A", "354_A", "357_A", "367_A", "369_A", "377_A"))
emm.mdd.list=data.frame("subject" =
                            c("134_A", "144_A", "147_A", "160_A", "161_A", "300_A", "304_A", "316_A", "317_A", "320_A", "330_A", "332_A", "336_A",
                              "339_A", "356_A", "358_A", "359_A", "360_A", "361_A", "362_A", "137_A", "310_A", "325_A", "372_A", "373_A", "376_A"))                      


cat("####################################################################################################\n")
cat("### sgACC and amygdala paper\n")

cat("*** Number of HCL subjects common to sgacc and amygdala RSFC paper: ")
print(length(intersect(sgacc.hcl.list$subject, amygdala.hcl.list$subject)))

cat("*** Number of MDD subjects common to sgacc and amygdala RSFC paper: ")
print(length(intersect(sgacc.mdd.list$subject, amygdala.mdd.list$subject)))

cat("*** Number of HCL in the amygdala paper but not in the sgacc paper: ")
print(length(setdiff(amygdala.hcl.list$subject, sgacc.hcl.list$subject)))

cat("*** Number of MDD in the amygdala paper but not in the sgacc paper: ")
print(length(setdiff(amygdala.mdd.list$subject, sgacc.mdd.list$subject)))


cat("*** Number of HCL in the sgacc paper but not in the amygdala paper: ")
print(length(setdiff(sgacc.hcl.list$subject, amygdala.hcl.list$subject)))

cat("*** Number of MDD in the sgacc paper but not in the amygdala paper: ")
print(length(setdiff(sgacc.mdd.list$subject, amygdala.mdd.list$subject)))


cat("####################################################################################################\n")
cat("### EMM and amygdala paper\n")

cat("*** Number of HCL subjects common to EMM and amygdala RSFC paper: ")
print(length(intersect(emm.hcl.list$subject, amygdala.hcl.list$subject)))

cat("*** Number of MDD subjects common to EMM and amygdala RSFC paper: ")
print(length(intersect(emm.mdd.list$subject, amygdala.mdd.list$subject)))

cat("*** Number of HCL in the amygdala paper but not in the EMM paper: ")
print(length(setdiff(amygdala.hcl.list$subject, emm.hcl.list$subject)))

cat("*** Number of MDD in the amygdala paper but not in the EMM paper: ")
print(length(setdiff(amygdala.mdd.list$subject, emm.mdd.list$subject)))


cat("*** Number of HCL in the EMM paper but not in the amygdala paper: ")
print(length(setdiff(emm.hcl.list$subject, amygdala.hcl.list$subject)))

cat("*** Number of MDD in the EMM paper but not in the amygdala paper: ")
print(length(setdiff(emm.mdd.list$subject, amygdala.mdd.list$subject)))
