rm(list=ls())
graphics.off()

library(MASS)
library(boot)

set.seed(1234567)

model=read.table(textConnection("
Grp subject age.in.years CDRS.t.score.A CDRS.t.score.C          mri
MDD     106     14.83242             80             62  0.028073812
MDD     111     14.17308             74             67 -0.048749052
MDD     112     17.27473             63             59  0.068577386
MDD     118     16.76648             53             54  0.047104210
MDD     136     16.44780             64             53  0.038103342
MDD     144     17.57692             63             53  0.310577482
MDD     147     16.88462             65             62 -0.009386584
MDD     158     16.80495             75             57  0.054385971
MDD     160     15.75549             63             52 -0.133139238
MDD     167     15.89560             57             57  0.043042418
MDD     301     14.91484             56             56  0.025857268
MDD     304     17.39011             82             69  0.100147866
MDD     313     16.79670             71             77  0.140405267
MDD     330     13.13462             75             77 -0.187526122
MDD     331     15.31044             69             81  0.193774804
MDD     335     16.64011             63             61 -0.195805609
MDD     336     14.15659             68             55  0.000000000
MDD     339     16.75824             81             63  0.036458120
MDD     351     17.85989             84             75  0.146064594
MDD     353     15.81319             79             84  0.193176582
MDD     364     15.25824             73             86 -0.032219738
MDD     366     13.93132             63             67  0.332895458
MDD     389     16.95055             75             58 -0.058629736"), header=TRUE)

model.formula=as.formula("CDRS.t.score.C ~ mri + age.in.years")


boot.huber <- function (data, indices, in.maxit=20) {

    df=data[indices, ]
    hb=rlm(model.formula, data=df, maxit=in.maxit)
    return(coefficients(hb))
}

voxel.boot = boot(data=model, statistic=boot.huber, R=100, in.maxit=1000)
print(voxel.boot)
    
bootstrappingBrikLabelSuffxes=c("booted.mean.beta", "bias", "booted.t.value", "ciLower", "ciUpper")

my.rlm=rlm(model.formula, data=model, maxit=20)
inRlmCoef=coefficients(summary(my.rlm))

nrows=length(rownames(inRlmCoef))
ncols=length(colnames(inRlmCoef))

indices=seq(ncols, nrows*ncols, by=ncols)

indices=c(indices, seq(indices[length(indices)] + grep("t.value", bootstrappingBrikLabelSuffxes), (nrows*ncols) + (nrows*length(bootstrappingBrikLabelSuffxes)), by=length(bootstrappingBrikLabelSuffxes)))


print (indices)

