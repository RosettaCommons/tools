
fullatom = T
global.data.loaded=F
Allscores.loaded=F
global.filtered.data=F
global.regressions.done=F
global.regressions.done=T

stop("exiting init.R normally -- proceed by hand",call.=F)

######## score plotting

source("jplots.R")
plot.scores(filtered.data)
stop()

######## regressions

source("jfit=coef.R")

glmbk = make.glmbk()

summary(glmbk)
a=summary(glmbk)$coefficients
b=summary(glmbk)$coefficients/glmbk$coefficients[3]
write.table(round(a,5),file="fit.bk.coefficients",quote=F,sep='\t')
write.table(round(b,5),file="fit.bk.coefficients.rep",quote=F,sep='\t')
 
glmbk.enrichments=round(calc.multi.enrichments(glmbk),4)
write.table(glmbk.enrichments,file="fit.bk.enrichments",quote=F,sep='\t')

plot.fit(glmbk,"fit.bk.ps")

stop()

################ old stuff
summary(glmbb.d)
summary(glmbb)
summary(glmbk)
summary(glmfa)

plot.fit(glmbb.d,"fit.bb.d.ps")
plot.fit(glmbb,"fit.bb.ps")
plot.fit(glmfa,"fit.fa.ps")
plot.fit(glmbk,"fit.bk.ps")

