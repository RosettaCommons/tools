
fullatom = T
global.data.loaded=F
Allscores.loaded=F
global.filtered.data=F
global.regressions.done=F
global.regressions.done=T

######## regressions

source("jfit=coef.R") # fitting
source("jE.R")        # enrichments
source("jz.R")        # zscores

glmbkelec = make.glmbkelec()

process.glm(glmbkelec,"fit.bkelec")
