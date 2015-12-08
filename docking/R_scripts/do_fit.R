
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

glmbk =   glm(good ~ 
              bk.atr + bk.rep + bk.sol + bk.hbbb + bk.hbsc + bk.dun + bk.pair 
              + gsolt + id, 
              data = mat,family="binomial")

process.glm(glmbk,"fit.bk")
