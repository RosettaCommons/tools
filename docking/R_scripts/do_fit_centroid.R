
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

print("performing regression")
myglm = glm(good ~ d.env + d.pair + contact + d.vdw
            + env + pair + hs + ss + sheet + cb + rsigma + hb + rg + id,
            data = mat, family="binomial")

process.glm(myglm,"fit.centroid")

