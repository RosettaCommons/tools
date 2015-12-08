
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


# good decoys should not be bumpy
bumpcut = array(0,length(pdbs))
names(bumpcut)=pdbs
for(j in 1:length(pdbs)) {
  mean.rep = mean(mat$bk.rep[mat$id==pdbs[j] &
                             mat$description=="output_decoy"])
  bumpcut[j] = mean.rep
}
mat$good = mat$good & mat$bk.rep < bumpcut[mat$id]



print("performing regression")
glmFA =  glm(good ~ 
             bk.atr + bk.rep + bk.sol + bk.hbbb + bk.hbsc + bk.dun + bk.pair 
             + d.Eatr + d.Erep + d.Eatr.lr + d.Erep.lr + gsolt + id,
             data = mat,family="binomial")

process.glm(glmFA,"fit.bumpcut")
