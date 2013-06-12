
fullatom = T
global.data.loaded=F
Allscores.loaded=F
global.filtered.data=F
global.regressions.done=F
global.regressions.done=T

######## score plotting

source("jplots.R")

global.score.list = list("score")
plot.scores(filtered.data)
#plot.targets(filtered.data)
