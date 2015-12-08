
fullatom = F
global.data.loaded=F
Allscores.loaded=F
global.filtered.data=F
global.regressions.done=F
global.regressions.done=T

######## score plotting

source("jplots.R")


i = 0
for (pdb in pdbs) {
  i = i+1
  # add centroid score
  all.data[[i]]$score4d = all.data[[i]]$d.env + all.data[[i]]$d.pair +
    pmax(all.data[[i]]$contact,-10) + all.data[[i]]$d.vdw +
      pmax(all.data[[i]]$d.fab)
}

global.score.list = append(global.score.list,"score4d")
global.score.list = list("d.env","d.pair","contact","d.vdw", "d.fab","score4d")

plot.scores(all.data)
#plot.targets(filtered.data)
