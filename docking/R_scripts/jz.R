
# JJG Feb 02

# the z-score: average decoy score minus native score, normalized by the
# standard deviation of the decoys
# note: we take the mean of the "native.score" in case a scorefile
# has multiple natives
zscore = function(native.score,decoy.scores) {
  decoy.avg = mean(decoy.scores)
  decoy.sd  = sd(decoy.scores)

  ans = ((decoy.avg - mean(native.score))/decoy.sd)
  if (!is.finite(ans)) { ans = NA }

  return(ans)
}


calc.zscores = function(dat=mat,scores=-fitted(glmbk),nlabel="native"){
                      #(dat=mat,scores=mat$score){
  z=array(0,length(pdbs))
  names(z)=pdbs
  for (pdb in pdbs) {
    n = scores[dat$id==pdb & dat$description==nlabel]
    d = scores[dat$id==pdb & dat$description=="output_decoy"]
    z[pdb] = zscore(n,d)
  }
  zdat=z
  z["avg"]  = mean(zdat,na.rm=T)
  z["bad"]  = sum(zdat<2) # number of low z-scores
  z["fail"] = sum(zdat<0) # number negative z-scores
  z["n"]    = length(zdat)-sum(is.na(zdat))
  return(z)
}




# the "winning native" scores better than *all* the decoys (returns T or F)
# (take the mean of the native.score in case there are multiple natives)
winning.native = function(native.score,decoy.scores) {
  as.logical(mean(native.score) < min(decoy.scores))
}




calc.winning.natives = function(dat=mat,scores=-fitted(glmbk),nlabel="native"){
                              #(dat=mat,scores=mat$score){
  w=array(0,length(pdbs))
  names(w)=pdbs
  for (pdb in pdbs) {
    n = scores[dat$id==pdb & dat$description==nlabel]
    d = scores[dat$id==pdb & dat$description=="output_decoy"]
    w[pdb] = winning.native(n,d)
  }
  wdat=w
  w["avg"]=mean(wdat,na.rm=T)
  w["bad"]=length(wdat)-sum(wdat,na.rm=T) # number of losers
  w["fail"]=w["bad"]                      # ditto
  w["n"]=length(wdat)-sum(is.na(wdat))
  return(w)
}



make.fit.summary = function(dat=mat,glm=glmbk){
  e = calc.multi.enrichments(glm)
  z.nat = calc.zscores(dat,-fitted(glm))
  w.nat = calc.winning.natives(dat,-fitted(glm))
  z.rep = calc.zscores(dat,-fitted(glm),"inp_repacked")
  w.rep = calc.winning.natives(dat,-fitted(glm),"inp_repacked")
  z.min = calc.zscores(dat,-fitted(glm),"inp_min")
  w.min = calc.winning.natives(dat,-fitted(glm),"inp_min")

  all = data.frame(cbind(e,z.nat,w.nat,z.rep,w.rep,z.min,w.min))
  return(all)
}


make.termwise.zscores = function(dat=mat,glm=glmbk,nlabel="native"){
  scores = attr(glm$terms,"term.labels")
  scores = scores[scores!="id"]
  pdbplus = c(pdbs,"avg")
  z=array(dim=c(length(pdbplus),length(scores)),dimnames=list(pdbplus,scores))
  print("make.termwise.zscores")
  for (pdb in pdbs) {
    print(pdb)
    for (scoreterm in scores) {
      ds = dat[dat$id==pdb & dat$description=="output_decoy",scoreterm]
      ns = dat[dat$id==pdb & dat$description==nlabel,scoreterm]
      z[pdb,scoreterm] = zscore(ns,ds)
    }
  }
  for (scoreterm in scores) {
    z["avg",scoreterm] = mean(z[pdbplus!="avg",scoreterm],na.rm=T)
  }
  return(z)
}
