
####### jerry's enrichment functions, modified #########

# JJG Jan 02

calc.percent.cutoff=function(x,perc){
  ntotal  = length(x)
  ntop    = round(ntotal*perc/100)     # number of things below cutoff
  cut     = x[order(x)][ntop]          # x at the cutoff
  return(cut)
}
  

# ie enrich(mat[24],mat[23])
# ie enrich(mat$rms,mat$score)
log.enrich.rms=function(rmsd,score,ptop=15,rmsd.cut=5){
    # arguments: ptop     = percentage of decoys that are 'top scoring'
    #            rmsd.cut = rms which defines decoys 'close' to native
    # result: enrichment of close decoys in top ptop percent of scores
  ntotal  = length(rmsd)               # number of decoys total
  close   = rmsd < rmsd.cut            # logical indicating close rmsds
  nclose  = sum(close)                 # number of 'close' decoys
  ntop    = round(ntotal*ptop/100)     # number of 'top scoring' decoys
  sc.cut  = score[order(score)][ntop]  # score at the cutoff
  top     = score < sc.cut             # logical indicating top scores
  both    = close & top                # logical indicating top & close decoys
  nboth   = sum(both)                  # number of top & close decoys
  enrichment = nboth*ntotal/ntop/nclose
  if (nclose==0) enrichment=0
  if (ntop==0) enrichment=-1
  return(log(enrichment))
}

# ie enrich(mat[24],mat[23])
# ie enrich(mat$rms,mat$score)
log.enrichn=function(rmsd,score,ptop=15,pclose=5){
    # arguments: ptop   = percentage of decoys that are 'top scoring'
    #            pclose = percentage of decoys that are 'close' to native
    # result: enrichment of close decoys in top sperc percent of scores
  ntotal  = length(rmsd)               # number of decoys total
  nclose  = round(ntotal*pclose/100)   # how many decoys are 'close'
  rmsd.cut= rmsd[order(rmsd)][nclose]  # rms at the cutoff
  close   = rmsd < rmsd.cut            # logical indicating close rmsds
  ntop    = round(ntotal*ptop/100)     # how many decoys are 'top scoring'
  sc.cut  = score[order(score)][ntop]  # score at the cutoff
  top     = score < sc.cut             # logical indicating top scores
  both    = close & top                # logical indicating top & close decoys
  nboth   = sum(both)
  enrichment = nboth*ntotal/ntop/nclose
  if (nclose==0) enrichment=0
  if (ntop==0) enrichment=-1
  return(log(enrichment))
}

# this is jerry's old function
# ie enrich(mat[24],mat[23])
# ie enrich(mat$rms,mat$score)
enrichnold=function(rmsd,score,sperc=15,rperc=5,rd=3){
    # arguments: sperc = percentage of decoys that are 'top scoring'
    #            rperc = percentage of decoys that are 'close' to native
    # result: enrichment of close decoys in top sperc percent of scores
  selr=round(length(rmsd)*rperc/100) # how many decoys are 'close'
  sels=round(length(rmsd)*sperc/100) # how many decoys are 'top scoring' 
  best=order(rmsd)[1:selr]           # the indices of the close decoys
  rmsd.cut=rmsd[order(rmsd)][selr]   # rms at the cutoff
  which.sel=order(score)[1:sels]     # the indices of the top scoring decoys 
  both=match(best,which.sel)         # the indices of close, top scoring decoys
  both=both[both!="NA"]              # without the crap    
  zz=round(length(both)/(length(rmsd)*(sperc/100)*(rperc/100)),rd)
  return(zz)
}

log.enrich=function(rmsd,score,perc=15) {
    # here perc is sperc and rperc in enrichn
  return (log.enrichn(rmsd,score,perc,perc))
}

####### end jerry's enrichment functions #########


calc.enrichments = function(glm=glmbk,rperc=rmsfrac*100,sperc=scfrac*100){
  
  # plots
  i=0
  e=array(0,length(pdbs))
  names(e)=pdbs
  sccut=array(0,length(pdbs))
  names(sccut)=pdbs
  for(pdb in pdbs) {
    i=i+1
    # get mat
    x = mat$rms[mat$id==pdb & mat$description=="output_decoy"]
    y = fitted(glm)[mat$id==pdb & mat$description=="output_decoy"]
    # enrichment
    e[pdb]=log.enrichn(x,-y,sperc,rperc)
#    # cutoffs, for reference
#    sccut[pdb]=-calc.percent.cutoff(-y,sperc)
#    rmscut[pdb]=calc.percent.cutoff(x,rperc)
  }
#  print(sccut)
#  print(rmscut)
  edat=e
  # enrichments are in log form
  e["avg"]=mean(edat[is.finite(edat)])
#  e["logavg"]=exp(mean(log(edat[edat!=0])))
  e["bad"]=sum(edat<0)
  e["fail"]=sum(is.infinite(edat))
  e["n"]=length(edat)
  return(e)
}

calc.enrichments.rms = function(glm=glmbk,rms.cut=5,sperc=scfrac*100){
  
  # plots
  i=0
  e=array(0,length(pdbs))
  names(e)=pdbs
  for(pdb in pdbs) {
    i=i+1
    # get mat
    x = mat$rms[mat$id==pdb]
    y = fitted(glm)[mat$id==pdb]
    # enrichment
    sperc=15
    e[pdb]=enrich.rms(x,-y,15,rms.cut)
  }
  e["avg"]=mean(e,na.rm=T)
  e["fail"]=sum(e<1)
  e["totalfail"]=sum(e==0)
  e["n"]=length(e)-3
  return(e)
}



calc.multi.enrichments = function(glm=glmbk,spercs=scpercs){
#  spercs=c(10,1,0.1)
  rperc=rmsfrac*100
  many.e = sapply(spercs,calc.enrichments,glm=glm,rperc=rperc)
  dimnames(many.e)[[2]]=paste("E",spercs,"pct",sep="")
  return(many.e)
}

