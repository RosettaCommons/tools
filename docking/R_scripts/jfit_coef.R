# set up size correctly
#options("object.size" = 2e+09)

# alpha proteins no need for hs, ss, rsigma, sheet


# fullatom = F
# Allscores.loaded = F

# --------------- prep logistic regression --------------

# load data from file, if we haven't already
if (!Allscores.loaded) {
  global.dataset = system("cd ../ ; basename $PWD",T)
  print("loading Allscores")
  mat = read.table("Allscores",header=T)
  Allscores.loaded = T
}
mat$id = as.factor(mat$id)

print("setting up...") #--------------------------------

# get pdb names
pdbs = unique(mat$id)
pdbs = as.character(pdbs)
pdbs = sort(pdbs)

# set up cutoff fractions
rmsfrac = 0.05 # cut rms at top 5% -- definition of "good"

# enrichment cutoffs 
scfrac = 0.15          # look for scores in the top 15% (enrichment default)
scpercs = c(10,1,0.1)  # look at 10, 1, and 0.1% top scores in enrichments

# calculate the 5% cutoff RMS for each PDB
# also, find any incomplete PDBs and mark for removal
rmscut = array(0,length(pdbs))
names(rmscut)=pdbs
bad.pdbs = 0 # this is a list of bad pdb numbers
for(j in 1:length(pdbs)) {
  print(paste("sizing up pdb ",pdbs[j]," --",
              sum(mat$id==pdbs[j] & mat$description=="output_decoy"),
              " decoys"))
  rmsd = mat$rms[mat$id==pdbs[j] & mat$description=="output_decoy"]
  sel = round(length(rmsd)*rmsfrac)
  if (sel==0) { # detect bad pdbs -- those without decoys
    rmscut[j] = 0
    bad.pdbs=c(bad.pdbs,j)
    print(paste("bad pdb: ",pdbs[j]," decoys: ",length(rmsd)))
  }
  else { # these are good pdbs
    rmscut[j] = rmsd[order(rmsd)][sel]
  }
}

# get rid of bad-pdb lines
mat=mat[rmscut[mat$id]!=0,]

# data filtering
mat=mat[mat$contact<9.5,]  # must be in contact

# the 'good' column tells whether the decoy is within the top 5%
mat$good=as.numeric(mat$rms<=rmscut[mat$id])
good.decoys=(1:dim(mat)[1])[mat$good==1] # a list of decoy numbers

# remove bad pdbs in lists
if (sum(bad.pdbs)) {
  pdbs = pdbs[-bad.pdbs]
  rmscut = rmscut[-bad.pdbs]
}
print("Processing these PDBs: ")
print(pdbs)



# regressions--------------------------------------------------

# all scores
if (!global.regressions.done) {
  # 2/8/2: I used to do all the regressions at once to compare; however,
  # most data sets are too large now, and you can only do one regression
  # per session due to memory limits.  Therefore, this section is no longer
  # regularly used                     
  print("producing linear model")
  print ("backbone regression") 
  glmbb.d <- glm(good ~ d.pair + d.env + contact + d.vdw + d.fab + id, 
                 data = mat,family="binomial")
  glmbb <- glm(good ~ d.pair + d.env + contact + d.vdw + d.fab + 
               env + pair + hs + ss + sheet + cb + rsigma + hb + rg + id, 
               data = mat,family="binomial")
  print("glmbb and glmbb.d available")
  if (fullatom) 
    {
      print ("fullatom regression") 
#   glm1 <- glm(good ~ d.pair + d.env + contact + d.vdw +
# 	      env + pair + vdw + hs + ss + sheet + cb + rsigma + hb + 
# 	      rg + co + rama + pc + pc.viol + dipolar +
# 	      bk.atr + bk.rep + bk.sol + bk.hbbb + bk.hbsc + bk.dun +
# 	      bk.prob + bk.pair + id, 
# 	      data = mat,family="binomial")
      glmfa <- glm(good ~ d.pair + d.env + contact + d.vdw + d.fab +
                   env + pair + hs + ss + sheet + cb + rsigma + hb + rg + 
                   bk.atr + bk.rep + bk.sol + bk.hbbb +
                   bk.hbsc + bk.dun + bk.pair 
                   + gsolt + id, 
                   data = mat,family="binomial")
      glmbk <- glm(good ~ 
                   bk.atr + bk.rep + bk.sol + bk.hbbb +
                   bk.hbsc + bk.dun + bk.pair 
                   + gsolt + id, 
                   data = mat,family="binomial")
      glm.10d1 <- glm(good ~ d.fab +
                      bk.rep + bk.hbbb + bk.hbsc + bk.pair 
                      + gsolt + id, 
                      data = mat,family="binomial")
      glm.10d2 <- glm(good ~ bk.rep + bk.sol + bk.hbbb + bk.hbsc + bk.dun
                      + gsolt + id,
                      data = mat,family="binomial")
      glm.10d3 <- glm(good ~ bk.atr +
                      bk.rep + bk.sol + bk.hbbb + bk.hbsc + bk.dun
                      + gsolt + id,
                      data = mat,family="binomial")
      print("glmfa and glmbk and glm.10d1 available")
    }
  global.regressions.done=T
}

# this function is also out of date...the glm call is now directly
# in the do=fit*.R scripts (4/29/2)
make.glmbk = function(){
  print ("fullatom regression") 
  glm(good ~ 
      bk.atr + bk.rep + bk.sol + bk.hbbb + bk.hbsc + bk.dun + bk.pair 
      + gsolt + id, 
      data = mat,family="binomial")
}




# plotting ------------------------------------------------------------


plot.fit = function(glm=glmbk,filename="fit.ps",dat=mat){
# to plot:
  
# only one:
# plot(dat$rms[dat$id=="1a32"],fitted(lm1)[dat$id=="1a32"])

    # page setup
    plots.per.line = 5
    lines.per.page = 4
    postscript(filename, horizontal = T,paper="letter")
      par(mfrow = c(lines.per.page,plots.per.line), 
          oma=c(0,0,5,0),mar=c(2.5,2.5,2,1.5),mgp=c(1.5,0.5,0))
    print(filename)

    # plots
    i=0
    for(pdb in pdbs) {
      i=i+1
      print(pdb)
      # get dat
      x = dat$rms[dat$id==pdb]
      y = fitted(glm)[dat$id==pdb]
      plot(x,y,xlab="rms",ylab="fitted",main=pdb,pch=".")
      # rms cutoff points
      rmss = dat$id==pdb & dat$rms<=rmscut[pdb] & dat$description=="output_decoy"
      xlow = dat$rms[rmss]
      ylow = fitted(glm)[rmss]
      points(xlow,ylow,col="red",pch=".")
      # natives
      natives = dat$id==pdb & dat$description=="native"
      x.nat = dat$rms[natives]
      y.nat = fitted(glm)[natives]
      points(x.nat,y.nat,col="pink",pch=4,cex=0.5) # native: red/pink X
      
      reps = dat$id==pdb & dat$description=="nat_repacked"
      x.rep = dat$rms[reps]
      y.rep = fitted(glm)[reps]
      points(x.rep,y.rep,col="green",pch=18,cex=0.5) # native-repacked: grn diamond

      ireps = dat$id==pdb & dat$description=="inp_repacked"
      x.irep = dat$rms[ireps]
      y.irep = fitted(glm)[ireps]
      points(x.irep,y.irep,col="cyan",pch=3,cex=0.5) # inp-repacked: cyan +
      
      mreps = dat$id==pdb & dat$description=="inp_min"
      x.mrep = dat$rms[mreps]
      y.mrep = fitted(glm)[mreps]
      points(x.mrep,y.mrep,col="blue",pch=22,cex=0.5) # inp-minimized:blue square
      
      # title and newpage, if necessary
      if ((i%%(plots.per.line*lines.per.page))==1) {  # add title
        ti = paste("log regr. fit vs rms",global.dataset,"5% cutoff",sep=" ")
        mtext(ti,line=2,side=3,outer=T,cex=1.5)
        mtext(glm$call,line=0,side=3,outer=T,cex=0.7)
      }
      # enrichment
      #scpercs=c(10,1,0.1) # score percent cutoffs -- now global
      lcols=c("darkgreen","purple","blue")
      decoys = dat$id==pdb & dat$description=="output_decoy"
      x.decoys = dat$rms[decoys]
      y.decoys = fitted(glm)[decoys]
      for(j in 1:3) {
        scperc=scpercs[j]
        lcol=lcols[j]
        lines(c(rmscut[pdb],rmscut[pdb]),c(min(y),max(y)),col="red")
        sc.cut=-calc.percent.cutoff(-y.decoys,scperc)
        lines(c(0,max(x)),c(sc.cut,sc.cut),col=lcol)
        pdb.enrich=round(log.enrichn(x.decoys,-y.decoys,scperc,rmsfrac*100),3)
        enrich.text=paste("logE =",pdb.enrich)
        mtext(enrich.text,side=3,line=-4+j,adj=0.95,cex=0.6,col=lcol)
      }
    }
    graphics.off()
}

#---------------------

jwrite.table = function(x,file="tmp",...) {
  system(paste('echo ',file,' > tmp'))
  system('basename $(dirname $PWD) >> tmp')
  system('date >> tmp')
  system('echo ---------------------------- >> tmp')
  system('echo -en "\t" >> tmp')
  system(paste('mv tmp ',file))
  write.table(x,file=file,quote=F,sep='\t',append=T,...)
}



#--------------------- a subroutine to process the fit ------------------

process.glm = function(glm=glmbk,label="fit.bk"){

  print(paste("process.glm -- ",label))
  # spit out the coefficients to a file
  summary(glm)
  coefs = summary(glm)$coefficients
  rcoefs = coefs # rescale coefficients by 10*repulsive score
  rcoefs[,"Estimate"] = rcoefs[,"Estimate"]/glm$coefficients["bk.rep"]/10
  rcoefs[,"Std. Error"] = rcoefs[,"Std. Error"]/glm$coefficients["bk.rep"]/10
  jwrite.table(round(coefs ,5),file=paste(label,".coefficients",sep=""))
  jwrite.table(round(rcoefs,5),file=paste(label,".coefficients.rep1",sep=""))

  # get the enrichments, z-scores, and winners, and spit to a file
  glm.fit.summary = make.fit.summary(mat,glm)
  jwrite.table(round(glm.fit.summary,2),file=paste(label,".summary",sep=""))

  termwise.zscores.nat = round(make.termwise.zscores(mat,glm),3)
  termwise.zscores.rep = round(make.termwise.zscores(mat,glm,"inp_repacked"),3)
  termwise.zscores.min = round(make.termwise.zscores(mat,glm,"inp_min"),3)
  jwrite.table(termwise.zscores.nat,file=paste(label,".ztable",sep=""))
  jwrite.table(termwise.zscores.rep,file=paste(label,".ztable.repacked",sep=""))
  jwrite.table(termwise.zscores.min,file=paste(label,".ztable.minimized",sep=""))
  
  # make a pretty file with coefficients and the summary
  system(paste("enscript ",label,".coefficients.rep1 ",
               label,".summary -G -o tmp.ps",sep=""))
  system(paste("psnup -2 tmp.ps > ",label,".summary.ps",sep=""))
  system("rm tmp.ps")
  
  # make a plot
  plot.fit(glm,paste(label,".ps",sep=""))

}
