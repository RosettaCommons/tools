
#--------------------------------------------------------------
# load data into frame

#global.data.loaded = F
#fullatom = T

if (!global.data.loaded) {

  global.dataset = system("cd ../ ; basename $PWD",T)

  pdbs = as.list(system("ls *sc |sed 's/.sc//g'",T))
#  pdbs = list("1AHW","1BVK","1DQJ","1MLC","1WEJ")

  all.data = pdbs
  names(all.data) = pdbs

#  if (fullatom) extension = ".fasc" else extension = ".sc"
  extension = ".sc"

  i = 0
  for (pdb in pdbs) {
    
    i = i+1
    filename = paste(pdb,extension,sep="")
    print(paste("loading ",filename))
    all.data[[i]] = read.table(filename,header=T) 
  }
  global.data.loaded = T
}


# removed: rama, pc, pc.viol,dipolar
global.score.list = list("score","d.env","d.pair","contact","d.vdw", "d.fab",
			 "env","pair","vdw","hs","ss",
			 "sheet","cb","rsigma","hb","rg")

# removed: bk.prob
if (fullatom) global.score.list = 
  append(global.score.list,list("bk.tot","bk.atr","bk.rep","bk.sol","bk.hbsc",
				"bk.hbbb","bk.dun","bk.pair","gsolt",
                                "d.Eatr","d.Erep","d.Eatr.lr","d.Erep.lr"))
  

#global.score.list = list("reweight")

#--------------------------------------------------------------
# Makex one plot of x versus y, label the native
# additional plot arguments (title, xlab) can be passed through
decoyplot =  function(x,y,d="output_decoy",...){
  
# plot x versus y
# d is the description to pick out native and native repacked

  all = (d=="output_decoy"|
         d=="native"|d=="nat_repacked"|d=="nat_min"|d=="nat_mcm"|
         d=="inp_repacked"|d=="inp_min"|d=="inp_mcm")
  x.all = x[all]
  y.all = y[all]
  d.all = d[all]

  # drop the worst 5% from the range
  y.cutoff=y.all[order(y.all)[round(length(y.all)*0.95)]]
  subset = y.all<=y.cutoff | d.all!="output_decoy"
  if (sum(subset)<2) { # prevent total failures
    print("Empty plot subset")
    subset = T
  }
  x.all = x.all[subset]
  y.all = y.all[subset]
  
  
					# decoys in open blue circles
  plot(x[d=="output_decoy"],y[d=="output_decoy"],col="black",pch=".",
       xlim=range(x.all),ylim=range(y.all),...) 
    
  x.nat = x[d=="native"]
  y.nat = y[d=="native"]
  points(x.nat,y.nat,col="red",pch=4,cex=0.5) # native: red X
  
  x.rep = x[d=="nat_repacked"]
  y.rep = y[d=="nat_repacked"]
  points(x.rep,y.rep,col="green",pch=18,cex=0.5) # native-repacked: grn diamond
  
  x.irep = x[d=="inp_repacked"]
  y.irep = y[d=="inp_repacked"]
  points(x.irep,y.irep,col="cyan",pch=3,cex=0.5) # inp-repacked: cyan +

  x.mrep = x[d=="inp_min"]
  y.mrep = y[d=="inp_min"]
  points(x.mrep,y.mrep,col="blue",pch=22,cex=0.5) # inp-minimized:blue square

  x.mrep = x[d=="inp_mcm"|d=="nat_mcm"]
  y.mrep = y[d=="inp_mcm"|d=="nat_mcm"]
  points(x.mrep,y.mrep,col="violet",pch=25,cex=0.5) # inp-mcm:violet triangle
}

#--------------------------------------------------------------

# plot all scores for a single target
# input: mydata = the table read in from the .sc file
#                 e.g.: all.data[["1a2y"]] or all.data[[1]]
plot.target = function(mydata,title=""){

  par(mfrow = c(3, 5),oma=c(0,0,5,0))

  page.title = paste(title," ",global.dataset)

  score.list = global.score.list

  for (score.i in score.list) {

    x = mydata$rms
    d = mydata$description
    
    y = mydata[score.i]
    
    print(paste(" ",score.i))

    decoyplot(x,y[[1]],d,main=score.i,xlab="rms",ylab=score.i)
    mtext(page.title,line=1,side=3,outer=T,cex=1.5)
    
  }

}

#--------------------------------------------------------------

# compare all scores 
plot.score = function(all.data,score.i,title.prefix=""){

  page.title = paste(title.prefix,global.dataset,score.i)
#  page.title = global.dataset

  print(score.i)

  for (i in 1:length(all.data)) {

    pdb.name = names(all.data)[[i]]

    print(paste(" ",pdb.name))

    x = all.data[[i]]$rms
    d = all.data[[i]]$description
    
    y = all.data[[i]][score.i]
    
    decoyplot(x,y[[1]],d,main=pdb.name,xlab="rms",ylab=score.i)
    mtext(page.title,line=1,side=3,outer=T,cex=1.5)
    
  }
}


#--------------------------------------------------------------

# plot all scores for each target individually
# output in one big postscript
plot.targets.onefile = function(dataset=all.data,filename="targets.ps",title.prefix="") {

#  x11(,11,8.5)
#  par(ask=T)

  postscript(filename,horizontal=T,paper="letter")
  par(ask=F)
    
  for (pdb in pdbs) {
#    plotname = paste(pdb,".ps",sep="")
    print(pdb)
#    plot.title = paste(title.prefix,dataset[[i]])
    plot.target(dataset[[pdb]],pdb)

  }

  dev.off()

}

#--------------------------------------------------------------

# plot all scores for each target individually
# in individual files
plot.targets = function(dataset=all.data,title.prefix="") {

#  x11(,11,8.5)
#  par(ask=T)
 
  for (pdb in pdbs) {
    plotname = paste(pdb,".ps",sep="")
    print(plotname)
    postscript(plotname,horizontal=T,paper="letter")
    par(ask=F)
    plot.title = paste(title.prefix,pdb)
    plot.target(dataset[[pdb]],plot.title)
    dev.off()
  }
}

#--------------------------------------------------------------

# plot all targets for each score individually
plot.scores.onefile = function(dataset=all.data,filename="scores.ps",title.prefix=""){

##choose one:
#  x11(,11,8.5) # on-screen
  postscript(filename,horizontal=T,paper="letter")

  par(mfrow = c(3, 5),oma=c(0,0,5,0))
  
#   score.list = list("score","d.env","d.pair","contact","d.vdw",
# 		    "env","pair","vdw","hs","ss",
# 		    "sheet","cb","rsigma","hb","rg",
# 		    "co","rama","pc","pc.viol","dipolar",
# 		    "bk.tot","bk.atr","bk.rep","bk.sol","bk.hbsc",
# 		    "bk.hbbb","bk.dun","bk.prob","bk.pair")
  
  score.list = global.score.list

  for (score.i in score.list) plot.score(dataset,score.i,title.prefix)

  dev.off()

}


#--------------------------------------------------------------

# plot all targets for each score individually
plot.scores = function(dataset=all.data,filebase="scores",title.prefix=""){

##choose one:
#  x11(,11,8.5) # on-screen
  
#   score.list = list("score","d.env","d.pair","contact","d.vdw",
# 		    "env","pair","vdw","hs","ss",
# 		    "sheet","cb","rsigma","hb","rg",
# 		    "co","rama","pc","pc.viol","dipolar",
# 		    "bk.tot","bk.atr","bk.rep","bk.sol","bk.hbsc",
# 		    "bk.hbbb","bk.dun","bk.prob","bk.pair")
  
  score.list = global.score.list

  plots.per.line = 5
  lines.per.page = 4
  scores.per.file = 1
  nplots = scores.per.file

  postscript()
  
  for (score.i in score.list) {
    if (nplots>=scores.per.file) {
      dev.off()
      filename = paste(filebase,score.i,"ps",sep=".")
      postscript(filename,horizontal=T,paper="letter")
      par(mfrow = c(lines.per.page,plots.per.line), 
          oma=c(0,0,5,0),mar=c(2.5,2.5,2,1.5),mgp=c(1.5,0.5,0))
      nplots = 0
    }
    plot.score(dataset,score.i,title.prefix)
    nplots = nplots + 1
  }
  dev.off()

}


#--------------------------------------------------------------
# filter out certain data
filter.data = function(indata){

  criteria = indata$d.vdw < 10 & indata$contact < 10
  outdata = indata[criteria,]

  print(paste("points removed: ", dim(indata)[1] - sum(criteria)))

  return(outdata)

}


if (!global.filtered.data) {
  filtered.data = lapply(all.data,filter.data)
  global.filtered.data=T
}


#--------------------------------------------------------------
#other diagnostics

#rms of the 12 best structures
smallest=function(dataset){sort(dataset$rms)[65:77]}
#starting point rms of the 12 best structures
smallest.start=function(dataset){
  cutoff=sort(dataset$rms)[77]
  return(dataset$st.rm[dataset$rms < cutoff & dataset$rms != 0.0])
}

# mean of rms and starting rms of the 12 best structures
smallest.mean=function(dataset){mean(sort(dataset$rms)[65:77])}
smallest.start.mean=function(dataset){
  cutoff=sort(dataset$rms)[77]
  return(mean(dataset$st.rm[dataset$rms < cutoff & dataset$rms != 0.0]))
}

# ratio of the above
improvement=function(dataset){return(smallest.start.mean(dataset)-smallest.mean(dataset))}

#lapply(all.data,improvement)
#lapply(all.data,smallest.start.mean)
#lapply(all.data,smallest.mean)
#lapply(all.data,smallest.start)
#lapply(all.data,smallest)


#--------------------------------------------------------------
## selecting a single dataset

pdbdata = function(pdb.name) {
  sel = pdbs==pdb.name
  return (all.data[sel][[1]])
}

showpoint = function(pdb,n) {
  point = all.data[pdb][[1]][n,]
#  print(point$filename)
  return(point)
}

singleplot = function(pdb,score.i,idn=0,data=all.data) {
  mydata = data[pdb][[1]]
  x = mydata$rms
  d = mydata$description
  y = mydata[score.i]

  title = paste(pdb,score.i)
  print(title)

  decoyplot(x,y[[1]],d,main=title,xlab="rms",ylab=score.i)

  if (idn!=0) {
    points = identify(x,y[[1]],n=idn)
    showpoint(pdb,points[1])
  }
}

#singleplot("1brs","hb")
#all.data$"1brs"[167,]



rms.histo = function(mydata,lowonly=F,...) {
  
  d = mydata$description
  x = mydata$rms[d=="output_decoy"]
  y = mydata$st.rm[d=="output_decoy"]

  if (lowonly) brks = c(0:6,100)*2
  else brks = c(0:(max(range(x),range(y))/2+2))*2
  
  rmhist = hist(as.vector(x),plot=F,breaks=brks) 
  sthist = hist(as.vector(y),plot=F,breaks=brks) 

  df = t(data.frame(start=sthist$counts[1:(length(sthist$counts)-lowonly)],
                    final=rmhist$counts[1:(length(sthist$counts)-lowonly)]))
  
  if (lowonly) colnames(df)= brks[-c(1,length(brks))]-1
  else colnames(df)= brks[-1]-1

  barplot(df,beside=T,legend=rownames(df),xlab="rms",ylab="counts",...)
}

rms.histos = function (dataset=all.data,printme=F,filename="rms.histo.ps",...) {

  if (printme) postscript(filename,horizontal=T,paper="letter")
  par(mfrow = c(2,3),oma=c(0,0,5,0))
  page.title = global.dataset

  for (i in 1:length(dataset)) {

    pdb.name = names(dataset)[[i]]

    rms.histo(dataset[[i]],main=pdb.name,...)

  }

  mtext(page.title,line=1,side=3,outer=T,cex=1.5)
  if (printme) dev.off()
}

rms.histos.print = function(dataset=all.data,...) {
  rms.histos(dataset,printme=T,filename="rms.histo.ps",...)
  rms.histos(dataset,printme=T,filename="rms.histo.low.ps",lowonly=T,...)
}




#----------------------------


# #how to load from a different directory

#   global.dataset = system("cd ../ ; basename $PWD",T)

#   pdbs = as.list(system("ls *sc |sed 's/.sc//g'",T))
# #  pdbs = list("1AHW","1BVK","1DQJ","1MLC","1WEJ")

#   all.data = pdbs
#   names(all.data) = pdbs

# #  if (fullatom) extension = ".fasc" else extension = ".sc"
#   extension = ".sc"
#   path="../../calrong-4-26BK/scorefiles/"

#   i = 0
#   for (pdb in pdbs) {
    
#     i = i+1
#     filename = paste(path,pdb,extension,sep="")
#     print(paste("loading ",filename))
#     all.data[[i]] = read.table(filename,header=T) 
#   }
#   global.data.loaded = T

# how to plot two sets to compare

#--------------------------------------------------------------
# Makex one plot of x versus y, label the native
# additional plot arguments (title, xlab) can be passed through
decoyplot2 =  function(x,y,d="output_decoy",...){
  
# plot x versus y
# d is the description to pick out native and native repacked

					# decoys in open blue circles
  points(x[d=="output_decoy"],y[d=="output_decoy"],col="plum1",pch=".",...) 
    
  x.nat = x[d=="native"]
  y.nat = y[d=="native"]
  points(x.nat,y.nat,col="plum1",pch=4,cex=0.5) # native: red X
  
  x.rep = x[d=="nat_repacked"]
  y.rep = y[d=="nat_repacked"]
  points(x.rep,y.rep,col="pink",pch=18,cex=0.5) # native-repacked: grn diamond
  
  x.irep = x[d=="inp_repacked"]
  y.irep = y[d=="inp_repacked"]
  points(x.irep,y.irep,col="lightskyblue",pch=3,cex=0.5) # inp-repacked: cyan +

  x.mrep = x[d=="inp_min"]
  y.mrep = y[d=="inp_min"]
  points(x.mrep,y.mrep,col="blue4",pch=22,cex=0.5) # inp-minimized:blue square
}

plottwo = function(set1,set2,...){

  decoyplot(set1[[1]]$rms,set1[[1]]$score,set1[[1]]$description,...)
  decoyplot2(set2[[1]]$rms,set2[[1]]$score,set2[[1]]$description,...)

}
