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

#--------------------------------------------------------------

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

#--------------------------------------------------------------

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


