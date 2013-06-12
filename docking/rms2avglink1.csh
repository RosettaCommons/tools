#!/bin/csh -f

#
# fix headers on jeffs rms files
#
#awk '\\
#  { \\
#    if( NR == 1 && $1 == "pdbs" ) { \\
#      printf $2; \\
#      for(i=3;i<=NF;i++) printf " "$i; \\
#      printf "\n"; \\
#    } else { print $0 } \\
#   }' >\
#tmp.rms
awk '{print $0}' >tmp.rms

if ( $1 == "" ) then
  set rmscutoff=5
else
  set rmscutoff=$1
endif

#set pdb=$(basename $(dirname $PWD))
#set heyruns=$(basename $(dirname $(dirname $PWD)))


echo RMS cutoff for cluster definitions: $rmscutoff > tmp.err


#
#  run the R script
#
cat << END | R --vanilla  --silent > tmp.R.debug 

library(stats)

rms <- read.table("tmp.rms",row.names = 1,header=T)
rms.dist <- as.dist(rms)

#
#  run average linkage heirarchical clustering
#
rms.hclust <- hclust( rms.dist ,"average")

#
#  take all clusters at a distance of 5 rms
#
rms.cu <- sort( cutree( rms.hclust , h=$rmscutoff ) )


write.table( rms.cu , "tmp.clusters", quote=F , col.names=F )

#
#  output a graph
#
if (length(names(rms))<501) {
    postscript("cluster.ps",horizontal=T)
    if (substr(names(rms)[1],1,3)=="Set") {
	# this is a meta-cluster run
	rms.labels = as.vector(paste(substr(names(rms),4,5),"-",substr(names(rms),14,14),sep=""))
	target.number = substr(names(rms),6,6)[1]
	my.title = paste("METACLUSTER, ",date(),sep="")
    } else {
	# this is a regular run
	rms.labels = as.vector(paste(substr(names(rms),0,0),substr(names(rms),0,0),sep=""))
	target.number = substr(names(rms),6,6)[1]
	my.title = paste("Target ",target.number,", top",length(names(rms))," ",date(),sep="")
    }
    plot(rms.hclust,labels=rms.labels,main=my.title) #,cex=0.5
    rms.rect = rect.hclust(rms.hclust,h=$rmscutoff)
    dev.off()
}

END

#
#  send R output to stdout
#
cat tmp.clusters


#\rm tmp.clusters tmp.rms tmp.R.debug







