#!/usr/bin/env Rscript
# kink_plots.R
# plot geometry of the H3 kink for a set of modeled antibody structures
# J. Gray June 2013

# make a rectangle around the desired range and count points within
rect.range = function(xmin,xmax,x,ymin,ymax,y,placement="bottomleft") {
	rect(xmin,ymin,xmax,ymax, border="blue")
	count = sum( x>xmin & x<xmax & y>ymin & y<ymax )
	fraction = round ( count / length(x), 3 )
	legend(placement,legend=fraction)
}

##### qbase #####

plot.qbase = function(scfile,top10,native,...){
	plot( scfile$kink_qbase,scfile$kink_q, col="black", pch=".",ylim=c(0,10),
		xlim=c(-180,180),...)
	points( top10$kink_qbase, top10$kink_q, col="green", pch=18, cex=0.5) # blue diamond
	points(native$qbase,native$q, col="red", pch=4) # red X
	rect.range(-10,70,scfile$kink_qbase,6.5,7.75,scfile$kink_q)
}

plot.native.qbase = function(native){
	plot( native$qbase,native$q, col="red", pch=".",ylim=c(0,10),
	  	 xlim=c(-180,180),main="XTALS")
	rect.range(-10,70,native$qbase,6.5,7.75,native$q)
}

# HBonds
#rect.HB = function() { rect(-10,2,70,4, border="blue") } #0.523 0.698 = 30+/-40

##### bb HB #####
plot.bbHB = function(scfile,top10,native,...){
	plot( scfile$kink_qbase,scfile$kink_bb_HB, col="black", pch=".",ylim=c(0,10),
		xlim=c(-180,180),...)
	points( top10$kink_qbase, top10$kink_bb_HB, col="green", pch=18, cex=0.5) # blue diamond
	points(native$qbase,native$bbHBdist, col="red", pch=4) # red X
	rect.range(-10,70,scfile$kink_qbase,2,4,scfile$kink_bb_HB)
}

plot.native.bbHB = function(native){
	plot( native$qbase,native$bbHBdist, col="red", pch=".",ylim=c(0,10),
	  	 xlim=c(-180,180),main="XTALS")
	rect.range(-10,70,native$qbase,2,4,native$bbHBdist)
}

##### RD HB #####
plot.RDHB = function(scfile,top10,native,...){
	plot( scfile$kink_qbase,scfile$kink_RD_HB, col="black", pch=".",ylim=c(0,10),
		xlim=c(-180,180),...)
	points( top10$kink_qbase, top10$kink_RD_HB, col="green", pch=18, cex=0.5) # blue diamond
	points(native$qbase,native$HBdist, col="red", pch=4) # red X
	rect.range(-10,70,scfile$kink_qbase,2,4,scfile$kink_RD_HB)
}

plot.native.RDHB = function(native){
	plot( native$qbase,native$HBdist, col="red", pch=".",ylim=c(0,10),
	  	 xlim=c(-180,180),main="XTALS")
	rect.range(-10,70,native$qbase,2,4,native$bbHBdist)
	}


##### Trp HB #####
plot.TrpHB = function(scfile,top10,native,...){
	plot( scfile$kink_qbase,scfile$kink_Trp_HB, col="black", pch=".",ylim=c(0,10),
		xlim=c(-180,180),...)
	points( top10$kink_qbase, top10$kink_Trp_HB, col="green", pch=18, cex=0.5) # blue diamond
	points(native$qbase,native$W_HBdist, col="red", pch=4) # red X
	rect.range(-10,70,scfile$kink_qbase,2,4,scfile$kink_Trp_HB)
}

plot.native.TrpHB = function(native){
	plot( native$qbase,native$W_HBdist, col="red", pch=".",ylim=c(0,10),
	  	 xlim=c(-180,180),main="XTALS")
	rect.range(-10,70,native$qbase,2,4,native$W_HBdist)
}


##### Score v RMSD #####
plot.scoreVrmsd = function(scfile,top10,...){
	# correct for q constraint
	#sc = scfile$total_score - scfile$dihedral_constraint
	sc = scfile$unconstr_score
	sc.max=sc[order(sc)[round(length(sc)*0.95)]]  # cut top 5% of scores from plot
	plot( scfile$H3_RMS,sc, col="black", pch=".",
		xlim=c(0,18),ylim=c(min(sc),sc.max),...)
	points( top10$H3_RMS, top10$unconstr_score, col="green", pch=18, cex=0.5) # blue diamond
	rect.range(par()$xaxp[1],3,top10$H3_RMS,
		       min(sc), max(top10$unconstr_score),top10$unconstr_score,
		       "bottomright")
}



##### Score v qbase #####
plot.scoreVqbase = function(scfile,top10,...){
	# correct for q constraint
	#sc = scfile$total_score - scfile$dihedral_constraint
	sc = scfile$unconstr_score
	sc.max=sc[order(sc)[round(length(sc)*0.95)]]  # cut top 5% of scores from plot
	plot( scfile$kink_qbase,sc, col="black", pch=".",
		xlim=c(-180,180),ylim=c(min(sc),sc.max),...)
	points( top10$kink_qbase, top10$unconstr_score, col="green", pch=18, cex=0.5) # blue diamond
	rect.range(-10,70,top10$kink_qbase,min(top10$unconstr_score)-1, max(top10$unconstr_score),top10$unconstr_score,"topleft")
}



##### load data #####
load.directory = function(direc=""){
	olddir = getwd()
	setwd(direc)

	scfile = read.table("remodel_h3.fasc.sort",header=T)
	names(scfile)
	top10 = read.table("remodel_h3.fasc.sort.top10",header=T)
	native = read.table("kink_geom.dat",header=T)
	names(native)

	setwd(olddir)
	return(list(scfile=scfile,top10=top10,native=native))
}



setup.dev = function(dtype="pdf"){
	if (dtype=="quartz") {
		quartz(height=11,width=8.5)
	}
	else { # default to pdf
		pdf(file="kink_plots.pdf",height=11,width=8.5,paper="letter")
   	}
	par(mfrow = c(5,3),oma=c(0,0,5,0))
	par(mar=c(3,3,2,1))
	par(mgp=c(2,1,0))
	return(dev.cur())
}


pageofplots = function(page.title,plot.native,plot.decoys,pdbs,natives) {
	par(mfrow = c(5,3))
	print(page.title)
	plot.native(natives)
	for (pdb in pdbs) {
		print(pdb)
		target = load.directory(pdb)
		plot.decoys(target$scfile,target$top10,target$native,main=pdb)
	}
	mtext(page.title,line=1,side=3,outer=T,cex=1.5)
}

pageoffunnels = function(page.title,plot.decoys,pdbs) {
	par(mfrow = c(5,3))
	plot.new() # to skip native pane
	print(page.title)
	for (pdb in pdbs) {
		print(pdb)
		target = load.directory(pdb)
		plot.decoys(target$scfile,target$top10,main=pdb)
	}
	mtext(page.title,line=1,side=3,outer=T,cex=1.5)
}



main = function() {
	dir = system("basename $PWD",T)
	pdbs = list("2R56","2VXS","2VXV","3AAZ","3G04","3GI9","3GNM",
				"3HNT","3HZV","3I9G","3IJH","3P0Y","3V6O")
	#pdbs = list("2R56","2VXS","2VXV","3AAZ","3G04","3GI9","3GNM",
	#			"3HNT","3HZV","3I9G","3P0Y")
	natives = read.table("/Users/jeff/git/Rosetta/tools/antibody/antibody_database/kink_geom.PDB.dat",header=T)

	setup.dev()

 	pageofplots(paste(dir,"qbase"),plot.native.qbase,plot.qbase,pdbs,natives)
 	pageofplots(paste(dir,"RD sc"),plot.native.RDHB,plot.RDHB,pdbs,natives)
 	pageofplots(paste(dir,"RD bb"),plot.native.bbHB,plot.bbHB,pdbs,natives)
 	pageofplots(paste(dir,"Trp"),plot.native.TrpHB,plot.TrpHB,pdbs,natives)
	pageoffunnels(paste(dir,"score v rmsd"),plot.scoreVrmsd,pdbs)
	pageoffunnels(paste(dir,"score v qbase"),plot.scoreVqbase,pdbs)

	print("finished")
	dev.off()
}

if ( ! interactive() ) { main() }

