
plot.one = function(scfile,top10,native,title="",...){

	plot( scfile$kink_qbase,scfile$kink_bb_HB, col="black", pch=".",ylim=c(0,10), 	
		xlim=c(-180,180),...)
	points( top10$kink_qbase, top10$kink_bb_HB, col="pink", pch=18, cex=0.5) # blue diamond
	points(native$qbase,native$bbHBdist, col="red", pch=4) # red X

	rect(-10,2,70,4, border="blue")
	#0.523 0.698 = 30+/-40
}


scfile = read.table("remodel_h3.fasc.sort",header=T)
names(scfile)
top10 = read.table("remodel_h3.fasc.sort.top10",header=T)

native = read.table("kink_geom.dat",header=T)
names(native)


#par(mfrow = c(3, 5),oma=c(0,0,5,0))
page.title = "HB"
plot.one(scfile,top10,native,"HB",main="2R56")
mtext(page.title,line=1,side=3,outer=T,cex=1.5)
