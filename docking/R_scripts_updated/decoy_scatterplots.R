#--------------------------------------------------------------
# Make one plot of x versus y, label the native
# additional plot arguments (title, xlab) can be passed through
# d is the description to pick out native and native repacked

decoyplot =  function(x,y,d="output_decoy",...){
    

    all = (d=="output_decoy" | d=="native" |
        d=="nat_repacked" | d=="nat_min" | d=="nat_mcm" |
        d=="inp_repacked" | d=="inp_min" | d=="inp_mcm")

    x.all = x[all]
    y.all = y[all]
    d.all = d[all]

    #TODO: FLAG CONTROLLED
        # drop the worst 5% from the range
        y.cutoff=y.all[order(y.all)[round(length(y.all)*0.95)]]
        subset = y.all<=y.cutoff | d.all!="output_decoy"
 	
        # prevent total failures
        if (sum(subset)<2) {
            print("Empty plot subset")
	        subset = T
        }
        x.all = x.all[subset]
        y.all = y.all[subset]
    
    # decoys in black dots
    x.dec = x[d=="output_decoy"]
    y.dec = y[d=="output_decoy"]
    plot(x.dec, y.dec, col="black", pch=".", xlim=range(x.all), ylim=range(y.all),...) 
 
    # native: red X
    x.nat = x[d=="native"]
    y.nat = y[d=="native"] 
    points(x.nat,y.nat,col="red",pch=4,cex=0.5)
  
    # native-repacked: grn diamond
    x.rep = x[d=="nat_repacked"]
    y.rep = y[d=="nat_repacked"]
    points(x.rep,y.rep,col="green",pch=18,cex=0.5)
  
    # inp-repacked: cyan +
    x.irep = x[d=="inp_repacked"]
    y.irep = y[d=="inp_repacked"]
    points(x.irep,y.irep,col="cyan",pch=3,cex=0.5)

    # inp-minimized:blue square
    x.imin = x[d=="inp_min"]
    y.imin = y[d=="inp_min"]
    points(x.imin, y.imin, col="blue",pch=22,cex=0.5)

    # inp-mcm:violet triangle
    x.imcm = x[d=="inp_mcm"|d=="nat_mcm"]
    y.imcm = y[d=="inp_mcm"|d=="nat_mcm"]
    points(x.imcm, y.imcm, col="violet", pch=25, cex=0.5) 
}

#--------------------------------------------------------------
# plot all scores for a single target
# input: mydata = the table read in from the score file
# e.g.: all.data[["1a2y"]] or all.data[[1]]

plot.target = function(mydata, score.list, x_name, oneFilePerScoreTerm, title ){
    for (score.i in score.list) {
        x = mydata[x_name]
        print( score.i)
        y = mydata[score.i]
        d = mydata$description
       	#d = mydata$type 
        if(oneFilePerScoreTerm){
	        dev.off()
	        plotname = paste(score.i,"ps",sep=".")
            postscript(plotname,horizontal=T,paper="letter")
        }
        if ( length( score.list ) > 1 ){
            plot.title = score.i
        } else {
            plot.title = title
        }
        decoyplot(x[[1]],y[[1]],d,main=plot.title,xlab=x_name,ylab=score.i)
    }
}

#--------------------------------------------------------------
# plot all scores for each target in individual files

plot.targets = function(dataset, pdbs, score.list, x_name, oneFilePerTarget, oneFilePerScoreTerm, title.prefix, plotname) {
    par(ask=F)
    
    determine_plots_per_page( pdbs, score.list, oneFilePerTarget, oneFilePerScoreTerm )  
    page.title = paste(title.prefix)

    if(!oneFilePerTarget){
        postscript(plotname,horizontal=T,paper="letter")
    } else { first = T }

    
    for (pdb in pdbs) {
        if ( length( pdbs ) > 1 ){
            plotname = paste(pdb,".ps",sep="")
        }
        if(oneFilePerTarget){
            if ( first == F ){
                mtext(page.title,line=1,side=3,outer=T,cex=1.5)
            } else { first = F }
            dev.off()
            postscript(plotname,horizontal=T,paper="letter")
        }

        plot.title = paste(title.prefix,pdb)
        plot.target(dataset[[pdb]], score.list, x_name, oneFilePerScoreTerm, plot.title)
    }
    
    mtext(plot.title,line=1,side=3,outer=T,cex=1.5)
    dev.off()
}

#--------------------------------------------------------------
# plot a score for multiple targets in one file

plot.score_type = function(dataset, pdbs, score.list, x_name, title.prefix, plotname) {
    par(ask=F)
    determine_plots_per_page( pdbs, score.list, F, T )
    for (score.i in score.list) {
        if ( length( score.list ) > 1 ){
            plotname = paste(score.i,".ps",sep="")
        }
        postscript(plotname,horizontal=T,paper="letter")
        for (pdb in pdbs) {
            x = dataset[[pdb]][x_name]
            y = dataset[[pdb]][score.i]
            d = dataset[[pdb]]$description
		    #d = dataset[[pdb]]$type
            decoyplot(x[[1]],y[[1]],d,main=score.i,xlab=x_name,ylab=score.i)
            
        }
        dev.off()
    }
}

#--------------------------------------------------------------
# function to set up the plot matrix

determine_plots_per_page = function(pdbs, score.list, oneFilePerTarget, oneFilePerScoreTerm ) {
    # initialize number of rows and columns to 1
    num_rows = 1
    num_cols = 1
    
    num_scores = length( score.list )
    num_targets = length( pdbs )
    if ( oneFilePerTarget ){
        if ( !oneFilePerScoreTerm ){
            # Use num_scores to set num_rows and num_cols
            
            # Quick hack to get a desired functionality
            if ( num_scores > 1 ){
                num_rows = 3
                num_cols = 5 
            }
        }
    } else{
        # we need to know how many targets there are
        
        if ( oneFilePerScoreTerm ){
            # Use num_targets to set num_rows and num_cols
        } else{
            # Use num_scores and num_targets to set num_rows and num_cols
        }
    }
    par(mfrow = c(num_rows, num_cols),oma=c(0,0,5,0))
    
}
    
#--------------------------------------------------------------
# wrapper function to call the correct internal function

plot.decoy_data = function(dataset,pdb.list, score.list, x_name, oneFilePerTarget=T, oneFilePerScoreTerm=F, title.prefix="", plotname="scores.ps") {
    if ( !oneFilePerTarget && oneFilePerScoreTerm ){
        plot.score_type( dataset, pdb.list, score.list, x_name, title.prefix, plotname )
    }
    else{
        plot.targets( dataset, pdb.list, score.list, x_name, oneFilePerTarget, oneFilePerScoreTerm, title.prefix, plotname )
    }   
}
    
    
    
