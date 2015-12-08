fa_score.list = list( "total_score", "dslf_ca_dih", "dslf_cs_ang",
    "dslf_ss_dih", "dslf_ss_dst", "fa_atr",	"fa_dun", "fa_intra_rep", 
    "fa_pair", "fa_rep", "fa_sol", "hbond_bb_sc", "hbond_lr_bb", "hbond_sc",	
    "hbond_sr_bb", "omega", "p_aa_pp", "pro_close", "rama", "ref" )

cen_score.list = list("score", "d.env", "d.pair", "contact", "d.vdw", "d.fab",
	"env", "pair", "vdw", "hs", "ss", "sheet", "cb", "rsigma", "hb", "rg")

# attempt to set a default set of scores to plot
if (fullatom){
    score.list = fa_score.list 
}else{
    score.list = cen_score.list
}

# attempt to detect all scorefiles to collect data from
pdbs = list()
if (fullatom){ 
    extension = "*.fasc"
} else{ 
    extension = "*.sc" 
}
command = paste( "ls", extension, "| awk '{split($0, a, \".\"); print a[1]}'" )
pdbs = as.list( system(command, T))

#--------------------------------------------------------------
# load data in score files for each pdb of interest

decoy_data.load = function( pdbs ) {
    all.data = pdbs
    names(all.data) = pdbs
    if (fullatom) extension = ".fasc" else extension = ".sc"

    for (pdb in pdbs) {
        filename = paste(pdb,extension,sep="")
        print(paste("loading ",filename))
        all.data[[pdb]] = read.table(filename,header=T, , skip=1) 
    }
    return( all.data )
}

