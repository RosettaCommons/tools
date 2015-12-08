x_name="rms"
fullatom = T

oneFilePerTarget=T
oneFilePerScoreTerm=F
title.prefix=""

# Ignored if one file per target is being used
plotname="scores.ps"

script_path = Sys.getenv("rosetta_R_scripts")
if( nchar( script_path[[1]] ) > 0 ){
	script_path = paste( script_path[[1]], "/", sep="" )
}
 
source( paste( script_path[[1]], "decoy_data.R", sep="" ))
source( paste( script_path[[1]], "decoy_scatterplots.R", sep=""))

#score.list = list()
score.list = list( "total_score" )
print( score.list )
#pdbs = list()
#pdbs = list( "1BZQ" )

dataset = decoy_data.load( pdbs )
plot.decoy_data(dataset, pdbs, score.list, x_name, oneFilePerTarget, oneFilePerScoreTerm, title.prefix, plotname)

