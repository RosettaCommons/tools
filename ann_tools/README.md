#ANN Tools


These scripts are designed to aid in the training, optimization and analysis of Rosetta vHTS neural nets. 

##Prerequisites

* mol files produced by the Rosetta vHTS ligand docking protocol
* A copy of the BCL (http://bclcommons.vueinnovations.com/)
* python

##Scripts and Files

* default_features.object
   * A file describing all the features produced by the BCL and computable by Rosetta.  This file will be used as the starting point for network training and optimization
   * features which start "MiscProperty" are computed by Rosetta and must be found in every record in the input molfile. All other features are computed internally by the BCL
