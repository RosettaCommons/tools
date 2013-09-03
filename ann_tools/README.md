#ANN Tools


These scripts are designed to aid in the training, optimization and analysis of Rosetta vHTS neural nets. 

##Prerequisites

* mol files produced by the Rosetta vHTS ligand docking protocol.
   * For the fastest performance, 
* A copy of the BCL (http://bclcommons.vueinnovations.com/)
* python

##BCL data input caveats

The BCL only supports V2000 SDF files as input.  It handles files that are uncompressed, or compressed with gzip or bzip2.  It is extremely picky about the SDF file specification and will not parse files which deviate even slightly from the published specification.  Files produced by Rosetta tightly adhere to this specification and will be parsed with no additional processing.  Files from other sources (including ChEMBL and the PDB) may require additional postprocessing.  Running an input file through babel is genearl 

##Scripts and Files

* default_features.object
   * A file describing all the features produced by the BCL and computable by Rosetta.  This file will be used as the starting point for network training and optimization
   * features which start "MiscProperty" are computed by Rosetta and must be found in every record in the input molfile. All other features are computed internally by the BCL
   
* sort_sdfs.py
   * sort a directory of sdfs by the Rosetta interface score (interface_delta_X)
   * both directories must exist
   * records are sorted in ascending order by score. 
   * To run this script:
   
   ```
   sort_sdfs.py input_dir/ output_dir/
   ```
