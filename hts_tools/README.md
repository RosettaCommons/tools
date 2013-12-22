HTS Tools
=========

These scripts are designed to aid in the setup and running of vHTS jobs in Rosetta.  They form the basis of a pipeline for processing input SDF files into params files suitable for input into Rosetta.

These tools are almost certainly overkill if you have < 100 ligands to process. They were written with >1 million compound libraries in mind. 

Prerequisites
-------------

* python
* a copy of ```molfile_to_params.py```
* a directory of sdf files
   * one file for each ligand and its conformers
   * each conformer must have the ligand name in the 1st row of the sdf record
   * each conformer must have hydrogens and 3D coordinates
* a list file
   * The list file maps filenames to ligands.  It is a CSV format file that is in a format like this:
   
   ```
      ligand_id,filename
      string,string
      ligand_1,path/to/ligand1
      ligand_2,path/to/ligand2
   ```
   * The filenames in this file should be relative to the working directory.  using absolute paths is not recommended
* CORINA (optional)
   * if you want to add hydrogens using the included script you will need a copy of corina
   
Scripts
-------

* sdf_split_organize.py 
   * This is a helper script for producing a directory of sdf files.  It takes as input a single sdf file, and will split that file into multiple files, each file containing all the conformers for 1 ligand.  Different ligands must have different names in the sdf records, and all conformers for one ligand must have the same name.  output filenames are based on the sha1 hash of the input filename, and are placed in a directory hashed structure. Thus, a ligand with the name "Written by BCL::WriteToMDL,CHEMBL29197" will be placed in the path ```/41/412d1d751ff3d83acf0734a2c870faaa77c28c6c.mol```.  The script will also output a list file in the format described in the prerequisites section above. 
   * To run the script:
   
   ```
      sdf_split_organize.py input.sdf output_dir/ output_file_list.csv
   ```
* setup_screening_project.py
   * This is the first script you use, it takes as input a csv list file and outputs an sqlite3 database that will act as both the input and output for the remaining scripts
   * by default this script will check that every file in the input list exists.  This can be very time consuming and can be disabled with ```--no-verify``` 
   * to run the script:
   
   ```
      setup_screening_project.py input.csv output.db3
   ``` 
* add_activity_tags_to_database.py
   * This script will add activity data to the database. It requires a tag file, in the following format:
   
   ```
      ligand_id,tag,value
      string,string,float
      ligand_1,foo,1.5
      ligand_2,bar,-3.7
   ```
   * ligand IDs should match the ligand IDs in the database. only ligand IDs that exist in both the database and the tag file will be processed.
   * to run the script:
   
   ```
      add_activity_tags_to_database.py database.db3 tag_file.csv
   ```
* add_hydrogens.py
   * This script will add hydrogens to files listed in the database using corina.  It will skip files that have already been processed by corina, and does not modify the 3D structure
   * By default the script will process all files. --only_tagged will cause the script to only process files with activity data
   * The script is multithreaded.  use -j to specify the number of threads you want to use.
   * molfiles must be decompressed prior to use
   * to run the script:
   
   ```
      add_hydrogens.py -j 4 database.db3
   ```
   
* make_params.py
   * This script will make param files for every ligand with activity tag data.
   * The script is multthreaded. use -j to specify the number of threads you want to use. 
   * The script splits param, conformer and pdb files into 255 directories by the first 2 characters of the sha1 of the ligand name
   * The length of the ligand name is determined by the number of input ligands.  if more than 46656 ligands are processed, ligand names will be > 3 characters
   * Ligands that cannot be processed by molfile_to_params are skipped.
   * the tag name and tag value are written as string and numeric properties respectively to the params file
   * the database is updated with the paths to the newly created ligands.  If one ligand has more than one activity tag, two params files will be generated and referenced in the database.  This allows for the possibility of a ligand having measured activity in more than one system.
   * to run the script:
   
   ```
      make_params.py -j 4 --database path/to/Rosetta/main/database --path_to_params path/to/molfile_to_params.py database.db3 output_dir/
   ```
   
* add_ligands_to_job_file.py
   * This script will add ligands with existing params files to a job file formatted for use with -in:file:screening_job_input
   * The input job file must have a "group_name" block for each record, the group_names must match the tags in the database
   * To run the script:
   
   ```
      add_ligands_to_job_file.py database.db3 input_screening_file.js output_screening_file.js
   ```

* get_descriptor_data.py
   * This script gets descriptor information from data output to a database using the screening features reporter
   * This script requires that sqlalchemy be installed
   * Currently this script only makes mysql engine connections, but can be easily modified to support other backends
   * To run the script:
   
   ```
      get_descriptor_data.py --batch-description=batch_description --host=host.com --username=user --database-name=db_name descriptor_data.js
   ```

* prepare_sdfs_for_bcl.py
   * given a ligand database and a json file produced by get_descriptor_data.py, output an sdf file with every field in the descriptor json file as a property
   * if multiple conformers are generated, only 1 conformer for each ligand will be used.
   * to run the script:
   
   ```
      prepare_sdfs_for_bcl.py database.db3 descriptor_data.js ligands_for_bcl.sdf
   ```
   
* clean_pdb.py
   * remove all residues except those specified and the 20 canonical AAs.
   * compute the center of a ligand (even if it is being removed) and add it to a JSON map
   * to run the script:
   
   ```
      clean_pdb.py --chains=A,B --ligands_to_keep=CL,ABC --ligand_for_center=LIG input.pdb output.pdb center_file.js
   ```
   
* make_evenly_grouped_jobs.py
   * Given a directory of params files prepared with make_params.py, and a directory of pdb files, produce n job files balanced so that each group of jobs contains no more than x jobs each.  Each job file is constructed to be as evenly sized as possible. n and x should be tuned to the memory and walltime requirements of your cluster. 
   * Each params file must have a "system_name" tag, and this tag must be a substring of the corresponding pdb files
      * For example, if the "system_name" tag is "1UBI", the script will look for "1UBI*.pdb" in the pdb directory.
   * The params files necessary to run the job file are included in the job file, removing the need for -in:file:extra_res_fa usage
   * job files can be loaded into rsetta with -in:file:screening_job_inputter
   * files output will be in the form output_prefix_n.js
   
   ```
   make_evenly_grouped_jobs.py params_dir/ pdb_dir/ output_prefix
   ```
   
* make_startfrom_files.py
   * Given a directory of input pdb files and a center file produced by clean_pdb.py (above), make a startfrom file for the Rosetta StartFrom mover

   ```
      make_startfrom_files.py center_file.js input_models/ start_from.js
   ```