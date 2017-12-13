#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   create_antibody_db.py
## @brief  Script to create database for RosettaAntibody
## @author Jeliazko Jeliazkov

## @details download non-redundant Chothia Abs from SAbDab
## Abs are downloaded by html query (is there a better practice)?
## Abs are Chothia-numbered, though we use Kabat to define CDRs.
## After download, trim Abs to Fv and extract FR and CDR sequences.
## For the purpose of trimming, truncate the heavy @112 and light @109

def download_antibody_pdbs():
    """
    Function queries SAbDab for the latest list of non-redundant,
    chothia-numbered, sub 3-Angstrom antibody structures.
    
    http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/DBrowser.php?submitted=true&nonredundant=true&ab_ident=99%25&only_complete=False&incomplex_nr=*&agtype=either&ag_ident=99%25&aglengthnr=50&resolution=3.0&submitted=true&nonredundant_selected=true
    """
    # configure options for SAbDab query
    options = { 'submitted' = 'true', 
                'nonredundant' = 'true',
                'ab_ident' = '99%25',# 25 need for html
                'only_complete' = 'false', # might need capital F
                'incomple_nr' = '*',
                'agtype' = '*',
                'ag_ident' = '99%25',
                'aglengthnr' = '50',
                'resolution' = '3.0',
                'nonredundant_selected' = 'true' }

def create_antibody_db():
    return

if __name__ == "__main__":
    create_antibody_db()

