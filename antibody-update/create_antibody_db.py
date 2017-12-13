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
## Some directories (antibody_database, info, etc...) are hard coded.

import urllib2
import re
import os
import bz2

def download_antibody_pdbs():
    """
    Function queries SAbDab for the latest list of non-redundant,
    chothia-numbered, sub 3-Angstrom antibody structures.
    """
    # base url for SAbDab query
    url_base = 'http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/'
    url_query_base = url_base + 'DBrowser.php?'

    # configure options for SAbDab query
    options = { 'submitted' : 'true', 
                'nonredundant' : 'true',
                'ab_ident' : '99%25',# 25 need for html
                'only_complete' : 'false', # might need capital F
                'incomple_nr' : '*',
                'agtype' : '*',
                'ag_ident' : '99%25',
                'aglengthnr' : '50',
                'resolution' : '3.0',
                'nonredundant_selected' : 'true' }

    # generate url from query
    sabdab_url =  url_query_base + '&'.join(['{0}={1}'.format(x,y) for x,y in options.items()])

    # open page
    webpage = urllib2.urlopen(urllib2.Request(sabdab_url))
    
    # read html and extract summary file, and download chothia PDBs
    webpage_html = webpage.read() # PDBs in table, summary at end

    # use regex to extract PDBs
    chothia_pdb_re = re.compile('data/entries/\w\w\w\w/structure/chothia/\w\w\w\w.pdb')
    summary_re = re.compile('a href=.* download=\"summary.csv\"')
    pdb_urls = chothia_pdb_re.findall(webpage_html, re.MULTILINE)
    summary_url = summary_re.findall(webpage_html, re.MULTILINE) # should only find two
    summary_url = summary_url[0].split(" download")[0][8:-1] # drop unnecessary chars

    # download summary file, always overwrite
    with open("info/sabdab_summary.tsv", "w") as f:
        u = urllib2.urlopen(urllib2.Request(url_base + summary_url))
        f.write(u.read())

    # download PDBs, check for *.pdb or *.pdb.bz before writing
    for pdb in pdb_urls:
        fpath = "antibody_database/{}.pdb".format(pdb[-8:-4])
        fpath_bz = "antibody_database/{}.pdb.bz2".format(pdb[-8:-4])
        if not (os.path.isfile(fpath) or os.path.isfile(fpath_bz)):
            with open(fpath_bz, "w") as f:
                print "Downloading " +  pdb[-8:-4] + "..."
                u = urllib2.urlopen(urllib2.Request(url_base + pdb))
                f.write(bz2.compress(u.read()))
    
    return

def create_antibody_db():
    download_antibody_pdbs()
    return

if __name__ == "__main__":
    # check for execution in correct dir
    # or risk creation of dirs/files elsewhere
    create_antibody_db()

