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

def create_blast_db():
    """
    Using info files, extract sequences per length and setup blast databases.
    """
    return

def extract_info():
    """
    Several quality filters need info files. 
    Generate these from the structures + SAbDab summary file.
    """
    return

def parse_sabdab_summary(file_name):
    """
    SAbDab produces a rather unique summary file. 
    This function reads that file into a dict with the key being the
    4-letter PDB code.
    """
    # dict[pdb] = {col1 : value, col2: value, ...}
    sabdab_dict = {}

    with open(file_name, "r") as f:
        # first line is the header, or all the keys in our sub-dict
        header = f.readline().strip().split("\t")
        # next lines are data
        for line in f.readlines():
            split_line = line.strip().split("\t")
            td = {} # temporary dict of key value pairs for one pdb
            for k,v in zip(header[1:], split_line[1:]):
                # pdb id is first, so we skip that for now
                td[k] = v
            # add temporary dict to sabdab dict at the pdb id
            sabdab_dict[split_line[0]] = td

    return sabdab_dict

def truncate_antibody_pdbs():
    """
    We only use the Fv as a template, so this function loads each pdb
    and deletes excess chains/residues. We define the Fv, under the
    Chothia numbering scheme as H1-H112 and L1-L109.
    """
    # read SAbDab info for chain identities
    sabdab_dict = parse_sabdab_summary("info/sabdab_summary.tsv")
    # iterate over PDBs in antibody_database and truncate accordingly
    unique_pdbs = set([x[:4] for x in os.listdir("antibody_database") if x.endswith(".pdb") or x.endswith(".pdb.bz2")])
    for pdb in unique_pdbs:
        # try reading bzipped pdb, then regular pdb
        pdb_text = ""
        try:
            with open("antibody_database/" + pdb + ".pdb.bz2", "r") as f:
                pdb_text = bz2.decompress(f.readlines())
        except IOError:
            try:
                with open("antibody_database/" + pdb + ".pdb", "r") as f:
                    pdb_text = f.readlines()
            except IOError:
                sys.exit("Failed to open {} in antibody_database/ !".format(pdb))
        # should have pdb_text now
        if len(pdb_text) == 0: sys.exit("Nothing parsed for PDB {} !".format(pdb))
        # great success lets extract the correct chains, sometimes "NA"
        hchain = sabdab_dict[pdb]["Hchain"]
        lchain = sabdab_dict[pdb]["Lchain"]
        # keep writing here


    # chains are at position 21 in PDBs
    print unique_pdbs
    # resnums are at 22-26
    return

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

    # make info and antibody_database directories if not present
    if not os.path.isdir("info"): os.mkdir("info")
    if not os.path.isdir("antibody_database"): os.mkdir("antibody_database")

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
    #create_antibody_db()
    truncate_antibody_pdbs()

