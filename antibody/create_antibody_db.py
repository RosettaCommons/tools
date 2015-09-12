#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   antibody.py
## @brief  Script to create antibody.db file from old info/* files
## @author Sergey Lyskov

import sys


def read_formated_text_file(file_name):
    ''' return dict with first field as key and value equal to dict(field_name=value)
    '''
    info, legend = {}, ''
    for l in file(file_name):
        if l.startswith('# '): legend = l[2:].split()
        elif len(l)>8: info[l.split()[0]] =  dict( zip(legend, l.split() ) )

    return info, legend


def write_formated_text_file(file_name, data, keys):
    ''' Create plain text file ising data dict and fill it with keys columns
    '''

    with file(file_name, 'w') as f:

        column_lengths = {}
        for k in keys:
            l = len('# '+k) if k==keys[0] else len(k)  # adjusting for legend prefix
            for p in data:
                if k in data[p]: l = max(l, len(data[p][k]) )
            column_lengths[k] = l


        f.write(' '.join(['{0:^{1}}'.format('# '+k if k==keys[0] else k, column_lengths[k]) for k in keys]) + '\n')

        sorted_pdbs = sorted(data.keys()) #, key=lambda x: x.partition('_chothia.pdb')[0][3:] )

        for p in sorted_pdbs:
            value = data.get(k, '-')

            f.write( ' '.join([ '{0:^{1}}'.format(data[p].get(k, '-'), column_lengths[k]) for k in keys]) + '\n')





def main(args):
    ''' Script to create antibody.db file from old info/* files
    '''


    cdr_info, cdr_legend = read_formated_text_file('info/cdr_info')


    frh_info, frh_legend = read_formated_text_file('info/frh_info')
    frl_info, frl_legend = read_formated_text_file('info/frl_info')

    # altering data representation
    for p in cdr_info.keys():
        if p in frh_info: cdr_info[p]['FRH'] = frh_info[p]['FRH']
        if p in frl_info: cdr_info[p]['FRL'] = frl_info[p]['FRL']

        new_key = p.partition('_chothia.pdb')[0][3:]
        cdr_info[new_key] = cdr_info[p]
        del cdr_info[p]

        cdr_info[new_key]['pdb'] = new_key

        for k in cdr_info[new_key]:
            cdr_info[new_key][k] = '-' if cdr_info[new_key][k]=='none' else cdr_info[new_key][k]


    write_formated_text_file('info/antibody.db', cdr_info, 'pdb resolution BioType date LightType StructSouce H1 H2 H3 L1 L2 L3 FRH FRL'.split())



if __name__ == "__main__": main(sys.argv)
