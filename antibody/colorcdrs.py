#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# @file   colorcdr.py
# @brief  This script can be run within the PyMOL python interpreter to provide
#         the 'color_cdr' function. The intended use is to run this script
#         from a user's ~/.pymolrc
# @author Brian D. Weitzner (brian.weitzner@gmail.com)

from __future__ import print_function

from collections import namedtuple
from pymol import cmd


# different numbering schemes should be added here.
schemes = {'chothia': {'l1': [24, 34], 'l2': [50, 56], 'l3': [89, 97],
                       'h1': [26, 35], 'h2': [50, 56], 'h3': [95, 102]
                       },
           'aho': {'l1': [24, 42], 'l2': [57, 72], 'l3': [107, 138],
                   'h1': [24, 42], 'h2': [57, 69], 'h3': [107, 138]},
           }

Selection = namedtuple('Selection', ['name', 'selection'])


def _get_selections(num_scheme):
    '''I should totally put a useful docstring here'''
    sel = []

    # general selections
    sel.append(Selection('bb', 'name CA+C+O or (not resn PRO and name N)'))
    sel.append(Selection('bbn', '(bb and not name CA) or elem h*'))
    sel.append(Selection('water', 'resname HOH'))

    # antibody-specific selections
    sel.append(Selection('light', 'chain L'))
    sel.append(Selection('heavy', 'chain H'))

    sel.append(Selection('ab', 'light or heavy'))
    sel.append(Selection('antigen', 'not ab'))

    # cdr selections
    # default to chothia if necessary
    if num_scheme not in schemes.keys():
        print('{} is not a valid numbering scheme. '.format(num_scheme) +
              'Valid schems are: {}. '.format(schemes.keys()) +
              'Defaulting to "chothia".')
        num_scheme = 'chothia'

    scheme = schemes[num_scheme]

    # helper functions to make string substitutions more straightforward
    def chain(cdr_name):
        return 'light' if cdr_name.startswith('l') else 'heavy'

    def stem_res_no(start_stop):
        return [start_stop[0] - 2, start_stop[0] - 1, start_stop[1] + 1,
                start_stop[1] + 2]

    # the actual selections happen here
    for cdr, range in sorted(scheme.iteritems()):
        # select the CDR loop
        sel.append(Selection(cdr, '{} and resi {}-{}'.format(chain(cdr),
                                                             *range)))
        # select the CDR loop stem residues
        sel.append(Selection('{}s'.format(cdr),
                             '{} and (resi {}-{} or resi {}-{})'.format(
                                 chain(cdr), *stem_res_no(range))))

    for letter in ['h', 'l']:
        cdrs = [x for x in scheme.keys() if x.startswith(letter)]

        sel.append(Selection('cdr{}'.format(letter),
                             '{} or {} or {}'.format(*cdrs)))

        sel.append(Selection('{}framework'.format(chain(letter)),
                             '{} and not cdr{}'.format(chain(letter), letter)))

        sel.append(Selection('{}s'.format(letter),
                             '{}s or {}s or {}s'.format(*cdrs)))

    sel.append(Selection('cdr', 'cdrl or cdrh'))

    # Alternate selection strings for heavy and light frameworks are below.
    # I am not entirely sure how these differ from those above and what
    # advantage, if any, they provide
    # select lightframework, light and (resi 4-6 or resi 10-23 or resi \
    # 35-38 or resi 45-49 or resi 57-66 or resi 71-88 or resi 98-105)
    # select heavyframework, heavy and (resi 4-6 or resi 10-25 or resi \
    # 36-39 or resi 46-49 or resi 66-94 or resi 103-110)

    sel.append(Selection('framework', 'lightframework or heavyframework'))
    sel.append(Selection('stem', 'ls or hs'))

    # TODO: right now this maps to chothia numbering only. Fix that.
    if num_scheme == 'chothia':
        sel.append(Selection('kinkplus', 'heavy and resi 94+100A-103'))
        sel.append(Selection('kinkhbond', 'heavy and resi 94+101'))

    return sel

def colorcdrs(numbering='chothia', paratope=False, epitope=False, group=False):
    '''
DESCRIPTION

    "colorcdrs" creates named atom selections for framework regions, CDR
    loops based on a canonical numbering scheme and colors them.

USAGE

    colorcdrs [numbering_scheme [, paratope [, epitope [, group ]]]]

ARGUMENTS

    numbering_scheme = the name of the numbering scheme used by the antibody.
                       valid options are 'chothia' and 'aho'.

    paratope = boolean?

    epitope = boolean?

    group = boolean?

NOTES

    If no parameters are passed to the "colorcdrs" command, the chothia
    numbering scheme will be used.

EXAMPLES

    colorcdrs
    colorcdrs aho

PYMOL API

    cmd.colorcdrs(string numbering='chothia', bool paratope=False,
                  bool epitope=False, bool group=False)

    '''

    selections = _get_selections(numbering)

    for s in selections:
        cmd.select(s.name, s.selection)

    sele = '(all)'
    cmd.hide(representation='lines', selection=sele)
    cmd.show(representation='cartoon', selection=sele)
    cmd.cartoon(type='loop', selection=sele)

    cmd.bg_color(color='white')

    cmd.color('green', selection='antigen')
    cmd.color('yellow', selection='light')
    cmd.color('blue', selection='heavy')
    cmd.color('brightorange', selection='cdrl')
    cmd.color('cyan', selection='cdrh')
    cmd.color('greencyan', selection='h3')
    cmd.color('orange', selection='l3')

    cmd.disable(name='framework')
    cmd.disable(name='stem')

    cmd.zoom(selection=sele)
    cmd.hide(representation='everything', selection='not ab')
    cmd.zoom(selection='ab')
    # cmd.select(name='doc', selection='resname DOC')

cmd.extend('colorcdrs', colorcdrs)
