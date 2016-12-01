#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# @file   colorcdrs.py
# @brief  This script can be run within the PyMOL python interpreter to provide
#         the 'colorcdrs' function. The intended use is to run this script
#         from a user's ~/.pymolrc
# @author Brian D. Weitzner (brian.weitzner@gmail.com)
# @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from __future__ import print_function

from collections import namedtuple
from pymol import cmd


# different numbering schemes should be added here.
schemes = {'chothia': {'L1': [24, 34], 'L2': [50, 56], 'L3': [89, 97],
                       'L4': [66, 71],
                       'H1': [26, 35], 'H2': [50, 56], 'H3': [95, 102],
                       'H4': [71, 78]
                       },
           'aho': {'L1': [24, 42], 'L2': [57, 72], 'L3': [107, 138],
                   'L4': [82, 89],
                   'H1': [24, 42], 'H2': [57, 69], 'H3': [107, 138],
                   'H4': [82, 89]
                   },
           }

Color = namedtuple('Color', ['color', 'selection'])
Selection = namedtuple('Selection', ['name', 'selection'])
GroupedSelection = namedtuple('GroupedSelection', ['group', 'name',
                                                   'selection'])


# Color mapping

colors_a = [

    Color('green', 'antigen'),
    Color('sand', 'light'),
    Color('lightteal', 'heavy'),
    Color('magenta', 'epitope'),

    Color('yellow',  'L1'),
    Color('orange',  'L2'),
    Color('salmon',  'L3'),
    Color('yellow',  'L4'),

    Color('cyan',    'H1'),
    Color('slate',   'H2'),
    Color('magenta', 'H3'),
    Color('cyan',    'H4'),

        ]

colors_b = [Color('green', 'antigen'),
            Color('yellow', 'light'),
            Color('blue', 'heavy'),
            Color('brightorange', 'cdrL'),
            Color('cyan', 'cdrH'),
            Color('greencyan', 'H3'),
            Color('orange', 'L3'),
            ]


def _get_selections(num_scheme, neighbor_dis=5.0):
    """Function to get a list of Selections and GroupedSelections for the
    antibody/antigen complex.

    :param num_scheme: dict
    :param neighbor_dis: int
    :rtype: list
    """

    # helper functions to make string substitutions more straightforward
    def chain(cdr_name):
        return 'light' if cdr_name.startswith('L') else 'heavy'

    def stem_res_no(start_stop):
        return [start_stop[0] - 2, start_stop[0] - 1, start_stop[1] + 1,
                start_stop[1] + 2]

    sel = []

    # general selections
    sel.append(Selection('bb', 'name CA+C+O or (not resn PRO and name N)'))
    sel.append(Selection('bbn', '(bb and not name CA) or elem h*'))
    sel.append(Selection('water', 'resname HOH'))

    # antibody-specific selections
    sel.append(Selection('ab', 'chain L or chain H'))
    sel.append(Selection('antigen', 'not ab'))

    for letter in ['H', 'L']:
        sel.append(GroupedSelection('{}_group'.format(chain(letter)),
                                    chain(letter), 'chain {}'.format(letter)))

    # cdr selections
    # default to chothia if necessary
    if num_scheme not in schemes.keys():
        print('{} is not a valid numbering scheme. '.format(num_scheme) +
              'Valid schems are: {}. '.format(schemes.keys()) +
              'Defaulting to "chothia".')
        num_scheme = 'chothia'

    scheme = schemes[num_scheme]

    # The actual selections happen here.
    # The order of selections matter, as we select based on previous selections
    for cdr, range in sorted(scheme.iteritems()):
        # select the CDR loop
        print(cdr)
        sel.append(GroupedSelection('{}_group'.format(cdr), cdr,
                                    '{} and resi {}-{}'.format(chain(cdr),
                                                               *range)))

        # CDR Epitope
        sel.append(GroupedSelection('{}_group'.format(cdr),
                                    '{}_epitope'.format(cdr),
                                    'br. ({} around {}) and not '
                                    'ab'.format(cdr, neighbor_dis)))

        # CDR Frame
        sel.append(GroupedSelection('{}_group'.format(cdr),
                                    '{}_frame'.format(cdr),
                                    'br. ({} around {}) and not '
                                    'antigen'.format(cdr, neighbor_dis)))

        # select the CDR loop stem residues
        sel.append(GroupedSelection('{}_group'.format(cdr),
                                    '{}_stem'.format(cdr),
                                    '{} and (resi {}-{} or resi '
                                    '{}-{})'.format(chain(cdr),
                                                    *stem_res_no(range))))

    # Chain group with CDRs
    for letter in ['H', 'L']:
        cdrs = [x for x in scheme.keys() if x.startswith(letter)]

        sel.append(GroupedSelection('{}_group'.format(chain(letter)),
                                    '{}_cdrs'.format(chain(letter)),
                                    '{} or {} or {}'.format(*cdrs)))

        sel.append(GroupedSelection('{}_group'.format(chain(letter)),
                                    '{}_framework'.format(chain(letter)),
                                    '{} and not '
                                    '{}_cdrs'.format(chain(letter),
                                                     chain(letter))))

        sel.append(GroupedSelection('{}_group'.format(chain(letter)),
                                    '{}_stem'.format(chain(letter)),
                                    '{}_stem or {}_stem or '
                                    '{}_stem'.format(*cdrs)))

    sel.append(Selection('cdrs', ' or '.join(scheme.keys())))

    sel.append(Selection('epitope', 'br. (ab around {}) and not '
                         'ab'.format(neighbor_dis)))

    sel.append(Selection('epitope_cdrs', ' or '.join(
        ['{}_epitope'.format(cdr) for cdr in scheme.keys()])))

    sel.append(Selection('paratope', 'br. (epitope around {}) and '
                         'not antigen'.format(neighbor_dis)))

    # Alternate selection strings for heavy and light frameworks are below.
    # I am not entirely sure how these differ from those above and what
    # advantage, if any, they provide
    # select lightframework, light and (resi 4-6 or resi 10-23 or resi \
    # 35-38 or resi 45-49 or resi 57-66 or resi 71-88 or resi 98-105)
    # select heavyframework, heavy and (resi 4-6 or resi 10-25 or resi \
    # 36-39 or resi 46-49 or resi 66-94 or resi 103-110)

    sel.append(Selection('framework', 'light_framework or heavy_framework'))
    sel.append(Selection('stem', 'light_stem or heavy_stem'))

    # TODO: right now this maps to chothia numbering only. Fix that.
    if num_scheme == 'chothia':
        sel.append(GroupedSelection('H3_kink', 'kinkplus',
                                    'heavy and resi 94+100A-103'))
        sel.append(GroupedSelection('H3_kink', 'kinkhbond',
                                    'heavy and resi 94+101'))

    return sel


def colorcdrs(numbering='chothia', group=True, neighbor_distance=5.0,
              classic_coloring=False):
    """
DESCRIPTION

    "colorcdrs" creates named atom selections for framework regions, CDR
    loops based on a canonical numbering scheme and colors them.

USAGE

    colorcdrs [numbering_scheme [, group [, neighbor_distance
              [, classic_coloring ]]]]

ARGUMENTS

    numbering_scheme = the name of the numbering scheme used by the antibody.
                       valid options are 'chothia' and 'aho'.

    group = Should we group the selections?

    neighbor_distance = Distance to measure neighbor/interface residues

    classic_coloring = Use the classic coloring scheme (colors less elements)

NOTES

    If no parameters are passed to the "colorcdrs" command, the chothia
    numbering scheme will be used.

EXAMPLES

    colorcdrs
    colorcdrs aho

PYMOL API

    cmd.colorcdrs(string numbering='chothia', bool group=True,
                 float neighbor_distance=5.0, bool classic_coloring=False)

    """

    selections = _get_selections(numbering, neighbor_distance)

    for s in selections:
        cmd.select(s.name, s.selection)

    if group:
        for s in selections:
            if type(s) == GroupedSelection:
                cmd.group(s.group, s.name)

    sele = '(all)'

    # cmd.hide(representation='lines', selection=sele)
    # cmd.show(representation='cartoon', selection=sele)
    # cmd.cartoon(type='loop', selection=sele)

    if classic_coloring:
        colors = colors_b
    else:
        colors = colors_a

    for c in colors:
        cmd.color(c.color, c.selection)

    # cmd.disable(name='framework')
    # cmd.disable(name='stem')

    # cmd.hide(representation='everything', selection='not ab')
    cmd.center(selection='epitope paratope')
    # cmd.select(name='doc', selection='resname DOC')

cmd.extend('colorcdrs', colorcdrs)
