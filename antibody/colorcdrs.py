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
# Color mapping
colors = [Color('green', 'antigen'),
          Color('sand', 'light'),
          Color('lightteal', 'heavy'),
          Color('magenta', 'epitope'),
          Color('yellow', 'L1'),
          Color('orange', 'L2'),
          Color('salmon', 'L3'),
          Color('yellow', 'L4'),
          Color('cyan', 'H1'),
          Color('slate', 'H2'),
          Color('magenta', 'H3'),
          Color('cyan', 'H4'),
          ]

# Custom type used to simplify calls to the PyMOL API
Selection = namedtuple('Selection', ['name', 'selection'])


def _get_selections(num_scheme, dis=5.0):
    """Function to get a list of Selections and GroupedSelections for the
    antibody/antigen complex.

    :param num_scheme: str
    :param dis: float
    :rtype: (list, dict)
    """

    # helper functions to make string substitutions more straightforward
    def chain(cdr_name):
        return 'light' if cdr_name.startswith('L') else 'heavy'

    def stem_res_no(start_stop):
        return [start_stop[0] - 2, start_stop[0] - 1, start_stop[1] + 1,
                start_stop[1] + 2]

    sel = []
    grp = {}

    # general selections
    sel.append(Selection('bb', 'name CA+C+O or (not resn PRO and name N)'))
    sel.append(Selection('bbn', '(bb and not name CA) or elem h*'))
    sel.append(Selection('water', 'resname HOH'))

    # antibody-specific selections
    sel.append(Selection('ab', 'chain L or chain H'))
    sel.append(Selection('antigen', 'not ab'))

    for letter in ['H', 'L']:
        c = chain(letter)
        sel.append(Selection(c, 'chain {}'.format(letter)))
        grp['{}_group'.format(c)] = [c]

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
    around = 'byres ({} around {}) and not {}'
    for cdr, range in sorted(scheme.iteritems()):
        # select the CDR loop
        sel.append(Selection(cdr, '{} and resi {}-{}'.format(chain(cdr),
                                                             *range)))

        # CDR Epitope
        sel.append(Selection('{}_epitope'.format(cdr),
                             around.format(cdr, dis, 'ab')))

        # CDR Frame
        sel.append(Selection('{}_frame'.format(cdr),
                             around.format(cdr, dis, 'antigen')))

        # select the CDR loop stem residues
        sel.append(Selection('{}_stem'.format(cdr), '{} and (resi {}-{} or '
                             'resi {}-{})'.format(chain(cdr),
                                                  *stem_res_no(range))))

        grp['{}_group'.format(cdr)] = [cdr,
                                       '{}_epitope'.format(cdr),
                                       '{}_frame'.format(cdr),
                                       '{}_stem'.format(cdr)
                                       ]
    # Chain group with CDRs
    r = ' or '
    s = scheme.keys()
    for letter in ['H', 'L']:
        cdrs = [x for x in s if x.startswith(letter)]
        c = chain(letter)

        sel.append(Selection('{}_cdrs'.format(c), r.join(cdrs)))

        sel.append(Selection('{}_framework'.format(c),
                             '{0} and not {0}_cdrs'.format(c)))

        sel.append(Selection('{}_stem'.format(c),
                             r.join([x + '_stem' for x in cdrs])))

        grp['{}_group'.format(c)].extend(['{}_cdrs'.format(c),
                                           '{}_framework'.format(c),
                                           '{}_stem'.format(c)
                                           ])

    sel.append(Selection('cdrs', r.join(s)))
    sel.append(Selection('epitope', around.format('ab', dis, 'ab')))
    sel.append(Selection('epitope_cdrs', r.join([x + '_epitope' for x in s])))
    sel.append(Selection('paratope', around.format('epitope', dis, 'antigen')))

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
        sel.append(Selection('kinkplus', 'heavy and resi 94+100A-103'))
        sel.append(Selection('kinkhbond', 'heavy and resi 94+101'))

        grp['H3_kink'] = ['kinkplus', 'kinkhbond']

    return sel, grp


def colorcdrs(numbering='chothia', cartoon=False, neighbor_distance=5.0):
    """
DESCRIPTION

    "colorcdrs" creates named atom selections for framework regions, CDR
    loops based on a canonical numbering scheme and colors them.

USAGE

    colorcdrs [numbering_scheme [, cartoon [, neighbor_distance ]]]

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

    selections, groups = _get_selections(numbering, neighbor_distance)

    for s in selections:
        cmd.select(s.name, s.selection)

    for group, names in sorted(groups.items(), key=lambda x: x[0].lower()):
        for name in names:
            cmd.group(group, name)

    sele = '(all)'

    if cartoon:
        cmd.hide(representation='lines', selection=sele)
        cmd.show(representation='cartoon', selection=sele)
        cmd.cartoon(type='loop', selection=sele)

    # cmd.bg_color(color='white')

    for c in colors:
        cmd.color(c.color, c.selection)

    # cmd.disable(name='framework')
    # cmd.disable(name='stem')

    # cmd.hide(representation='everything', selection='not ab')
    cmd.zoom(selection='ab')

cmd.extend('colorcdrs', colorcdrs)
