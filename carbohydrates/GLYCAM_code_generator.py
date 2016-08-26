#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

"""File:  GLYCAM_code_generator.py

Brief:  This Python script auto-generates a .codes file of alternative three-
letter codes for use in Rosetta when loading in PDB files that use the GLYCAM
naming scheme.

Details:  See glycam.org/docs/forcefield/glycam-naming-2/ for more information.

To properly add to this program, be sure to add new codes to both the lists and
dictionaries for the linkage codes and sugar codes.

Author:  Morgan Nance
Labonte <JWLabonte@jhu.edu> (PEP8/style tweaks, GLYCAM code additions)

"""

#Imports
from __future__ import with_statement
import os
from carbohydrate_data import *


# filename preparation
working_dir = os.getcwd() + '/'
filename = working_dir + "glycam.codes"
file_header = "# A list of 3-letter GLYCAM codes for monosaccharide " \
              "residue and their\n" \
              "# Rosetta equivalents, including default Rosetta HETNAM " \
              "descriptions, which\n" \
              "# may be overridden by conflicting information from LINK " \
              "records.\n\n" \
              "# To allow these codes, use the options flag\n" \
              "# -alternate_3_letter_codes glycam.codes\n\n" \
              "# See DeMarco & Woods (2007) Glycobiology, 18(6), 426-40 for " \
              "more information.\n\n" \
              "# CODES ARE CASE-SENSITIVE!\n\n" \
              "# Code  Rosetta Code  Default HETNAM\n" \

# add each line to this list, which will then be dumped into a file
file_contents = []

# add the file header
file_contents.append(file_header)

# loop over alpha and beta codes
for anomer_type in alpha_and_beta_codes:

    # loop over all capital one letter sugar codes
    for one_letter_sugar_code in one_letter_sugar_code_list:

        # loop over the keys of the dictionaries to produce each combination of
        # sugar
        for one_letter_linkage_code in one_letter_linkage_code_list:

            ## put together the three letter code for the sugar
            # "linkage code" + "one letter sugar code" + "anomer type"
            three_letter_code = one_letter_linkage_code + \
                                one_letter_sugar_code + anomer_type

            # get the appropriate Rosetta code
            Rosetta_code = one_letter_sugar_code_dict[one_letter_sugar_code]

            # determine if this sugar has a terminal code
            if one_letter_linkage_code in terminal_values_dictionary.keys():
                terminal_code = terminal_values_dictionary[
                                               one_letter_linkage_code.upper()]
            elif one_letter_linkage_code == '0':
                terminal_code = "# terminal"
            else:
                terminal_code = ''


            ## put together the Default HETNAM
            # get the linkage direction
            default_hetnam = one_letter_linkage_code_dict[
                                                       one_letter_linkage_code]

            # determine if sugar is alpha or beta
            if anomer_type.upper() == 'A':
                default_hetnam += "-alpha"
            else:
                default_hetnam += "-beta"

            # determine if sugar is D or L
            if one_letter_sugar_code.isupper():
                default_hetnam += "-D"
            else:
                default_hetnam += "-L"

            # add the three (or more) letter Rosetta code
            default_hetnam += "-%s" %Rosetta_code

            # finally, determine if the sugar is in pyranose or furanose form
            if anomer_type.isupper():
                default_hetnam += 'p'
            else:
                default_hetnam += 'f'


            ## determine how many spaces you need between the Rosetta code and
            ## the Default HETNAM
            num_of_spaces = 14 - len(Rosetta_code)
            spaces = ' ' * num_of_spaces


            ## put together the row to be writen into the file
            row = "  %s   %s%s%s  %s\n" %(three_letter_code, Rosetta_code,
                                          spaces, default_hetnam, terminal_code)
            file_contents.append(row)



# write data to glycam.codes file
with open(filename, 'wb') as fh:
    fh.writelines(file_contents)
