#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

"""File:  carbohydrate_data.py

Brief: This code contains the necessary data for GLYCAM code generation.

Details:  See glycam.org/docs/forcefield/glycam-naming-2/ for more information.

To properly add to this program, be sure to add new codes to both the lists and
dictionaries for the linkage codes and sugar codes.

Author:  Morgan Nance
Labonte <JWLabonte@jhu.edu> (PEP8/style tweaks, GLYCAM code additions)

"""

# linkage code dictionary
# { "one letter linkage code" : "Default HETNAM linkage" }
one_letter_linkage_code_list = ['0', '1', '2', '3', '4', '5', '6', 'Z', 'Y',
                                'O', 'X', 'W', 'N', 'V', 'U', 'M', 'T', 'L',
                                'S', 'R', 'K', 'Q', 'J', 'P', 'I']

one_letter_linkage_code_dict = {'0' : "->4)", '1' : "->1)", '2' : "->2)",
                                '3' : "->3)", '4' : "->4)", '5' : "->5)",
                                '6' : "->6)", 'Z' : "->2)", 'Y' : "->2)",
                                'O' : "->2)", 'X' : "->2)", 'W' : "->3)",
                                'N' : "->3)", 'V' : "->3)", 'U' : "->4)",
                                'M' : "->5)", 'T' : "->2)", 'L' : "->2)",
                                'S' : "->2)", 'R' : "->2)", 'K' : "->2)",
                                'Q' : "->3)", 'J' : "->3)", 'P' : "->2)",
                                'I' : "->2)"}

# one letter code dictionary
# iterate over keys to get the one letter sugar codes
# { "one letter GLYCAM code" : "Rosetta code" }
one_letter_sugar_code_list = ['A', 'D', 'R', 'X', 'N', 'E', 'L', 'G', 'K',
                              'I', 'M', 'T', 'C', 'P', 'B', 'J', 'F', 'Q',
                              'H', 'O', 'a', 'd', 'r', 'x', 'n', 'e', 'l',
                              'g', 'k', 'i', 'm', 't', 'c', 'p', 'b', 'j',
                              'f', 'q', 'h', 'o']

one_letter_sugar_code_dict = {'A' : "Ara", 'D' : "Lyx", 'R' : "Rib",
                              'X' : "Xyl", 'N' : "All", 'E' : "Alt",
                              'L' : "Gal", 'G' : "Glc", 'K' : "Gul",
                              'I' : "Ido", 'M' : "Man", 'T' : "Tal",
                              'C' : "Fru", 'P' : "Psi", 'B' : "Sor",
                              'J' : "Tag", 'F' : "Fuc", 'Q' : "Qui",
                              'H' : "Rha", 'O' : "Gal", 'Z' : "Glc",
                              'U' : "Ido", 'V' : "Gal", 'Y' : "Glc",
                              'W' : "Man", 'S' : "Neu", 'a' : "Ara",
                              'd' : "Lyx", 'r' : "Rib", 'x' : "Xyl",
                              'n' : "All", 'e' : "Alt", 'l' : "Gal",
                              'g' : "Glc", 'k' : "Gul", 'i' : "Ido",
                              'm' : "Man", 't' : "Tal", 'c' : "Fru",
                              'p' : "Psi", 'b' : "Sor", 'j' : "Tag",
                              'f' : "Fuc", 'q' : "Qui", 'h' : "Rha",
                              'o' : "Gal", 'z' : "Glc", 'u' : "Ido",
                              'v' : "Gal", 'y' : "Glc", 'w' : "Man",
                              's' : "Neu"}

# alpha and beta codes
# making a list only for iteration purposes
alpha_and_beta_codes = ['A', 'B', 'a', 'b']


# terminal values dictionary
# { "one letter linkage code" : "terminal values" }
terminal_values_dictionary = {'Z' : "# 2,3", 'Y' : "# 2,4", 'O' : "# 2,5",
                              'X' : "# 2,6", 'W' : "# 3,4", 'N' : "# 3,5",
                              'V' : "# 3,6", 'U' : "# 4,6", 'M' : "# 5,6",
                              'T' : "# 2,3,4", 'L' : "# 2,3,5",
                              'S' : "# 2,3,6", 'R' : "# 2,4,6",
                              'K' : "# 2,5,6", 'Q' : "# 3,4,6",
                              'J' : "# 3,5,6", 'P' : "# 2,3,4,6",
                              'I' : "# 2,4,5,6"}
