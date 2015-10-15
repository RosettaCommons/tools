#!/usr/bin/python

########################################################################
# see glycam.org/docs/forcefield/glycam-naming-2/ for more information #
########################################################################


########################################################################
# to properly add to this program, be sure to add new codes to both    #
# the lists and dictionaries for the linkage codes and sugar codes     #
########################################################################


import os


# filename preparation
working_dir = os.getcwd() + '/'
filename = working_dir + "glycam.codes"
file_header = "# A list of 3-letter GLYCAM codes for monosaccharide residues \
 and their\n# Rosetta equivalents, including default Rosetta HETNAM descriptions, \
which\n# may be overridden by conflicting information from LINK records.\n\n# \
To allow these codes, use the options flag\n# -alternate_3_letter_codes \
glycam.codes\n\n# See DeMarco & Woods (2007) Glycobiology, 18(6), 426-40 \
for more information.\n\n# CODES ARE CASE-SENSITIVE!\n\n# Code  Rosetta Code  Default HETNAM\n"

# linkage code dictionary
# { "one letter linkage code" : "Default HETNAM linkage" }
one_letter_linkage_code_list = [ '0', '1', '2', '3', '4', '5', '6', 'Z', 'Y', 'O', \
'X', 'W', 'N', 'V', 'U', 'M', 'T', 'L', 'S', 'R', 'K', 'Q', 'J', 'P', 'I' ]

one_letter_linkage_code_dict = { '0' : "->4)", '1' : "->1)", '2' : "->2)", '3' : "->3)", \
'4' : "->4)", '5' : "->5)", '6' : "->6)", 'Z' : "->2)", 'Y' : "->2)", 'O' : "->2)", \
'X' : "->2)", 'W' : "->3)", 'N' : "->3)", 'V' : "->3)", 'U' : "->4)", 'M' : "->5)", \
'T' : "->2)", 'L' : "->2)", 'S' : "->2)", 'R' : "->2)", 'K' : "->2)", 'Q' : "->3)", \
'J' : "->3)", 'P' : "->2)", 'I' : "->2)" }


# one letter code dictionary
# iterate over keys to get the one letter sugar codes
# { "one letter GLYCAM code" : "Rosetta code" }
one_letter_sugar_code_list = [ 'A', 'D', 'R', 'X', 'N', 'E', 'L', 'G', 'K', \
'I', 'M', 'T', 'C', 'P', 'B', 'J', 'F', 'Q', 'H', 'a', 'd', 'r', 'x', 'n', \
'e', 'l', 'g', 'k', 'i', 'm', 't', 'c', 'p', 'b', 'j', 'f', 'q', 'h' ]

one_letter_sugar_code_dict = { 'A' : "Ara", 'D' : "Lyx", 'R' : "Rib", 'X' : "Xyl", \
'N' : "All", 'E' : "Alt", 'L' : "Gal", 'G' : "Glc", 'K' : "Gul", 'I' : "Ido", \
'M' : "Man", 'T' : "Tal", 'C' : "Fru", 'P' : "Psi", 'B' : "Sor", 'J' : "Tag", \
'F' : "Fuc", 'Q' : "Qui", 'H' : "Rha", 'a' : "Ara", 'd' : "Lyx", 'r' : "Rib", 'x' : "Xyl", \
'n' : "All", 'e' : "Alt", 'l' : "Gal", 'g' : "Glc", 'k' : "Gul", 'i' : "Ido", \
'm' : "Man", 't' : "Tal", 'c' : "Fru", 'p' : "Psi", 'b' : "Sor", 'j' : "Tag", \
'f' : "Fuc", 'q' : "Qui", 'h' : "Rha" }

# alpha and beta codes
# making a list only for iteration purposes
alpha_and_beta_codes = [ 'A', 'B', 'a', 'b' ]


# terminal values dictionary
# { "one letter linkage code" : "terminal values" }
terminal_values_dictionary = { 'Z' : "# 2,3", 'Y' : "# 2,4", 'O' : "# 2,5", \
'X' : "# 2,6", 'W' : "# 3,4", 'N' : "# 3,5", 'V' : "# 3,6", 'U' : "# 4,6", \
'M' : "# 5,6", 'T' : "# 2,3,4", 'L' : "# 2,3,5", 'S' : "# 2,3,6", 'R' : "# 2,4,6", \
'K' : "# 2,5,6", 'Q' : "# 3,4,6", 'J' : "# 3,5,6", 'P' : "# 2,3,4,6", 'I' : "# 2,4,5,6" }


with open( filename, 'wb' ) as fh:
    # write the file header
    fh.write( file_header )
    
    # loop over alpha and beta codes
    for anomer_type in alpha_and_beta_codes:
                
        # loop over all capital one letter sugar codes
        for one_letter_sugar_code in one_letter_sugar_code_list:
            
            # loop over the keys of the dictionaries to produce each combination of sugar
            for one_letter_linkage_code in one_letter_linkage_code_list:
            
                ## put together the three letter code for the sugar
                # "linkage code" + "one letter sugar code" + "anomer type"
                three_letter_code = one_letter_linkage_code + one_letter_sugar_code + anomer_type
                
                # get the appropriate Rosetta code
                Rosetta_code = one_letter_sugar_code_dict[ one_letter_sugar_code ]
                
                # determine if this sugar has a terminal code
                if one_letter_linkage_code in terminal_values_dictionary.keys():
                    terminal_code = terminal_values_dictionary[ one_letter_linkage_code.upper() ]
                elif one_letter_linkage_code == '0':
                    terminal_code = "# terminal"
                else:
                    terminal_code = ''
                
                    
                ## put together the Default HETNAM
                # get the linkage direction
                default_hetnam = one_letter_linkage_code_dict[ one_letter_linkage_code ]
                
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
                
                    
                ## determine how many spaces you need between the Rosetta code and the Default HETNAM
                num_of_spaces = 14 - len( Rosetta_code )
                spaces = ' ' * num_of_spaces
                
                
                ## put together the row to be writen into the file
                row = "  %s   %s%s%s  %s\n" %( three_letter_code, Rosetta_code, spaces, default_hetnam, terminal_code )
                fh.write( row )

            fh.write( "\n" )
