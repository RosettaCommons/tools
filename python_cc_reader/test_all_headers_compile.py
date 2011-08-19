# :noTabs=true:

from test_compile import test_compile
from code_utilities import scan_compilable_files, regex_subset

from optparse import OptionParser, IndentedHelpFormatter

import re, sys


def main(args):
    '''Script for testing that each include compile on its own.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS]")
    parser.set_description(main.__doc__)

    parser.add_option("--verbose", action="store_true", dest='verbose', default=False,
      help="Run script in verbose mode. (off by default)"
    )
    (options, args) = parser.parse_args(args)


    includes = scan_compilable_files()
    re_hh_header  = re.compile("\S*\.hh$")
    re_hpp_header = re.compile( "\S*\.hpp$")

    all_files = includes.keys()

    hh_headers = regex_subset( all_files, re_hh_header )
    hpp_headers = regex_subset( all_files, re_hpp_header )

    hh_headers.extend( hpp_headers )
    headers = hh_headers

    headers.sort()

    all_compile = True
    for header in headers :
        if not test_compile( header, verbose=options.verbose ) :
            # Re-run one more time with verbose option to get exact command line was used to compile
            test_compile( header, verbose=True )
            print header, "fails to compile on its own with error message:", file("out.log").read() + '\n'
            all_compile = False
        else :
            pass
            #print header, "compiles on its own"

    if all_compile :
       sys.exit( 0 )
    else :
       sys.exit( 1 )


if __name__ == "__main__":
    main(sys.argv)
