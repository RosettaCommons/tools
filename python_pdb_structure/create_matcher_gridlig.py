from pdb_structure import *
from optparse import OptionParser

# output, sent to stdout, will look something like this
#NAME: gridlig
#BASE: 20.586 1.231 67.446
#SIZE: 105   122   147
#LENGTH: 0.500   0.500   0.500

def initialize_options_parser() :
        parser = OptionParser()
        parser.add_option( "-s", "--pdb-structure", dest="struct", type="string", help="PDB file to find the bounding box for", default="" )
        parser.add_option( "-x", "--extension", dest="extension", type="float", help="distance to increase the matcher's bounding box beyond the bounding box of the input pdb", default="5" );
        return parser

def gridlig_bounding_box_lines( pose, extra ) :
    low, hi = xyz_limits_for_pdb( pose )
    low.x_ -= extra
    low.y_ -= extra
    low.z_ -= extra
    hi.x_ += extra
    hi.y_ += extra
    hi.z_ += extra
    lines = []
    lines.append( "NAME: gridlig" + "\n" )
    lines.append( "BASE: %7.3f %7.3f %7.3f" % ( low.x_, low.y_, low.z_ ) + "\n" )
    lines.append( "SIZE: %7d %7d %7d" % ( (( hi.x_ - low.x_ )*2 + 1 ), (( hi.y_ - low.y_)*2 + 1 ), ((hi.z_ - low.z_)*2+1)) + "\n" )
    lines.append( "LENGTH:  0.50    0.50    0.50" + "\n" )
    return lines

if __name__ == "__main__" :
    parser = initialize_options_parser()
    (options, args) = parser.parse_args()
    pose = pdbstructure_from_file( options.struct )
    lines = gridlig_bounding_box_lines( pose, options.extension )
    for line in lines :
        print(line, end=' ')
