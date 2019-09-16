from . import blargs

def is_int( val ) :
    try :
        int( val )
        return True
    except ValueError :
        return False

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "input_file" ).required()
        p.str( "output_file" ).required()
    lines = open( input_file ).readlines()
    newlines = []
    count = 0
    for line in lines :
        if len( line ) > 3 :
            if line[ 0:3 ] == "---" :
                count += 1
                cols = line.split()
                if is_int( cols[1] ) :
                    newlines.append( "--- %d %s\n" % (count, " ".join( cols[2:] )))
                else :
                    newlines.append( "--- %d %s\n" % (count, " ".join( cols[1:] )))
                continue
        newlines.append( line )
    open( output_file, "w" ).writelines( newlines )
