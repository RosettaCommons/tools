import blargs

class XMLToken :
    def __init__( self ) :
        self.line_start = -1
        self.position_start = -1
        self.line_end = -1
        self.position_end = -1
        self.contents = ""
        self.whitespace = False
        self.deleted = False
    @staticmethod
    def tok_for_single_pos( line_number, position_on_line ) :
        tok = XMLToken()
        tok.line_start = line_number
        tok.line_end = line_number
        tok.position_start = position_on_line
        tok.position_end = position_on_line
        return tok
    def uninitialized( self ) :
        return self.line_start == -1

class Tag :
    def __init__( self ) :
        self.name = ""
        self.attributes = []
        self.tokens = []
        self.closed = False
        self.new_tokens_since_parsing = False

    def remove_token( self, token ) :
        #print "remove token", token.contents, ",".join( [ x.contents for x in self.tokens ] )
        self.tokens.remove( token )
        self.new_tokens_since_parsing = True

    def add_token( self, token, pos ) :
        # you cannot replace either the first or the last token
        # which you shouldn't need to, because these are the "<" and
        # ">" symbols.
        #print "add token before: ", "".join( [ x.contents for x in self.tokens ] ), pos, token.contents
        assert( pos != 0 and pos != len(self.tokens) )
        self.tokens.insert( pos, token )
        #print "add token after: ", "".join( [ x.contents for x in self.tokens ] )
        self.new_tokens_since_parsing = True

    def replace_tokens_in_list( self, token_list ) :
        #slow, don't use this
        if not self.new_tokens_since_parsing : return token_list
        for i in xrange( len( token_list ) ) :
            # look for the first token and then the last token
            if self.tokens[0] is token_list[ i ] :
                for j in xrange( i+1, len( token_list )) :
                    if self.tokens[-1] is token_list[ j ] :
                        token_list = token_list[ :i ] + self.tokens + token_list[ j+1: ]
                        return token_list
        assert( False )

class Element :
    def __init__( self ) :
        self.sub_elements = []
        self.name = ""
        self.tags = [] # perhaps only a single tag
    def reconstitute_token_list( self, old_tok_list, new_tok_list, tok_pos ) :
        assert( len( self.tags ) == 1 or len( self.tags ) == 2 )

        # append the old tokens to the new list until we reach the first token
        # that belongs to the tag; then insert all of that tag's tokens
        for which_tag in xrange( len( self.tags ) ) :
            tag_first_tok = self.tags[ which_tag ].tokens[0]
            for i in xrange( tok_pos, len( old_tok_list ) ) :
                if old_tok_list[ i ] is tag_first_tok :
                    new_tok_list.extend( self.tags[ which_tag ].tokens )
                    tag_last_tok = self.tags[ which_tag ].tokens[-1]
                    for j in xrange( i+1, len( old_tok_list ) ) :
                        if old_tok_list[ j ] == tag_last_tok :
                            tok_pos = j+1
                            if which_tag == 0 :
                                # recurse!
                                for sub_element in self.sub_elements :
                                    tok_pos, new_tok_list = sub_element.reconstitute_token_list( old_tok_list, new_tok_list, tok_pos )
                            if which_tag + 1 == len( self.tags ) :
                                return tok_pos, new_tok_list
                            else :
                                break
                    # only reachable if we're not on the last tag
                    assert( which_tag + 1 < len( self.tags ) )
                    break
                else :
                    new_tok_list.append( old_tok_list[ i ] )

        unreachable = False
        assert( unreachable ) # we should not reach here

def tokenize_lines( lines ) :
    in_tag = False
    in_block_comment = False
    in_double_quote = False
    in_single_quote = False
    in_whitespace = False
    tokens = []

    curr_token = XMLToken()
    for i,line in enumerate(lines) :
        for j,symb in enumerate( line ) :
            #print i, j, symb, "in tag", in_tag, "in block comment", in_block_comment, "in quote", in_quote, "in whitespace", in_whitespace
            if in_block_comment :
                if symb == ">" and j > 2 and line[j-2:j+1] == "-->" :
                    assert( not curr_token.uninitialized() )
                    curr_token.line_end = i
                    curr_token.position_end = j
                    tokens.append( curr_token )
                    curr_token = XMLToken()
                    in_block_comment = False
            elif symb == "<" :
                in_whitespace = False
                if j+3 < len(line) and line[j:j+4] == "<!--" :
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1]) - 1
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_block_comment = True
                    curr_token.line_start = i
                    curr_token.position_start = j
                else :
                    in_tag = True
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1]) - 1
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    tokens.append( XMLToken.tok_for_single_pos( i, j ))
            elif in_tag :
                if symb == ">" :
                    in_whitespace = False
                    in_tag = False
                    end_of_element = j > 0 and line[ j-1:j+1 ] == "/>"
                    if end_of_element :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 1 else i-1
                            curr_token.position_end = j-2 if j != 1 else len(lines[i-1]) - 1
                            if curr_token.position_end >= curr_token.position_start :
                                tokens.append( curr_token )
                            curr_token = XMLToken()
                        curr_token.line_start = i
                        curr_token.position_start = j-1
                        curr_token.line_end = i
                        curr_token.position_end = j
                        tokens.append( curr_token )
                        curr_token = XMLToken()
    
                    else :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 0 else i-1
                            curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                            tokens.append( curr_token )
                            curr_token = XMLToken()

                        # look at the previous two tokens; if the last token was whitespace
                        # and the one before was "/", then make this token a "/>"
                        if len(tokens) > 2 and tokens[-1].whitespace and contents_for_token( lines, tokens[-2] ) == "/" :
                            #print "FOUND IT"
                            curr_token.line_start = i; curr_token.position_start = tokens[-2].position_start
                            curr_token.line_end = i; curr_token.position_end = j

                            tokens.pop(); # destroy the last two tokens
                            tokens.pop();

                            tokens.append( curr_token )
                            curr_token = XMLToken()
                        else :
                            #if len(tokens) > 2 :
                            #    print "closing; last two? '" + contents_for_token( lines, tokens[-2] ) + "' and '" + contents_for_token( lines, tokens[-1] ) + "'"
                            tokens.append( XMLToken.tok_for_single_pos( i, j ))
                elif in_double_quote or in_single_quote:
                    if ( symb == '"' and in_double_quote ) or ( symb == "'" and in_single_quote ) :
                        curr_token.line_end = i
                        curr_token.position_end = j
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                        in_double_quote = False
                        in_single_quote = False
                elif symb == '"' or symb == "'" :
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    if symb == '"' : in_double_quote = True
                    else           : in_single_quote = True
                    curr_token.line_start = i
                    curr_token.position_start = j
                elif symb == " " or symb == "\t" or symb == "\n" :
                    if in_whitespace :
                        pass
                    else :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 0 else i-1
                            curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                            tokens.append( curr_token )
                            curr_token = XMLToken()
                        curr_token.line_start = i
                        curr_token.position_start = j
                        curr_token.whitespace = True
                        in_whitespace = True
                elif symb == "=" :
                    in_whitespace = False
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    tokens.append( XMLToken.tok_for_single_pos( i, j ))
                else :
                    if in_whitespace :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_whitespace = False
                    if curr_token.uninitialized() :
                        curr_token.line_start = i
                        curr_token.position_start = j
            else :
                # outside of a tag; just nab all the text and put it into a single token
                # until we get to something interesting
                if curr_token.uninitialized() :
                    curr_token.line_start = i
                    curr_token.position_start = j
    #end for
    if not curr_token.uninitialized() :
        curr_token.line_end = len( lines ) - 1
        curr_token.position_end = len( lines[-1] ) - 1
        tokens.append( curr_token )

    for i,tok in enumerate( tokens ) :
        tok.contents = contents_for_token( lines, tok )
        #print i, tok.contents

    for tok in tokens :
        # look for tokens that are "/" + whitespace + ">"
        # and change them to just be "/>"
        # because you cannot have whitespace at the close of a tag
        if len(tok.contents) > 2 and tok.contents[0] == "/" and tok.contents[-1] == ">" :
            all_whitespace = True
            for i in xrange(1,len(tok.contents)-1) :
                if tok.contents[i] == " " or tok.contents[i] == "\n" or tok.contents[i] == "\t" : continue
                all_whitespace = False
                break
            if all_whitespace :
                tok.contents = "/>"

    return tokens

def contents_for_token( lines, tok ) :
    tok_lines = []
    for lineno in xrange(tok.line_start, tok.line_end + 1 ) :
        if lineno == tok.line_start :
            if lineno == tok.line_end :
                tok_lines.append( lines[ lineno ][ tok.position_start:(tok.position_end+1) ] )
            else :
                tok_lines.append( lines[ lineno ][ tok.position_start: ] )
        elif lineno == tok.line_end :
            tok_lines.append( lines[ lineno ][ :( tok.position_end + 1 ) ] )
        else :
            tok_lines.append( lines[ lineno ] )
    return "".join( tok_lines )


def tokens_into_tags( tokens ) :
    # interpets tokens as XML tags.
    # marks whitespace tokens between attribute names, the equals sign, and their values
    # as deleted.
    tags = []
    in_tag = False
    in_attribute = False

    curr_tag = Tag()
    for i,tok in enumerate( tokens ) :

        if in_tag :
            curr_tag.tokens.append( tok )
            if curr_tag.name == "" :
                if tok.contents[0] == "/" :
                    curr_tag.name = tok.contents[1:]
                    curr_tag.closed = True
                elif tok.whitespace :
                    # you cannot have whitespace between the opening "<" and the tag name
                    tok.deleted = True
                else :
                    curr_tag.name = tok.contents
            elif in_attribute :
                #assert( len( tok.contents ) == 1 )
                if tok.whitespace :
                    tok.deleted = True
                else :
                    curr_tag.attributes[-1].append( tok )
                    in_attribute = False
            elif tok.contents == "=" :
                in_attribute = True
                # now we need to go backwards to find the first non-whitespace token before
                # the equals sign; mark all of the whitespace tokens in between as deleted.
                j = i-1
                while ( j >= 0 ) :
                    if not tokens[j].whitespace :
                        curr_tag.attributes.append( [ tokens[j] ] )
                        break
                    else :
                        tokens[j].deleted = True
                    j -= 1
            elif tok.contents[-1] == ">" :

                if tok.contents == "/>" :
                    curr_tag.closed = True
                #print i, tok.contents, ", ".join( [ x.contents for x in curr_tag.tokens ] )
                tags.append( curr_tag )
                curr_tag = Tag()
                in_tag = False
        elif tok.contents == "<" :
            in_tag = True
            curr_tag.tokens.append( tok )

    elements = []
    root = None
    for i,tag in enumerate( tags ) :
        #print i, tag.name
        if len( elements ) > 0 :
            if elements[-1].name != tag.name :
                elements.append( Element() )
                elements[-1].name = tag.name
                elements[-2].sub_elements.append( elements[-1] )
            elements[-1].tags.append( tag )
            #print "appending tag", tag.name, "to", elements[-1].name, len(elements[-1].tags)
        else :
            elements.append( Element() )
            elements[-1].tags.append( tag )
            elements[-1].name = tag.name
            root = elements[-1]
            #print "appending root tag"
        if tag.closed :
            elements.pop()
    return tags, root

# ok -- take a set of tags, and edit them so that the rhs of each attribute statement
# is surrounded by quotes
def surround_attributes_w_quotes( tags ) :
    for tag in tags :
        for attr in tag.attributes :
            rhs = attr[1]
            if rhs.contents[0] == '"' and rhs.contents[-1] == '"' : continue
            if rhs.contents[0] == "'" and rhs.contents[-1] == "'" : continue
            assert (( rhs.contents[0] != '"' and rhs.contents[0] != "'" ) or
                    rhs.contents[0] == rhs.contents[-1] )
            rhs.contents = '"' + rhs.contents + '"'

def replace_element_name_w_name_attribute( element, new_name ) :
    element.name = new_name
    for i,tag in enumerate( element.tags ) :
        old_name = tag.name
        tag.name = new_name
        # find the name token and remove it
        done = False
        for j,tok in enumerate( tag.tokens ) :
            if j == 0 : continue
            for char in tok.contents :
                if char != " " and char != "\t" and char != "\n" :
                    tag.remove_token( tok )
                    done = True
                    break
            if done : break
        #and now add a token for the name
        name_token = XMLToken()
        name_token.line_start = tag.tokens[0].line_start
        name_token.line_end = tag.tokens[0].line_end
        name_token.position_start = tag.tokens[0].position_start + 1
        name_token.position_end = tag.tokens[0].position_start + 1 # garbage
        name_token.contents = new_name if i == 0 else ("/" + new_name )
        tag.add_token( name_token,  1 )
        if i == 0 :
            attr_tokens = [ XMLToken(), XMLToken(), XMLToken(), XMLToken() ]
            spellings = [ " ", "name", "=", '"' + old_name + '"' ]
            for j,tok in enumerate( attr_tokens ) :
                tok.line_start = tok.line_end = name_token.line_start
                tok.position_start = tok.position_end = name_token.position_start # garbage
                tok.contents = spellings[j]
                tag.add_token( tok, 2+j )

def recursively_rename_subelements( root, element_name, new_subelement_name ) :
    for element in root.sub_elements :
        if element.name == element_name :
            for sfxn_element in element.sub_elements :
                replace_element_name_w_name_attribute( sfxn_element, new_subelement_name )
        rename_score_functions( element )

def rename_score_functions( root ) :
    recursively_rename_subelements( root, "SCOREFXNS", "ScoreFunction" )

def rename_fragments_from_frag_reader( root ) :
    # for elements beneath a "FRAGSETS" element:
    # 1. elements beneath FRAGMENTS element gets renamed as
    # "FragReader" and the old name becomes the "name" attribute.
    # 2. everything else beneath the "FRAGSETS" element
    # gets renamed "FragSet" and the old name becomes the name attribute
    # 3. The FRAGMENTS element must appear as the first child of the
    # FragSetLoader. It's also required??
    for element in root.sub_elements :
        if element.name == "FRAGSETS" :
            for fragset_element in element.sub_elements :
                if fragset_element.name == "FRAGMENTS" :
                    for fragreader_element in fragset_element.sub_elements  :
                        replace_element_name_w_name_attribute( fragreader_element, "FragReader" )
                else :
                    replace_element_name_w_name_attribute( fragset_element, "FragSet" )
            # TO DO: make sure that FRAGMENTS element becomes first child of the
            # FRAGSET element?? Shit, the reconstitute_token_list function assumes that
            # tokens don't change their order!
        rename_fragments_from_frag_reader( element )
                

def rename_monte_carlo_elements_from_monte_carlo_loader( root ) :
    # TO DO!!!
    # for elements beneath a "MONTECARLOS" element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "MonteCarlo" instead.
    recursively_rename_subelements( root, "MONTECARLOS", "MonteCarlo" )

def rename_interface_builders_from_interface_builder_loader( root ) :
    # TO DO!!
    # for elements beneath an INTERFACE_BUILDERS element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "InterfaceBuilder" instead.
    recursively_rename_subelements( root, "INTERFACE_BUILDERS", "InterfaceBuilder" )

def rename_movemaps_from_movemap_loader( root ) :
    # for elements beneath a MOVEMAP_BUILDERS element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "MoveMapBuilder" instead.
    recursively_rename_subelements( root, "MOVEMAP_BUILDERS", "MoveMapBuilder" )

def rename_ligand_areas_from_ligand_area_loader( root ) :
    # for elements beneath a LIGAND_AREAS element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "LigandArea" instead.
    recursively_rename_subelements( root, "LIGAND_AREAS", "LigandArea" )


def rename_bridge_chains_mover_to_bridge_chains( root ) :
    # "BridgeChainsMover" has been eliminated, now only "BridgeChains"
    # is acceptible.
    for element in root.sub_elements :
        if element.name == "BridgeChainsMover" :
            element.name = "BridgeChains"
            for i, tag in enumerate( element.tags ) :
                tag.name = "BridgeChains"
                done = False
                for j,tok in enumerate( tag.tokens ) :
                    if j == 0 : continue
                    for char in tok.contents :
                        if char != " " and char != "\t" and char != "\n" :
                            assert( tok.contents == "BridgeChainsMover" or tok.contents == "/BridgeChainsMover" )
                            tok.contents = "BridgeChains" if i == 0 else "/BridgeChains"
                            done = True
                            break
                    if done : break
        rename_bridge_chains_mover_to_bridge_chains( element )

def rename_dockdesign_to_ROSETTASCRIPTS( root ) :
    # dock_design is no longer an acceptible starting tag for a rosetta script
    for i,tag in enumerate( root.tags ) :
        #print "tag:", i, tag.name
        for j,tok in enumerate( tag.tokens ) :
            #print "tok:", j, tok.contents
            if j == 0 : continue
            if tok.deleted: continue
            assert( tok.contents == tag.name or tok.contents == "/" + tag.name )
            tok.contents = "ROSETTASCRIPTS" if i==0 else "/ROSETTASCRIPTS"
            break

def rename_RotamerBoltzmannFilter_threshold_subelements( root ):
	# RotamerBoltzmannWeights subtags get renamed Threshold and old name becomes "restype"
	pass

def print_element( depth, element ) :
    print "-" * depth, element.name
    for child in element.sub_elements :
        print_element( depth + 1, child )

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "input" )
        p.str( "output" )

    lines = open( input ).readlines()
    #newlines = add_attribute_quotes( lines )
    #open( output, "w" ).writelines( newlines )
    toks = tokenize_lines( lines )
    #for i,tok in enumerate( toks ) :
    #    for j,line in enumerate(tok.contents) :
    #        if len(tok.contents) == 1 :
    #            print "tok: %4d" % i, line
    #        elif j == 0 :
    #            print "tok: %4d" % i, line
    #        else :
    #            print "         ", line
    
    tags, element_root = tokens_into_tags( toks )
    #print_element( 0, element_root )
    
    surround_attributes_w_quotes( tags )
    modifications = [ rename_score_functions, rename_fragments_from_frag_reader,
                      rename_monte_carlo_elements_from_monte_carlo_loader,
                      rename_interface_builders_from_interface_builder_loader,
                      rename_movemaps_from_movemap_loader,
                      rename_ligand_areas_from_ligand_area_loader,
                      rename_bridge_chains_mover_to_bridge_chains,
                      rename_dockdesign_to_ROSETTASCRIPTS ]

    for modfunc in modifications :
        modfunc( element_root )
    
    dummy, new_toks =  element_root.reconstitute_token_list( toks, [], 0 )

    #print "rewritten version:"
    #print "".join( [ (x.contents if not x.deleted else "") for x in new_toks ] )
    open( output, "w" ).write( "".join( [ (x.contents if not x.deleted else "") for x in new_toks ] ) + "\n" )

    #for i,tag in enumerate( tags ) :
    #    print i, "".join( [ x.contents for x in tag.tokens ] )
