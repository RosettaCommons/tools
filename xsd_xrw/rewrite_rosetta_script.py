import blargs
import sys
#import copy

debug = False


class XMLToken :
    def __init__( self ) :
        self.line_start = -1
        self.position_start = -1
        self.line_end = -1
        self.position_end = -1
        self.contents = ""
        self.whitespace = False
        self.deleted = False
        self.index = 0
        self.in_tag = False
    @staticmethod
    def tok_for_single_pos( line_number, position_on_line, in_tag ) :
        tok = XMLToken()
        tok.line_start = line_number
        tok.line_end = line_number
        tok.position_start = position_on_line
        tok.position_end = position_on_line
        tok.in_tag = in_tag
        return tok
    def uninitialized( self ) :
        return self.line_start == -1

class Tag :
    def __init__( self ) :
        self.name = ""
        self.attributes = []
        self.tokens = []
        self.closed = False
        self.terminator = False # e.g. </SCOREFXNS>
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

def find_attribute_in_tag( tag, attribute_name ) :
    for attr in tag.attributes :
        if attr[0].contents == attribute_name :
            return attr
    return None

def name_token_of_tag( tag ) :
    for j,tok in enumerate( tag.tokens ) :
        if j == 0 : continue # skip the "<"
        for char in tok.contents :
            if char != "" and char != "\t" and char != "\n" :
                return tok
    return None

class Element :
    def __init__( self ) :
        self.sub_elements = []
        self.name = ""
        self.tags = [] # perhaps only a single tag
    def reconstitute_token_list( self, old_tok_list, new_tok_list, tok_pos ) :
        #print self.name, tok_pos, old_tok_list[ tok_pos ].contents
        assert( len( self.tags ) == 1 or len( self.tags ) == 2 )

        # append the old tokens to the new list until we reach the first token
        # that belongs to the tag; then insert all of that tag's tokens
        for which_tag in xrange( len( self.tags ) ) :
            tag_first_tok = self.tags[ which_tag ].tokens[0]
            for i in xrange( tok_pos, len( old_tok_list ) ) :
                if old_tok_list[ i ] is tag_first_tok :
                    if debug : print "Found token: ", tag_first_tok.contents, "for", self.tags[which_tag].name
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
                                if debug : print "returning", tok_pos
                                return tok_pos, new_tok_list
                            else :
                                break
                    # only reachable if we're not on the last tag
                    assert( which_tag + 1 < len( self.tags ) )
                    break
                else :
                    if debug : print self.name, "skipping", tok_pos,"\"" + old_tok_list[i].contents + "\" looking for \"" + tag_first_tok.contents + "\""
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
            #print i, j, symb, "in tag", in_tag, "in block comment", in_block_comment, "in double quote", in_double_quote, "in whitespace", in_whitespace
            if in_block_comment :
                if symb == ">" and j >= 2 and line[j-2:j+1] == "-->" :
                    assert( not curr_token.uninitialized() )
                    curr_token.line_end = i
                    curr_token.position_end = j
                    curr_token.in_tag = in_tag
                    tokens.append( curr_token )
                    curr_token = XMLToken()
                    in_block_comment = False
            elif symb == "<" :
                in_whitespace = False
                if j+3 < len(line) and line[j:j+4] == "<!--" :
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1]) - 1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_block_comment = True
                    curr_token.line_start = i
                    curr_token.position_start = j
                else :
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1]) - 1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_tag = True
                    curr_token.in_tag = in_tag
                    tokens.append( XMLToken.tok_for_single_pos( i, j, in_tag ))
            elif in_tag :
                if symb == ">" :
                    in_whitespace = False
                    end_of_element = j > 0 and line[ j-1:j+1 ] == "/>"
                    if end_of_element :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 1 else i-1
                            curr_token.position_end = j-2 if j != 1 else len(lines[i-1]) - 1
                            if curr_token.position_end >= curr_token.position_start :
                                curr_token.in_tag = in_tag
                                tokens.append( curr_token )
                            curr_token = XMLToken()
                        curr_token.line_start = i
                        curr_token.position_start = j-1
                        curr_token.line_end = i
                        curr_token.position_end = j
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    else :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 0 else i-1
                            curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                            curr_token.in_tag = in_tag
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

                            curr_token.in_tag = in_tag
                            tokens.append( curr_token )
                            curr_token = XMLToken()
                        else :
                            #if len(tokens) > 2 :
                            #    print "closing; last two? '" + contents_for_token( lines, tokens[-2] ) + "' and '" + contents_for_token( lines, tokens[-1] ) + "'"
                            curr_token.in_tag = in_tag
                            tokens.append( XMLToken.tok_for_single_pos( i, j, in_tag ))
                    in_tag = False

                elif in_double_quote or in_single_quote:
                    if ( symb == '"' and in_double_quote ) or ( symb == "'" and in_single_quote ) :
                        curr_token.line_end = i
                        curr_token.position_end = j
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                        in_double_quote = False
                        in_single_quote = False
                        #print "finished reading quote", i, j
                elif symb == '"' or symb == "'" :
                    in_whitespace = False
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        curr_token.in_tag = in_tag
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
                            curr_token.in_tag = in_tag
                            tokens.append( curr_token )
                            curr_token = XMLToken()
                        #print "Starting new non-whitespace token", i, j, symb
                        curr_token.line_start = i
                        curr_token.position_start = j
                        curr_token.whitespace = True
                        in_whitespace = True
                elif symb == "=" :
                    in_whitespace = False
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    curr_token.in_tag = in_tag
                    tokens.append( XMLToken.tok_for_single_pos( i, j, in_tag ))
                else :
                    if in_whitespace :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_whitespace = False
                    if curr_token.uninitialized() :
                        #print "Starting new non-whitespace token", i, j, symb
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
        curr_token.in_tag = in_tag
        tokens.append( curr_token )

    for i,tok in enumerate( tokens ) :
        tok.contents = contents_for_token( lines, tok )
        tok.index = i # save this for the purpose of reconstructing the
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
                    # March backwards and look to see if
                    # the first non-whitespace token behind this "/" is a "<"
                    j = i-1
                    while j >= 0 :
                        if not tokens[j].whitespace :
                            if tokens[j].contents == "<" :
                                curr_tag.terminator = True
                            break
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
                if debug : print i, tok.contents, ", ".join( [ x.contents for x in curr_tag.tokens ] )
                tags.append( curr_tag )
                curr_tag = Tag()
                in_tag = False
        elif tok.contents == "<" :
            in_tag = True
            curr_tag.tokens.append( tok )

    elements = []
    root = None
    for i,tag in enumerate( tags ) :
        if debug : print i, tag.name
        if len( elements ) > 0 :
            if elements[-1].name != tag.name :
                if tag.terminator :
                    print "Input script is not valid XML"
                    print "terminating tag \"" + name_token_of_tag( tag ).contents + "\" does not match the wrapping tag \"" + elements[-1].name + "\""
                    if len(elements) >= 2 and elements[-2].name == tag.name :
                        print "Perhaps you need to close the", elements[-1].name + " tag?"
                    sys.exit( 1 )
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

def replace_element_name_w_attribute( element, new_name, new_attribute ) :
    element.name = new_name
    for i,tag in enumerate( element.tags ) :
        old_name = tag.name
        tag.name = new_name
        # find the name token and remove it
        done = False
        for j,tok in enumerate( tag.tokens ) :
            if j == 0 : continue
            for char in tok.contents : # skip the "<"
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
            spellings = [ " ", new_attribute, "=", '"' + old_name + '"' ]
            for j,tok in enumerate( attr_tokens ) :
                tok.line_start = tok.line_end = name_token.line_start
                tok.position_start = tok.position_end = name_token.position_start # garbage
                tok.contents = spellings[j]
                tag.add_token( tok, 2+j )

def recursively_rename_subelements( root, element_name, new_subelement_name, new_attribute_name = "name" ) :
    for element in root.sub_elements :
        if element.name == element_name :
            for sfxn_element in element.sub_elements :
                replace_element_name_w_attribute( sfxn_element, new_subelement_name, new_attribute_name )
        recursively_rename_subelements( element, element_name, new_subelement_name, new_attribute_name )

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
                        replace_element_name_w_attribute( fragreader_element, "FragReader", "name" )
                else :
                    replace_element_name_w_attribute( fragset_element, "FragSet", "name" )
        rename_fragments_from_frag_reader( element )

def move_OUTPUT_as_last_child_of_ROSETTASCRIPTS( root, tokens ) :
    if root.name == "ROSETTASCRIPTS" :
        output_index = -1
        for i,elem in enumerate( root.sub_elements ) :
            if elem.name == "OUTPUT" :
                output_index = i
                break
        if output_index != -1 and output_index + 1 != len( root.sub_elements ) :
            out_elem = root.sub_elements[ output_index ]
            last_elem = root.sub_elements[ -1 ]
            prev_last = root.sub_elements[output_index-1].tags[-1].tokens[-1].index if output_index != 0 else root.tags[0].tokens[-1].index
            out_first = out_elem.tags[0].tokens[0].index
            out_last  = out_elem.tags[-1].tokens[-1].index
            next_first = root.sub_elements[output_index+1].tags[0].tokens[0].index
            next_last = root.sub_elements[output_index+1].tags[-1].tokens[-1].index
            #last_first = last_elem.tags[0].tokens[0].index
            #last_last  = last_elem.tags[-1].tokens[-1].index
            assert( len(root.tags) == 2 )
            root_first = root.tags[1].tokens[0].index
            print "OUTPUT REORDER:", prev_last+1, out_last+1, root_first
            print 0, out_first            
            print next_first, next_last+1 
            print out_last+1, next_first  
            print out_first, out_last+1   
            print next_last+1, "END"        
            #tokens = tokens[:prev_last+1] + tokens[(out_last+1):(root_first)] + tokens[prev_last+1:(out_last+1)] + tokens[(root_first):]
            tokens = tokens[ : out_first            ] + \
                     tokens[next_first: next_last+1 ] + \
                     tokens[out_last+1: next_first  ] + \
                     tokens[out_first: out_last+1   ] + \
                     tokens[next_last+1:            ]

            print "old subelement order:"
            for elem in root.sub_elements :
                print elem.name
            root.sub_elements.remove( out_elem )
            root.sub_elements.append( out_elem )

            print "new subelement order:"
            for elem in root.sub_elements :
                print elem.name
            
            renumber_tokens( tokens )

    for elem in root.sub_elements :
        tokens = move_OUTPUT_as_last_child_of_ROSETTASCRIPTS( elem, tokens )

    return tokens

def move_res_filter_as_first_child_of_OperateOnCertainResidues( root, tokens ) :
    # the first tag beneath the OperateOnCertainResidues (OperateOnResidueSubset) task operation needs
    # to be the ResFilter (ResidueSelector) and the second tag has to be the ResLvlTaskOperation
    # The valid ResLvlTaskOperations as of Nov 18 2016 are:

    # AddBehaviorRLT
    # DisallowIfNonnativeRLT
    # ExtraChiCutoffRLT
    # ExtraRotamersGenericRLT
    # IncludeCurrentRLT
    # PreserveCBetaRLT
    # PreventRepackingRLT
    # RestrictAbsentCanonicalAASRLT
    # RestrictToRepackingRLT

    if root.name == "OperateOnCertainResidues" or root.name == "OperateOnResidueSubset" :
        assert( len( root.sub_elements) == 2 or len( root.sub_elements) == 1 )
        res_lvl_task_operations = set( [ "AddBehaviorRLT",
                                         "DisallowIfNonnativeRLT",
                                         "ExtraChiCutoffRLT",
                                         "ExtraRotamersGenericRLT",
                                         "IncludeCurrentRLT",
                                         "PreserveCBetaRLT",
                                         "PreventRepackingRLT",
                                         "RestrictAbsentCanonicalAASRLT",
                                         "RestrictToRepackingRLT" ] )
        if len( root.sub_elements ) == 2 and root.sub_elements[0].name in res_lvl_task_operations :
            # print 0, root.sub_elements[0].tags[0].tokens[0].index 
            # print root.sub_elements[1].tags[0].tokens[0].index, ( root.sub_elements[1].tags[-1].tokens[-1].index + 1 )
            # print (root.sub_elements[0].tags[-1].tokens[-1].index+1), root.sub_elements[1].tags[0].tokens[0].index
            # print root.sub_elements[0].tags[0].tokens[0].index, ( root.sub_elements[0].tags[-1].tokens[-1].index + 1 )
            # print (root.sub_elements[1].tags[-1].tokens[-1].index+1), "end"
            # print
            tokens = tokens[ : root.sub_elements[0].tags[0].tokens[0].index ] + \
                tokens[ root.sub_elements[1].tags[0].tokens[0].index:( root.sub_elements[1].tags[-1].tokens[-1].index + 1 ) ] + \
                tokens[ (root.sub_elements[0].tags[-1].tokens[-1].index+1):root.sub_elements[1].tags[0].tokens[0].index ] + \
                tokens[ root.sub_elements[0].tags[0].tokens[0].index:( root.sub_elements[0].tags[-1].tokens[-1].index + 1 ) ] + \
                tokens[ (root.sub_elements[1].tags[-1].tokens[-1].index+1): ]
            root.sub_elements = [ root.sub_elements[1], root.sub_elements[0] ]
            renumber_tokens( tokens )
    
    for elem in root.sub_elements :
        tokens = move_res_filter_as_first_child_of_OperateOnCertainResidues( elem, tokens )
    return tokens

def move_fragments_as_first_child_of_fragset( root, tokens ) :
    # TO DO: make sure that FRAGMENTS element becomes first child of the
    # FRAGSET element?? Shit, the reconstitute_token_list function assumes that
    # tokens don't change their order!
    pass

def move_fragments_as_first_child_of_abscript( root, tokens ) :
    # TO DO
    # make sure that Fragments is the first sub-element of the Abscript Mover
    pass


def rename_monte_carlo_elements_from_monte_carlo_loader( root ) :
    # for elements beneath a "MONTECARLOS" element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "MonteCarlo" instead.
    recursively_rename_subelements( root, "MONTECARLOS", "MonteCarlo" )

def rename_interface_builders_from_interface_builder_loader( root ) :
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
                tok = name_token_of_tag( tag )
                assert( tok )
                assert( tok.contents == "BridgeChainsMover" or tok.contents == "/BridgeChainsMover" )
                tok.contents = "BridgeChains" if i == 0 else "/BridgeChains"

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
        tag.name = "ROSETTASCRIPTS"
    root.name = "ROSETTASCRIPTS"

def rename_report_to_db_children( root, tokens ):
    # All children of ReportToDB and TrajectoryReportToDB currently
    # have name "Feature". They now will have the name of the feature reporter
    # being parsed from them (which was formerly taken from the "name" attribute ).
    if root.name == "ReportToDB" :
        for elem in root.sub_elements :
            new_name = None
            for i,tag in enumerate( elem.tags ) :
                if i == 0 :
                    attr = find_attribute_in_tag( tag, "name" )
                    assert( attr )
                    new_name = attr[1].contents[1:-1] # trim the "s
                    attr[0].deleted = True
                    attr[1].deleted = True
                    tokens[attr[0].index+1].deleted = True
                    if tokens[attr[1].index+1].whitespace :
                        tokens[attr[1].index+1].deleted = True
                assert( new_name is not None )
                tag.name = new_name
                name_tok = name_token_of_tag( tag )
                assert( name_tok )
                name_tok.contents = new_name if i == 0 else ( "/" + new_name )
            elem.name = new_name
    for elem in root.sub_elements :
        rename_report_to_db_children( elem, tokens )

def give_all_stubsets_children_an_element_name( root ):
    # The children of the StubSets element, that is itself a subelement of mulitple different Movers,
    # need to be given the name "Add"
    if root.name == "StubSets" :
        for elem in root.sub_elements :
            for i,tag in enumerate( elem.tags ) :
                name_tok = name_token_of_tag( tag )
                assert( name_tok )
                name_tok.contents = "Add" if i == 0 else "/Add"
                tag.name = "Add"
            elem.name = "Add"
    for elem in root.sub_elements :
        give_all_stubsets_children_an_element_name( elem )

def give_parsed_protocol_children_an_element_name( root ):
    # The children of the PROTOCOLS and ParsedProtocol element
    # need to be given the name "Add" -- before, their names were ignored
    if root.name == "PROTOCOLS" or root.name == "ParsedProtocol" :
        for elem in root.sub_elements :
            for i,tag in enumerate( elem.tags ) :
                name_tok = name_token_of_tag( tag )
                assert( name_tok )
                name_tok.contents = "Add" if i == 0 else "/Add"
                tag.name = "Add"
            elem.name = "Add"
    for elem in root.sub_elements :
        give_parsed_protocol_children_an_element_name( elem )


def give_all_calculator_filter_children_an_element_name( root ):
    #Children of CalculatorFilter will either be called Var
    #(if they have the attribute filter or filter_name ) or 
    #Value (if they have the attribute value but not one of the other two). 
    #Subtags without any of these attributes are invalid.
    # test 720 not working properly!
    if root.name == "CalculatorFilter" :
        for elem in root.sub_elements :
            new_name = None
            for i,tag in enumerate( elem.tags ) :
                if i == 0 :
                    filter_attr = find_attribute_in_tag( tag, "filter" )
                    filter_name_attr = find_attribute_in_tag( tag, "filter_name" )
                    value_attr = find_attribute_in_tag( tag, "value" )
                    if filter_attr or filter_name_attr :
                        new_name = "Var"
                    elif value_attr :
                        new_name = "Value"
                    else :
                        # This shouldn't happen -- the input tag should have one of these three
                        assert( filter_attr or filter_name_attr or value_attr )
                assert( new_name )
                name_tok = name_token_of_tag( tag )
                name_tok.contents = new_name if i == 0 else ( "/" + new_name )
                tag.name = new_name
            elem.name = new_name
    for elem in root.sub_elements :
        give_all_calculator_filter_children_an_element_name( elem )

def give_all_combined_filter_children_an_element_name( root ):
    #All children of CombinedValue will now be named Add
    if root.name == "CombinedValue" :
        for elem in root.sub_elements :
            for i,tag in enumerate( elem.tags ) :
                name_tok = name_token_of_tag( tag )
                assert( name_tok )
                name_tok.contents = "Add" if i == 0 else "/Add"
                tag.name = "Add"
            elem.name = "Add"
    for elem in root.sub_elements :
        give_all_combined_filter_children_an_element_name( elem )

def give_all_generic_montecarlo_filters_an_element_name( root ) :
    # The children of the Filters element that is itself a child of the GenericMonteCarlo
    # element need to be given the name "AND"
    # This also applies to the classes that rely on the GenericMonteCarlo's structure:
    # GenericSimulatedAnnealer and EvolutionaryDynamics
    if root.name == "GenericMonteCarlo" or root.name == "GenericSimulatedAnnealer" or root.name == "EvolutionaryDynamics" :
        for elem in root.sub_elements :
            if elem.name != "Filters" : continue
            for subelement in elem.sub_elements :
                for i,tag in enumerate( subelement.tags ):
                    name_tok = name_token_of_tag( tag )
                    assert( name_tok )
                    name_tok.contents = "AND" if i == 0 else "/AND"
                    tag.name = "AND"
                subelement.name = "AND"
    for elem in root.sub_elements :
        give_all_generic_montecarlo_filters_an_element_name( elem )

def give_all_map_hotspot_Jumps_an_element_name( root ) :
    # The children of the Jumps element that is itself a child of the MapHotspot mover
    # need to be given the name "Jump"
    if root.name == "MapHotspot" :
        for elem in root.sub_elements :
            if elem.name != "Jumps" : continue
            for subelement in elem.sub_elements :
                for i,tag in enumerate( subelement.tags ) :
                    name_tok = name_token_of_tag( tag )
                    assert( name_tok )
                    name_tok.contents = "Jump" if i == 0 else "/Jump"
                    tag.name = "Jump"
                subelement.name = "Jump"
    for elem in root.sub_elements :
        give_all_map_hotspot_Jumps_an_element_name( elem )

def give_all__PlaceStub_or_PlaceSimultaneously__sub_subelements_the_name_Add( root ):
    # TO DO!!!
    #				<PlaceStub name=place_phe stubfile=native_phe_stub.pdb add_constraints=1 final_filter=hbond_ddg minimize_rb=1 hurry=1>
    #					 <DesignMovers>
    #						<Add mover_name=srsc/>
    #						<Add mover_name=des1 coord_cst_std=0.6/>
    #						<Add mover_name=des3 use_constraints=0/>
    #					</DesignMovers>
    #				</PlaceStub>
    #the tags inside DesignMovers (or, equivalently, NotifyMovers, StubMinimize, or StubSets) DO NOT have a name in parse_my_tag; the schema will call them Add
    #ALL subelements of PlaceStub and PlaceSimultaneously have this unnamed Add subelement
    pass

def give_all_dock_with_hotspots_HotspotFiles_an_element_name( root ) :
    # The children of the HotspotFiles element that is itself a child of the DockWithHotspotMover 
    # element ( and some others: SetupHotspotConstraintsLoopsMover, SetupHotspotConstraintsMover )
    # need to be given the name "HotspotFile"
    if root.name == "DockWithHotspotMover" or root.name == "SetupHotspotConstraintsLoop" or root.name == "SetupHotspotConstraintsMover" :
        for elem in root.sub_elements :
            if elem.name != "HotspotFiles" :
                #print "skipping", elem.name
                continue
            for subelement in elem.sub_elements :
                for i,tag in enumerate( subelement.tags ) :
                    name_tok = name_token_of_tag( tag )
                    assert( name_tok )
                    name_tok.contents = "HotspotFile" if i == 0 else "/HotspotFile"
                    tag.name = "HotspotFile"
                subelement.name = "HotspotFile"
    for elem in root.sub_elements :
        give_all_dock_with_hotspots_HotspotFiles_an_element_name( elem )

def rename_RotamerBoltzmannFilter_threshold_subelements( root ):
    # RotamerBoltzmannWeights subtags get renamed Threshold and old name becomes "restype"
    recursively_rename_subelements( root, "RotamerBoltzmannWeight", "Threshold", "restype" )

def rename_3mer_and_9mer_attributes_of_HybridizeMover( root ):
    # TO DO
    # attributes may not begin with a numeral, so the 3mers and 9mers attributes
    # of the HybridizeMover need to be changed to "three_mers" and "nine_mers" respectively
    if root.name == "Hybridize" :
        for elem in root.sub_elements :
            if elem.name != "Fragments" : continue
            for attr in elem.tags[0].attributes :
                if attr[0].contents == "3mers" :
                    attr[0].contents = "three_mers"
                if attr[0].contents == "9mers" :
                    attr[0].contents = "nine_mers"
    for elem in root.sub_elements :
        rename_3mer_and_9mer_attributes_of_HybridizeMover( elem )

def replace_raw_ampersand_w_and( tokens ) :
    # if there are raw ampersands in the "comments", then they need to be replaced with the
    # word and
    pass

def renumber_tokens( tokens ) :
    for i,tok in enumerate( tokens ) :
        tok.index = i

def turn_attributes_of_common_subtag_of_ModulatedMover_into_individual_subtags( root, tokens ) :
    # The common subtag of the ModulatedMover accepts any attribute as it currently stands
    # and that is a little redonk.
    # Take the "type" attribute out of the ModulatedMover, and that will become the new
    if root.name == "ModulatedMover" :
        common_element = None
        mover_attributes = [] # the set of tokens describing this
        mover_name = None
        for attr in root.tags[0].attributes :
            if attr[0].contents == "type" :
                mover_name = attr[1].contents[1:-1]
                for i in xrange( attr[0].index, attr[1].index+1 ) :
                    tokens[i].deleted = True
                if tokens[attr[1].index+1].whitespace :
                    tokens[attr[1].index+1].deleted = True
                break
        assert( mover_name is not None )
        for elem in root.sub_elements:
            if elem.name == "common" :
                common_element = elem
                for attr in elem.tags[0].attributes :
                    mover_attributes.append( ( attr[0].contents, attr[1].contents ) )
        for elem in root.sub_elements :
            if elem.name == "Interp" :
                seed_val = None
                key_val = None
                for attr in elem.tags[0].attributes :
                    if attr[0].contents == "start" :
                        seed_val = attr[1].contents
                    if attr[0].contents == "key" :
                        key_val = attr[1].contents
                    if attr[0].contents == "value" :
                        seed_val = attr[1].contents
                assert( seed_val is not None and key_val is not None )
                mover_attributes.append( ( key_val[1:-1], seed_val ) )
        new_mover_tag_line = [ "<", mover_name ]
        new_mover_tag_line.append( "".join( [ " " + x[0] + "=" + x[1] for x in mover_attributes ] ) )
        new_mover_tag_line.append( "/>" )
        new_mover_tag_line = [ "".join( new_mover_tag_line ) ]
        new_tokens = tokenize_lines( new_mover_tag_line )
        tags, new_mover_element = tokens_into_tags( new_tokens )

        first_root_subelement_token = root.sub_elements[0].tags[0].tokens[0].index
        tokens = tokens[:first_root_subelement_token] + new_tokens + tokens[first_root_subelement_token:]
        #for tok in tokens : print tok.contents,

        root.sub_elements.insert( 0, new_mover_element )
        # now, remove the old common element
        if common_element :
            for tag in common_element.tags :
                for tok in tag.tokens :
                    tok.deleted = True
        renumber_tokens( tokens )

    for elem in root.sub_elements :
        tokens = turn_attributes_of_common_subtag_of_ModulatedMover_into_individual_subtags( elem, tokens )
    return tokens

def turn_wild_ampersands_into_and( tokens ) :
    for tok in tokens:
        if tok.in_tag : continue
        print "comment:", tok.contents
        tok.contents = tok.contents.replace( "&&", "AND" )
        tok.contents = tok.contents.replace( "&", "and" )

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
    turn_wild_ampersands_into_and( toks )
    #print_element( 0, element_root )

    surround_attributes_w_quotes( tags )
    modifications = [ rename_score_functions, rename_fragments_from_frag_reader,
                      rename_monte_carlo_elements_from_monte_carlo_loader,
                      rename_interface_builders_from_interface_builder_loader,
                      rename_movemaps_from_movemap_loader,
                      rename_ligand_areas_from_ligand_area_loader,
                      rename_bridge_chains_mover_to_bridge_chains,
                      rename_dockdesign_to_ROSETTASCRIPTS,
                      give_all_stubsets_children_an_element_name,
                      give_all_calculator_filter_children_an_element_name,
                      give_all_combined_filter_children_an_element_name,
                      give_all_generic_montecarlo_filters_an_element_name,
                      give_all_map_hotspot_Jumps_an_element_name,
                      give_all_dock_with_hotspots_HotspotFiles_an_element_name,
                      rename_RotamerBoltzmannFilter_threshold_subelements,
                      give_parsed_protocol_children_an_element_name,
                      rename_3mer_and_9mer_attributes_of_HybridizeMover
                  ]

    for modfunc in modifications :
        modfunc( element_root )

    rename_report_to_db_children( element_root, toks )

    # big modification to the behavior of the ModulatedMover XML structure
    toks = turn_attributes_of_common_subtag_of_ModulatedMover_into_individual_subtags( element_root, toks )
    toks = move_res_filter_as_first_child_of_OperateOnCertainResidues( element_root, toks )

    #for tok in toks : print tok.contents,

    last_tok_index, new_toks =  element_root.reconstitute_token_list( toks, [], 0 )
    new_toks.extend( toks[last_tok_index:])
    print "How many tokens at the end?", last_tok_index, len( toks )
    print "\n".join( [ ( "remainder:" + x.contents ) for x in toks[ last_tok_index: ] ] )

    mostly_rewritten_version = "".join( [ (x.contents if not x.deleted else "") for x in new_toks ] ) + "\n"

    lines2 = [ x + "\n" for x in mostly_rewritten_version.split( "\n" ) ]
    toks2 = tokenize_lines( lines2 )
    tags, element_root = tokens_into_tags( toks2 )
    root_first_tok = element_root.tags[0].tokens[0]
    if ( root_first_tok.line_start != 0 or root_first_tok.position_start != 0 ) :
        lines2 = [ "<ROSETTASCRIPTS>" ] + lines2[:root_first_tok.line_start ] + lines2[root_first_tok.line_start + 1: ]

    #debug = True

    toks3 = tokenize_lines( lines2 )
    tags, element_root = tokens_into_tags( toks3 )
    toks3 = move_OUTPUT_as_last_child_of_ROSETTASCRIPTS( element_root, toks3 )
    last_tok_index, new_toks = element_root.reconstitute_token_list( toks3, [], 0 )
    new_toks.extend( toks3[last_tok_index:])

    #print "rewritten version:"
    #print 
    #open( output, "w" ).writelines( [ x + "\n" for x in lines2 ] )
    open( output, "w" ).write( "".join( [ (x.contents if not x.deleted else "") for x in new_toks ] ) )

    #for i,tag in enumerate( tags ) :
    #    print i, "".join( [ x.contents for x in tag.tokens ] )
