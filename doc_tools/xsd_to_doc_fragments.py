#!/usr/bin/env python
'''xsd_to_doc_fragments.py - Convert a RosettaScripts XSD into fragments of Markdown.

The current RosettaScript XSD schema file can be obtained by running rosetta_scripts with the `-parser::output_schema rosettascripts.xsd` option.

You can then run this script to (re)generate the Documentation Repository's XSD fragments.

    ./xsd_to_doc_fragments.py rosettascripts.xsd  $ROSETTA/documentation/scripting_documentation/RosettaScripts/xsd

'''

import sys, os
import xml.etree.ElementTree as ET
import codecs

TYPE_ALIASES={ 'xs:string':'string', 'xs:integer':'integer', 'xs:decimal':'real', 'rosetta_bool':'bool' }
# SECTIONS are actually the top level groups (xs:group) entries in the XSD (though not the nonce ones).
SECTIONS = [ 'mover', 'filter', 'packer_palette', 'task_operation', 'residue_selector', 'simple_metric', 'res_lvl_task_op', 'res_filter', 'constraint_generator', 'features_reporter', 'pose_selector', 'scoring_grid', 'pose_property_reporter', 'denovo_perturber', 'compound_architect_pairing_group', 'jump_selector', 'loop_definer' ]
# Nonce groups:
# ['loop_definer', 'layer_design_ss_layer', 'layer_design_ss_layer_or_taskop', 'denovo_architect', 'compound_architect_pairing_group', 'denovo_perturber', 'denovo_folder', 'rdf_function', 'pose_property_reporter', 'envclaim', 'scriptcm']

# These are the high-level groups with uniform subentries that aren't encapsulated by an xs:group -- This actually works with the SECTIONS parsing machinery.
SECTIONS += [ 'movemap_factory_loader_MOVE_MAP_FACTORIES_type', 'sfxn_loader_SCOREFXNS_type', 'interface_builder_INTERFACE_BUILDERS_type', 'ligand_area_builder_LIGAND_AREAS_type', 'move_map_builder_MOVEMAP_BUILDERS_type', ]

#'database_session_loader_DATABASE_SESSIONS_type', monte_carlo_loader_MONTECARLOS_type, rosetta_scripts_parser_IMPORT_type, rosetta_scripts_parser_OUTPUT_type

# SINGLE_TYPES are high-level groups where the tag is defined in-type. We load these directly
SINGLE_TYPES = [ ('rosetta_scripts_parser_OUTPUT_type','OUTPUT'), ("rosetta_scripts_parser_IMPORT_type","IMPORT"), ("database_session_loader_DATABASE_SESSIONS_type","DATABASE_SESSIONS"), ('monte_carlo_loader_MONTECARLOS_type',"MONTECARLOS"), ]

# These are type definitions to eliminate circular dependencies in documentation
COMMON_TYPES={ # typename:(pseudoname,docstring)
'rosetta_scripts_parser_ROSETTASCRIPTS_type':('ROSETTASCRIPTS','A full [[RosettaScripts]] protocol, as a subtag'),

'mover':('Mover Tag','Any of the [[RosettaScripts Mover|Movers-RosettaScripts]] tags'),
'filter':('Filter Tag','Any of the [[RosettaScripts Filters|Filters-RosettaScripts]] tags'),
'packer_palette':('PackerPalette Tag','Any of the [[RosettaScripts PackerPalettes|PackerPalette]] tags'),
'task_operation':('TaskOperation Tag','Any of the [[RosettaScripts TaskOperation|TaskOperations-RosettaScripts]] tags'),
'residue_selector':('Residue Selector Tag','Any of the [[ResidueSelectors]]'),
'simple_metric':('Simple Metric Tag', 'Any of the [[SimpleMetrics]]'),
'jump_selector':('Jump Selector Tag','Any of the [[JumpSelectors]]'),
'res_lvl_task_op':('ResidueLevelTaskOperation Tag','Any of the [[Residue Level TaskOperations]]'),
'res_filter':('ResFilter Tag','Any of the [[ResFilters|OperateOnCertainResiduesOperation#ResFilters]]'),
'features_reporter':('Features Reporter Tag', 'Any of the [[FeatureReporters]]'),
'constraint_generator':('Constraint Generator Tag','Any of the [[ConstraintGenerators]]'),
'pose_selector':('Pose Selectors Tag','Any of the [[Pose Selectors|RosettaScripts-MultiplePoseMover#pose-selectors]]'),
'pose_property_reporter':('Pose Property Reporter Tag','Any of the [[Pose Property Reporters|RosettaScripts-MultiplePoseMover#pose-property-reporters]]'),
'scoring_grid':('Scoring Grid Tag','Any of the [[ScoringGrids|RosettaScripts#scoringgrids]]'),
'rdf_function':('RDF function Tag','Any of the [[RDF Functions|ComputeLigandRDF]]'),

'denovo_perturber':('Perturber Tag','Any of the [[CompoundPerturbers|BuildDeNovoBackboneMover]]'), # Nonce, but added because of recursion
'compound_architect_pairing_group':('CompoundArchitect pairings Tags','Any of the [[CompoundArchitects|BuildDeNovoBackboneMover]]'), # Nonce, but added because of recursion
'denovo_architect':('DenovoArchitect pairings Tags','Any of the [[DenovoArchitects|BuildDeNovoBackboneMover]]'), # Nonce, but added because of recursion
}

ALL_ENTRIES = {}

def alias_type( typename ):
    if typename in TYPE_ALIASES:
        return TYPE_ALIASES[ typename ]
    return typename

def get_annotation(node, name1, name2=None):
    if node.tag != '{http://www.w3.org/2001/XMLSchema}annotation': return ''
    if len(node) < 1 or node[0].tag != '{http://www.w3.org/2001/XMLSchema}documentation':
        if name2 is not None:
            print "Error parsing annotation for entry", name1, "in", name2
        else:
            print "Error parsing annotation for entry", name1
        return ''
    return node[0].text

def parse_subelement( subelem, parentname ):
    #print "PARSING ELEMENT for ", parentname

    main_doc = ''
    tag_lines = []
    doc_lines = []

    if 'name' not in subelem.attrib:
        print 'Error: missing name when parsing subtag for ', parentname
    name = subelem.attrib['name']

    if len(subelem) >=2:
        print "Error: can't parse type for", name, "in", parentname
        return None

    if len( subelem ) == 0:
        if 'type' not in subelem.attrib:
            print "Error: can't parse type entry for", name, "in", parentname
            return None
        typename = subelem.attrib['type']
        if typename in COMMON_TYPES:
            #print 'DOING COMMON', typename , ' in ', parentname, name
            subtag_info, subdoc_line = COMMON_TYPES[ typename ]
            tag_lines.append( '<' + name + ' ... />' )
            doc_lines.append( 'Subtag **' + name +'**: ' + subdoc_line.strip() )
            return '', main_doc, tag_lines, doc_lines

        # Parse complex type in ALL_ENTRIES
        if typename not in ALL_ENTRIES:
            print "Error: can't parse type of ", typename, " for", name, "in", parentname
            return None
        ctype = ALL_ENTRIES[typename]
    elif len( subelem ) == 1:
        ctype = subelem[0]

    if ctype.tag != '{http://www.w3.org/2001/XMLSchema}complexType':
        print "Error getting type for", name, "in", parentname
        return None
    sct = parse_complextype( name, parentname, ctype )
    if sct is None:
        return None
    subtagname, submain_doc, subtag_lines, subdoc_lines = sct
    tag_lines.extend( subtag_lines )

    if len(subdoc_lines) != 0:
        if subtagname:
            doc_lines.append('')
            doc_lines.append('Subtag **' + subtagname + "**:   " + submain_doc.strip() )
            doc_lines.append('')
        doc_lines.extend(subdoc_lines)

    return '', main_doc, tag_lines, doc_lines

def parse_multi( node, parentname ):
    nodetype = node.tag.split('}')[-1]

    #print "PARSING", nodetype.upper(), " for ", parentname

    main_doc = ''
    tag_lines = []
    doc_lines = []

    for subelem in node:
        if subelem.tag == '{http://www.w3.org/2001/XMLSchema}choice':
            sc = parse_choice(subelem , parentname)
            if sc is None:
                continue
            subtagname, submain_doc, subtag_lines, subdoc_lines = sc
            tag_lines.extend( subtag_lines )

            if len(subdoc_lines) != 0:
                if subtagname:
                    doc_lines.append('')
                    doc_lines.append('Subtag **' + subtagname + "**:   " + submain_doc.strip() )
                doc_lines.extend(subdoc_lines)
        elif subelem.tag == '{http://www.w3.org/2001/XMLSchema}element':
            se = parse_subelement(subelem , parentname)
            if se is None:
                continue
            subtagname, submain_doc, subtag_lines, subdoc_lines = se
            tag_lines.extend( subtag_lines )

            if len(subdoc_lines) != 0:
                if subtagname:
                    doc_lines.append('')
                    doc_lines.append('Subtag **' + subtagname + "**:   " + submain_doc.strip() )
                doc_lines.extend(subdoc_lines)
        elif subelem.tag == '{http://www.w3.org/2001/XMLSchema}group':
            sg = parse_group( subelem, parentname )
            if sg is None:
                continue
            subtagname, submain_doc, subtag_lines, subdoc_lines = sg
            tag_lines.extend( subtag_lines )
            doc_lines.extend( subdoc_lines )
        else:
            print "Error parsing subtag of", nodetype, "for", parentname, '::', subelem.tag
            continue


    return '', main_doc, tag_lines, doc_lines

# These used to be separate, but when all said and done I realized they're identical in how we treat the documentation
# I keep the functions separate, in case we want to break them apart at some other time, but for now we'll alias
# them to the same function
parse_choice = parse_multi
parse_sequence = parse_multi
parse_all = parse_multi

def parse_group( node, parentname ):
    #print "PARSING GROUP for ", parentname

    if 'ref' not in node.attrib:
        print "Error parsing group element for ", parentname
        return None

    ref = node.attrib['ref']
    if ref in COMMON_TYPES:
        pseudoname, docline = COMMON_TYPES[ref]
        tagline = '<' + pseudoname + ' ... />'
        docline = '"' + pseudoname + '": ' + docline.strip()
        return '', '', [ tagline, ], ['',docline,]

    if ref not in ALL_ENTRIES:
        print "Error: can't parse group ", typename, " for in", parentname
        return None
    ctype = ALL_ENTRIES[ref]
    if ctype.tag != '{http://www.w3.org/2001/XMLSchema}group' or len(ctype) != 1:
        print "Error: issue parsing group", ref, 'in', parentname
        return None
    ctype = ctype[0]
    if ctype.tag != '{http://www.w3.org/2001/XMLSchema}choice':
        print 'Error: issue parsing group', ref, 'in', parentname
    print '... putting group', ref, 'in-line for', parentname
    sc = parse_choice( ctype, parentname )
    if sc is None:
        return None
    return sc

def parse_complextype( name, parentname, node ):
    main_doc = ''
    attributes = []
    subtags = []

    if parentname is not None:
        desig = parentname + '::' + name
    else:
        desig = name


    for gc in node:
        if gc.tag == '{http://www.w3.org/2001/XMLSchema}annotation':
            main_doc = get_annotation( gc, name ).replace( "*", "\*" )
        elif gc.tag == '{http://www.w3.org/2001/XMLSchema}attribute':
            if 'name' not in gc.attrib:
                if parentname is not None:
                    print "Error parsing attribute for entry", name, 'on', parentname
                else:
                    print "Error parsing attribute for entry", name
                continue
            attrib = dict( gc.attrib )
            docstring = ''
            if 'use' in gc.attrib and gc.attrib['use'] == 'required':
                docstring += '(REQUIRED) '
            for ggc in gc:
                if ggc.tag == '{http://www.w3.org/2001/XMLSchema}annotation':
                    docstring += get_annotation( ggc, gc.attrib['name'], desig ).strip()
            attrib['docstring'] = docstring
            attributes.append( attrib )
        elif gc.tag == '{http://www.w3.org/2001/XMLSchema}choice':
            st = parse_choice(gc, desig)
            if st is not None:
                subtags.append( st )
        elif gc.tag == '{http://www.w3.org/2001/XMLSchema}sequence':
            st = parse_sequence(gc, desig)
            if st is not None:
                subtags.append( st )
        elif gc.tag == '{http://www.w3.org/2001/XMLSchema}group':
            sg = parse_group(gc, desig)
            if sg is not None:
                subtags.append( sg )
        elif gc.tag == '{http://www.w3.org/2001/XMLSchema}all':
            sa = parse_all(gc, desig)
            if sa is not None:
                subtags.append( sa )
        else:
            print "Warning: skipping entry of type '"+ gc.tag+ "' in", name
            pass

    # Format the tags:

    tag_lines = []
    main_line = '<'+name
    #Pull the name subentry out front, if it exists
    for nameattrib in [a for a in attributes if a['name']=='name']:
        main_line += " "+nameattrib['name']+'="(&' + alias_type(nameattrib['type']) + ';)"' #Assume there's no default
    for attrib in attributes:
        if attrib['name'] == "name":
            continue
        attrib_text = attrib['name']+'="('
        if 'default' in attrib:
          attrib_text += attrib['default'] + ' '
        attrib_text += '&' + alias_type(attrib['type']) + ';)"'
        if len(main_line) + len(attrib_text) > 79: # Online Gollum has ~80 char window.
            tag_lines.append(main_line)
            main_line = "       " # 8 spaces, once the extra one gets added.
        main_line += " " + attrib_text
    if len(subtags) == 0:
        main_line += ' />'
    else:
        main_line += ' >'
    tag_lines.append(main_line)

    for subtagname, submain_doc, subtag_lines, subdoc_lines in subtags:
        tag_lines.extend( '    ' + line for line in subtag_lines )
    if len(subtags) != 0 :
        tag_lines.append( '</'+name+'>' )

    #Format the documentation entries

    doc_lines = []
    for attrib in attributes:
        # Skip entry for name
        if attrib['name'] == 'name': continue
        if 'docstring' in attrib and len(attrib['docstring']) != 0:
            doc_lines.append( '-   **' + attrib['name'] + "**: " + attrib['docstring'].strip() )

    for subtagname, submain_doc, subtag_lines, subdoc_lines in subtags:
        if len(subdoc_lines) == 0:
            continue
        doc_lines.append('')

        if subtagname:
            doc_lines.append('Subtag **' + subtagname + "**:   " + submain_doc)
            doc_lines.append('')
        doc_lines.extend(subdoc_lines)

    return name, main_doc, tag_lines, doc_lines

def process( name, node, outfilename ):

    main_doc = ''
    attributes = []
    subtags = []

    if node.tag != '{http://www.w3.org/2001/XMLSchema}complexType':
        print "Error: expected complex type."
        return

    _, main_doc, tag_lines, doc_lines = parse_complextype( name, None, node )

    with codecs.open(outfilename, 'w',encoding='utf8' ) as f:
        f.write( "<!-- THIS IS AN AUTOGENERATED FILE: Don't edit it directly, instead change the schema definition in the code itself. -->\n\n" )

        # TEMPORARY autogenerated labeling
        f.write( '_Autogenerated Tag Syntax Documentation:_\n' )
        f.write( '\n---\n' )

        f.write( main_doc.strip() + '\n\n')
        f.write( '```xml\n' )
        f.write( '\n'.join(tag_lines) )
        f.write( '\n```\n\n' )
        f.write( '\n'.join(doc_lines) )
        f.write( '\n' )

        f.write( '\n---\n' ) # TEMPORARY autogenerated labeling

def main( xsdfile, outdir ):

    xsd = ET.parse( xsdfile )
    root = xsd.getroot()

    if root.tag != '{http://www.w3.org/2001/XMLSchema}schema':
        print "ERROR: Malformed XSD schema file."
        sys.exit()

    entries = {}
    for child in root:
        if 'name' not in child.attrib:
            print "Error in reading", child.tag
            continue

        entries[ child.attrib['name'] ] = child

    global ALL_ENTRIES
    ALL_ENTRIES = entries

    # We do this odd two-stage pass to deal with name conflicts on case-insensitive filesystems
    to_process = {} # dictionary of lower-cased name:(xsd_type_name, display_name)
    for section in SECTIONS:
        if section not in entries:
            print "ERROR: can't find section for", section, " in XSD."
            continue
        sec_entry = entries[ section ]
        sec_choice = None
        for se in sec_entry:
            if se.tag != '{http://www.w3.org/2001/XMLSchema}choice':
                continue
            if sec_choice is not None:
                print "ERROR: malformed entry for section", section
                continue
            sec_choice = se
        if sec_choice is None:
            print "ERROR: malformed entry for section", section
            continue
        for child in sec_choice:
            if 'type' not in child.attrib or 'name' not in child.attrib:
                print "ERROR: malformed entry in section", section
                continue
            entry_type = child.attrib['type']
            entry_name = child.attrib['name']
            if entry_type not in entries:
                print "ERROR: can't find definition for", entry_type
                continue
            to_process.setdefault( entry_type.lower(), [] ).append( (entry_type, entry_name) )

    for single_type in SINGLE_TYPES:
        to_process.setdefault( single_type[0].lower(), [] ).append( single_type )

    for key in to_process:
        to_process[key].sort()
        if len( to_process[key] ) > 1:
            print "Case sensitivity conflict between: ", ' '.join( et for et, en in to_process[key] )
        for ii, (entry_type, entry_name) in enumerate( to_process[key] ):
            if ii == 0:
                outname = outdir + '/' + entry_type + '.md'
            else:
                outname = outdir + '/' + entry_type + "_" + str(ii) + '.md' # Do we need a better way of disambiguating this?
            process( entry_name, entries[entry_type], outname )


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print
        print "ERROR: Usage: ./create_rosetta_scripts_docs.py <xsd file> <directory to write docs to>"
        print
        print __doc__
        sys.exit()

    main(sys.argv[1], sys.argv[2])


