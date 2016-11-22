#!/usr/bin/env python

import sys, os
import xml.etree.ElementTree as ET
import codecs

TYPE_ALIASES={ 'xs:string':'string', }
SECTIONS = [ 'mover', ] # ['loop_definer', 'residue_selector', 'scoring_grid', 'task_operation', 'layer_design_ss_layer', 'layer_design_ss_layer_or_taskop', 'res_lvl_task_op', 'res_filter', 'mover', 'constraint_generator', 'denovo_architect', 'compound_architect_pairing_group', 'denovo_perturber', 'denovo_folder', 'rdf_function', 'pose_selector', 'pose_property_reporter', 'filter', 'features_reporter', 'envclaim', 'scriptcm']


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


def parse_choice( subtag, parentname ):
    if 'minOccurs' in subtag.attrib and subtag.attrib['minOccurs'] == '0':
        optional = True
    else:
        optional = False

    return None

def parse_complextype( name, parentname, node ):
    main_doc = ''
    attributes = []
    subtags = []

    for gc in node:
        if gc.tag == '{http://www.w3.org/2001/XMLSchema}annotation':
            main_doc = get_annotation( gc, name )
        elif gc.tag == '{http://www.w3.org/2001/XMLSchema}attribute':
            if 'name' not in gc.attrib:
                if parentname is not None:
                    print "Error parsing attribute for entry", name, 'on', parentname
                else:
                    print "Error parsing attribute for entry", name
                continue
            attrib = dict( gc.attrib )
            for ggc in gc:
                if ggc.tag == '{http://www.w3.org/2001/XMLSchema}annotation':
                    attrib['docstring'] = get_annotation( ggc, gc.attrib['name'], name )
            attributes.append( attrib )
        elif gc.tag == '{http://www.w3.org/2001/XMLSchema}choice':
            st = parse_choice(gc, name)
            if st is not None:
                subtags.append( st )
        else:
            #print "Warning: skipping entry of type '"+ gc.tag+ "' in", name
            pass

    # Format the tags:

    tag_lines = []
    main_line = '<'+name
    #Pull the name subentry out front, if it exists
    for nameattrib in [a for a in attributes if a['name']=='name']:
        main_line += " "+nameattrib['name']+'="(' + alias_type(nameattrib['type']) + ')"'
    for attrib in attributes:
        if attrib['name'] == "name":
            continue
        main_line += " "+attrib['name']+'="(' + alias_type(attrib['type'])
        if 'default' in attrib:
          main_line += "; " + attrib['default']
        main_line += ')"'
    if len(subtags) == 0:
        main_line += ' />'
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
            doc_lines.append( '-   ' + attrib['name'] + ": " + attrib['docstring'].strip() )

    for subtagname, submain_doc, subtag_lines, subdoc_lines in subtags:
        if len(doc_lines) == 0:
            continue
        doc_lines.append('')

        doc_lines.append('For subtag ' + subtagname + ": " + submain_doc)
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
        f.write( main_doc.strip() + '\n\n')
        f.write( '```\n' )
        f.write( '\n'.join(tag_lines) )
        f.write( '\n```\n\n' )
        f.write( '\n'.join(doc_lines) )
        f.write( '\n' )

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

    for section in SECTIONS:
        if section not in entries:
            print "ERROR: can't find section for", section, " in XSD."
            continue
        sec_entry = entries[ section ]
        if len(sec_entry) != 1 or sec_entry[0].tag != '{http://www.w3.org/2001/XMLSchema}choice':
            print "ERROR: malformed entry for section", section
            continue
        sec_entry = sec_entry[0]
        for child in sec_entry:
            if 'type' not in child.attrib or 'name' not in child.attrib:
                print "ERROR: malformed entry in section", section
                continue
            entry_type = child.attrib['type']
            entry_name = child.attrib['name']
            if entry_type not in entries:
                print "ERROR: can't find definition for", entry_type
                continue
            process( entry_name, entries[entry_type], outdir + '/' + entry_type + '.md' )

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print "ERROR: Usage: ./create_rosetta_scripts_docs.py <xsd file> <directory to write docs to>"
        sys.exit()

    main(sys.argv[1], sys.argv[2])


