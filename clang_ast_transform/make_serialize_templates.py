"""
This script uses the class data definition geratned by the AST tool
to auto-insert serialization templates into class files.
Templates for the serialize() function will be replaced unless
a template is annotated with "noauto." Separate load/save
templates will not be replaces.

Summary of changes this scritp performs:

1. Insert a serialize() stub into the class header file for all
non-templated classes. For templated classes, private classes and
structs, the full serialize() function is inserted.

2. Create a .tmpl.hh file and insert the implementation of the
serialize() implementation there.

3. For polymorphic classes, insert a registration stub into
the corresponding .cc file.

Notes:

1. Raw pointers and references will NOT be serialized, but
a commented will be added added instead.

2. COPs and CAPs and other const variables will be const_cast'ed.

@author Luki Goldschmidt <lugo@uw.edu>

Run: python3 src/core/io/serialization/make_serialize_templates.py src/core/io/serialization/class_data.def
"""

### TODO: load_and_construct arguments?

import os, sys, re
import blargs

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../" )
#print( sys.path )
import python_cc_reader.code_reader


dryrun = False # Don't change files

########################################################################
# Buffered file

class BufferedFile:
    def __init__(self, filename):
        self.filename = filename
        try:
            self.contents = self.readFile(self.filename)
        except:
            self.contents = ""
        self.relative_filename = filename[ filename.find('src/')+4 : ]
        self.dirty = False
        self.offsets = []

    def save(self):
        if self.dirty:
            self.writeFile(self.filename, self.contents)
            self.dirty = False

    def readFile(self, filename):
        with open(filename, "r") as f:
            return f.read()

    def writeFile(self, filename, contents):
        if dryrun:
            print ("Would write to %s:\n%s" % ( filename, contents ))
            return
        with open(filename, "w") as f:
            return f.write(contents)

    def getOffset(self, position):
        offset = 0
        for o in self.offsets:
            if o[0] <= position:
                offset += o[1]
        return offset

    def getTruePosition(self, position):
        return self.getOffset(position) + position

    def insertStub(self, position, contents, isOrigFilePosition):
        orig_position = position
        if isOrigFilePosition:
            position = self.getTruePosition(position)

        self.contents = self.contents[0:position] + contents + self.contents[position:]
        #print ( "new contents!", self.filename, self.contents )
        self.dirty = True

        self.offsets.append( [ orig_position, len(contents) ] )

    def deleteRange( self, pos1, pos2, isOrigFilePosition ) :
        if isOrigFilePosition :
            pos1 = self.getTruePosition(pos1)
            pos2 = self.getTruePosition(pos2)
        self.contents = self.contents[0:pos1] + self.contents[pos2:]
        self.dirty = True
        self.offsets.append( [ pos1, pos1-pos2 ] ) #APL -- I'm not certain this is correct

########################################################################
# Declaration node types

# the base class for some kind of declaration made in a header file
class Decl:
    def __init__(self, data):
        self.decl, self.name, self.namespace, self.location, self.filerange = data[0:5]
        self.filerange = self.filerange.split("-")
        self.filerange = [ int(self.filerange[0]), int(self.filerange[1]) ]
        self.filename = self.location[0:self.location.find(':')]

# a class or a struct
class RecordDecl(Decl):
    def __init__(self, data):
        super(RecordDecl, self).__init__(data[0:5])
        decl_data = data[5:]
        self.is_templated = decl_data[0] == '1'
        self.is_public = decl_data[1] == '0'
        self.is_polymorphic = decl_data[2] == '1'
        self.base_class_names = None
        self.has_serialization_methods = False
        self.has_save = False
        self.has_load = False
        self.has_load_and_construct = False
        self.fields = []
        self.constructors = []
        self.need_load_construct = False
        self.requires_protected_default_ctor = False
        #self.methods = []
        self.nReferenceVars = 0
        self.nConstVars = 0

# a data member of a class or a struct
class FieldDecl(Decl):
    def __init__(self, data):
        super(FieldDecl, self).__init__(data[0:5])
        decl_data = data[5:]
        self.vartype, self.fullvartype = decl_data
        self.classname = self.namespace
        self.name = self.name.split("::")[-1]

class ConstructorDecl(Decl):
    def __init__(self, data):
        super(ConstructorDecl, self).__init__(data[0:5])
        decl_data = data[5:]
        self.access = int(decl_data[0])
        self.is_public = self.access == 0
        self.is_default_ctor = decl_data[1] == '1'
        self.is_copy_ctor = decl_data[2] == '1'
        self.n_ctor_initializers = int(decl_data[3])
        self.n_min_arguments = int(decl_data[4])
        self.classname = self.namespace
        self.name = self.name.split("::")[-1]

# class MethodDecl(Decl) :
#         def __init__( self, data) :
#                 super(MethodDecl, self).__init__(data[0:5])
#                 decl_data = data[5:]
#                 self.access = int( decl_data[0] )
#                 self.is_public = self.access == 0
#                 self.is_templated = decl_data[1] == "1"
#                 # I don't know if I can obtain this! maybe also parse w/ HeavyCodeReader
#                 self.body_present = decl_data[2]
# 
#                 # all we're interested in is the name of the method
#                 # whether it's public or not, whether it's defined in the
#                 # header or in the cc, and its location in the file
#                 # so that we can add the serialization functions
#                 # beneath it in the header; we'll also take the last
#                 # public function that lives in the .cc and put the
#                 # serialization functions beneath it.
# 

########################################################################
# Translation unit, i.e. a header file containing one or more classes

class Unit:
    def __init__(self, filename):
        self.filename = filename
        self.basename = filename[0:filename.rfind('.')]
        self.hh = BufferedFile(self.basename + ".hh")
        self.srlz = BufferedFile(self.basename + ".srlz.hh")
        self.cc = BufferedFile(self.basename + ".cc")
        self.classes = {}
        self.class_decs_from_cr = {}
        #self.templates = []
        self.class_names = []
        self.polymorphic_types = []
        #self.need_cereal = False
        self.need_srlz_hh = False
        self.needs_access_fwd = False
        self.contains_polymorphic_class = False
        self.default_ctor = []
        self.cc_needs_cereal_base_class_hpp = False
        self.stl_containers = set([])
        self.uses_vector0 = False
        self.uses_vector1 = False
        self.uses_xyzVector = False
        self.parsed_cc_file = None

    def save(self):
        self.hh.save()
        self.srlz.save()
        self.cc.save()

    def preprocess(self):
        """
        Examine all the classes/structs in a single entire translation unit (i.e. a .hh file with one or more classes)
        """
        hcropts = python_cc_reader.code_reader.HCROptions()
        hcropts.defined_preprocessor_variables.append( "SERIALIZATION" )
        classdecs = python_cc_reader.code_reader.read_classes_from_header( self.filename, open( self.filename ).readlines() )
        #print ( "read", len(classdecs), "class in file", self.filename )
        ns = self.filename.rpartition("/")[ 0 ].partition("/source/src/")[2].replace("/","::")+"::"
        for cdec in classdecs :
            #print( ns + cdec.name )
            if (ns + cdec.name) not in self.classes :
                print( "Warning: Missing class:", ns + cdec.name, "from", self.filename, "when read with the HeavyCodReader" )
                continue
            self.class_decs_from_cr[ ns + cdec.name ] = cdec

        for name in self.classes:
            #if name.find("Factory") > 0 or name.find("Creator") > 0 or name.find("Singleton") > 0:
            #    print ("Skipping class %s" % name )
            #    continue
            # print( "For class", name, "we did not find any match to a class in the .hh file; perhaps these?" )
            # if name not in self.class_decs_from_cr :
            #      print( self.class_decs_from_cr.keys() )
            pycc_classdec = self.class_decs_from_cr[name] if name in self.class_decs_from_cr else None            
            self.preprocessObject(self.classes[name], pycc_classdec )
                        
    def addSerializationRoutines( self ) :
        for name in self.classes :
            self.addSerializationRoutinesForClass( self.classes[name] )
        self.save()

    def findSTLTypesInVartype( self, vartype ) :
        if vartype.find( "std::" ) == -1 : return

        container = vartype[ vartype.find( "std::" )+5: vartype.find("<", vartype.find( "std::" )+5) ]
        if container == "" : return
        if container == "shared_ptr" or container == "weak_ptr" :
            self.stl_containers.add( "polymorphic" )
        elif container == "pair" :
            self.stl_containers.add( "utility" )
        elif container == "basic_string" :
            self.stl_containers.add( "string" )
        elif container == "allocator" or container == "less" :
            pass
        else :
            self.stl_containers.add( container )
            #print( self.filename, "uses stl", container )

        remnant = vartype[ vartype.find("<")+1 : vartype.rfind( ">" ) ]
        #print( "remnant:", remnant )
        cols = remnant.split( "," )
        for col in cols :
            self.findSTLTypesInVartype( col )
                
                
    def preprocessObject(self, decl, pycc_decl ):
        """
        Process single object (class, struct, ...) in this unit
        """
        print("preprocessObject:", decl.name)
        # Check if this class already has save or load method
        if pycc_decl :
            for m in pycc_decl.functions :
                if m.name == "save" :
                    decl.has_save = True
                elif m.name == "load" :
                    decl.has_load = True
                elif m.name == "load_and_construct" :
                    decl.has_load_and_construct = True
                if decl.has_save and ( decl.has_load or decl.has_load_and_construct ) :
                    decl.has_serialization_methods = True
    
        if decl.has_serialization_methods :
            # we don't need to add any serialization methods to this class/struct
            print( " decl.has_serialization_methods")
            return

        # Check if this class has reference variables or raw pointers; if so, skip it
        # Check if it has constant variables -- these have to be initialized
        # in a constructor; flag the user to identify which constructor to call
        # in a load-and-construct stub
        for v in decl.fields:
            if v.fullvartype.find("&") >= 0:
                decl.nReferenceVars += 1
            if self.isConstNonPointer(v.fullvartype):
                if v.fullvartype[6:] in [ "double", "int", "float", "bool", "_Bool", "Size" ]:
                    # These can be trivially const_cast'ed
                    continue
                decl.nConstVars += 1

        # Make list of variable to serialize
        decl.membvars = []
        for v in decl.fields:
            vd = { 'enable': True, 'name': v.name, 'vartype': v.vartype, 'cast':'', 'constptr' : False, "fullvartype" : v.fullvartype }
            if self.isRawPtr(v.fullvartype):
                vd['enable'] = False
                vd['comment'] = "%s: %s" % ("reference" if v.vartype.find('&') >= 0 else "raw pointer", v.vartype)
                decl.need_load_construct = True
            elif v.fullvartype.find("ObjexxFCL") >= 0:
                vd['enable'] = False
                vd['comment'] = "ObjexxFCL: %s" % ( v.vartype )
            elif self.isConstNonPointer(v.fullvartype):
                # c = self.getNonConstCast(v.vartype, v.fullvartype)
                # if c != "":
                #     vd['cast'] = c
                # else:
                #     #vd['enable'] = False
                #     vd['comment'] = "const?"
                decl.need_load_construct = True
            elif self.isConstPointer( v.fullvartype ):
                vd['constptr'] = True
            if v.fullvartype.find("std::") >= 0 :
                # for the sake of the #includes we'll need in the .cc file,
                # look at all the stl datatypes used in this 
                self.findSTLTypesInVartype( v.fullvartype )
            if v.fullvartype.find( "xyzVector" ) >= 0 or v.fullvartype.find( "xyzMatrix" ) >= 0 :
                self.uses_xyzVector = True
            if v.fullvartype.find( "utility::vector1" ) >= 0 :
                self.uses_vector1 = True
            if v.fullvartype.find( "utility::vector0" ) >= 0 :
                self.uses_vector0 = True
            decl.membvars.append( vd )


        #decl.simplest_ctor = sorted(decl.constructors, key = lambda x: x.n_min_arguments)[0] if decl.constructors else None # fewest num of arguments
        decl.default_ctor = [ x for x in decl.constructors if x.is_default_ctor ]
        if decl.default_ctor :
            decl.default_ctor = decl.default_ctor[0]
        else :
            decl.default_ctor = None # necessary?

        decl.inline = decl.is_templated #or decl.decl == "struct" # apl modification -- put struct serialization methods into .cc, too
        decl.has_vars = [1 for x in decl.membvars if x['enable'] ] != []
        #decl.need_load_construct = False # decl.decl == "class" and ( decl.nReferenceVars > 0 or decl.nConstVars > 0 )

        # If this class is polymorphic, get base class names
        if decl.is_polymorphic:
            self.contains_polymorphic_class = True
            #decl.base_class_names = []
            #body = self.hh.contents[ self.hh.getTruePosition(decl.filerange[0]) : self.hh.getTruePosition(decl.filerange[1]) ]
            #body = body[ : body.find("{") -1 ]
            #if body.find("#") > 0:
            #    body = body[ : body.find("#") -1 ]
            #p = body.find(":")
            #baseclasses = body[p+1:] if p >= 0 else body
            #for n in baseclasses.split(","):
            #    n = self.strip(n.strip(), ["class", "public", "private", "protected", "virtual"] )
            #    if n.split("::")[-1] == "ReferenceCount": # let's not count ReferenceCount as polymorphic
            #        continue
            #    decl.base_class_names.append( n )
            if pycc_decl and pycc_decl.base is not None :
                self.cc_needs_cereal_base_class_hpp = True
                decl.base_class_names = pycc_decl.base.split(",") if pycc_decl.base.find(",") > 0 else [ pycc_decl.base ]
                #DEBUG
                #for b in decl.base_class_names :
                #    print( " base class for decl.name:", b )

        # Include cereal headers in .hh file if we're inlining things
        if decl.inline:
            self.need_srlz_hh = True
                        
    def addSerializationRoutinesForClass( self, decl ) :
        print("addSerializationRoutinesForClass:", decl.name)

        if decl.has_serialization_methods :
            # we don't need to add any serialization methods to this class/struct
            return

        if decl.inline :
            print(  "Oops! Cannot create serialization functions (now) for class", decl.name )
            return

        # Insert the serialization function declarations into the .hh file
        self.insertHeaderStub( decl )

        # Insert the serialization function implementations into the .cc file
        if not decl.inline:
            self.insertCCStub( decl )

    def isRawPtr(self, s):
        # Guess is this variable is a raw pointer / reference
        return s.find('*') >= 0 or s.find('&') >= 0

    def getContainedType(self, s):
        # Guess type contained in an STL container. Not super robust...
        if s.find('<') < 0:
            return s
        container = s[ 0 : s.find('<') ]
        contained = s[ s.find('<')+1 : s.rfind('>') -1 ]
        container = container.replace("class ", "")
        if container in [ "std::vector", "utility::pointer::vector0", "utility::pointer::vector1" ]:
            if contained.find(',') > 0:
                contained = contained[ 0 : contained.find(',') ]
        if container in [ "std::map" ]:
            p = contained.find(',', contained.find(',')+1)
            if p > 0:
                contained = contained[ 0 : p ]
        return contained

    def getContainerType( self, s ):
        if s.find( '<' ) < 0 :
            return ""
        return s[ 0 : s.find( '<' ) ]
        

    def isConstNonPointer(self, s):
        # Guess if this is a const variable
        return s[0:6] == 'const '

    def isConstPointer( self, s ) :
        container = self.getContainerType(s)
        #print( "isConstPointer container:", container )
        if container == "" :
            #print( "isConstPointer container:", False )
            return False
        if container == "class std::shared_ptr" or container == "class std::access_ptr" :
            contained = self.getContainedType(s)
            #print( "isConstPointer contained:", contained, contained.find('const ') >= 0 )
            return contained.find('const ') >= 0
        else :
            return self.isConstPointer( self.getContainedType(s) )

    def getNonConstCast(self, short, full):
        # Guess what the non-const equivalent would be
        if short[-3:] == "COP" or short[-3:] == "CAP":
            return ""

        c = short
        c = c.replace("const ", "")
        c = c.strip()

        if c != short:
            return "const_cast< %s & >( @@ )" % ( c )

        #c = full.replace("const ", "").strip()
        #if c != full:
        #    return "const_cast< %s & >( @@ )" % ( c )
        return ""
    def remove_trailing_underscore( self, varname ) :
        if varname[-1] == "_" :
            return varname[:-1]
        return varname

    def excise_templated_type( self, fulltype, templated_type ) :
        i = 0
        while i < len(fulltype) :
            i = fulltype.find( ", " + templated_type )
            if i == -1 : break
            bracket_count = 1
            j = i + len( ", " + templated_type + "<" )
            #print( "excising allocators:", fulltype, i, j )
            while j < len(fulltype) :
                #print ( " while", j, "fulltype[j]", fulltype[j] )
                if fulltype[j] == ">" :
                    bracket_count -= 1
                    if bracket_count == 0 :
                        j+=1
                        break
                elif fulltype[j] == "<" :
                    bracket_count += 1
                j+=1
            
            fulltype = fulltype[:i] + fulltype[j:]
            #print( "excised allocator:", fulltype, i )
        return fulltype

    def nonconstptr_type( self, fulltype ) :
        fulltype = fulltype.replace( "const ", "" )
        fulltype = fulltype.replace( "class ", "" )
        fulltype = fulltype.replace( "<", "< " )
        fulltype = self.excise_templated_type( fulltype, "std::allocator" )
        fulltype = self.excise_templated_type( fulltype, "std::less" )
        fulltype = fulltype.replace( ">", " >" )
        fulltype = fulltype.replace( "  >", " >" )
        print ("final type", fulltype )
        return fulltype

    def makeVarsArStub( self, decl, membvars, style, include_base = True ):
        """
        Make the contents of the serialization function, serializing all member variables
        """
        if not decl.has_vars and not decl.base_class_names :
            return ""

        stub = []
        if decl.base_class_names and style != "load_and_construct" :
            for base in decl.base_class_names :
                if base == "utility::pointer::ReferenceCount" : continue
                # assume no virtual inheritance
                stub.append( "arc( cereal::base_class< %s >( this ) );\n" % base )

        if style == "load_and_construct" :
            stub.append( "NOTE: The automatically generated load_and_construct stub here requires manual intervention\n" );
            stub.append( "construct( /* please choose the proper constructor arguments */ );\n" );

        for i, v in enumerate(membvars):
            nextline = []
            var_name = v['name']
            deserialize_varname = var_name
            if style == "load_and_construct":
                var_name = "construct->" + var_name
                deserialize_varname = var_name
            # cast_str = v['cast'].replace("@@", var_name) if v.get('cast') else ""

            if style == "load" and v[ 'constptr' ] :
                #print( "CONSTPTR:", var_name)
                deserialize_varname = "local_%s" % self.remove_trailing_underscore( v[ 'name' ] )
                stub.append( "%s %s;\n" % ( self.nonconstptr_type( v[ 'fullvartype' ] ), deserialize_varname ));

            if not v['enable']:
                nextline.append( "// " )
            nextline.append( "arc( " )

            if style == "serialize":
                nextline.append( "CEREAL_NVP( %s )" % v['name'] )
            elif style == "save":
                nextline.append( "CEREAL_NVP( %s )" % v['name'] )
            elif style == "load":
                nextline.append( deserialize_varname )
            elif style == "load_and_construct":
                nextline.append( deserialize_varname )

            # Add comma at end of line
            #if [ 1 for x in vars[i+1:] if x['enable'] ]:
            #    nextline.append( "" )
            nextline.append( " );" )

            # Add comment
            nextline.append( " // " + v['vartype'] )
            if v.get('comment'):
                nextline.append(  "; " + v['comment'] )
            nextline.append( "\n" )
            stub.append( "".join( nextline ))
            if v['constptr'] and ( style == "load" or style == "load_and_construct" ) :
                stub.append( "%s = %s; // copy the non-const pointer(s) into the const pointer(s)\n" % ( var_name, deserialize_varname ))

        return "".join( stub )

    #def makeFullStub( self, decl, style ):
    #    """
    #    Make the full stub for the serialize, save, load or load_and_construct functions
    #    """
    #    stub = ""
    #    if style in [ "serialize", "save", "load" ]:
    #        stub += self.makeVarsArStub( decl, decl.membvars, style )
    #    elif style == "load_and_construct":
    #        # Assume the first class n vars are for the c'tor... ouch
    #        ctor_vars = []
    #        more_vars = decl.membvars
    #
    #        """
    #        # Try to use non-default c'tor
    #        n_min_arguments = decl.simplest_ctor.n_min_arguments if decl.simplest_ctor else 0
    #        for v in decl.membvars[ 0 : n_min_arguments ]:
    #            vartype = self.strip(v['vartype'], [ "enum", "const", "class", "struct" ]).rstrip("& ")
    #            name = self.strip( v['name'].rstrip("_"), [ "enum", "const", "class", "struct" ] )
    #            stub += "%s %s;\n" % ( vartype, name )
    #            vc = v.copy()
    #            vc['name'] = name;
    #            ctor_vars.append( vc )
    #            more_vars = more_vars[ 1: ]
    #        """
    #        if ctor_vars:
    #            # TODO: how do you use cereal::base_class< T > in static load_and_construct()?
    #            stub += self.makeVarsArStub( decl, ctor_vars, "load", False     ) + "\n"
    #            args = ", ".join( [ v['name'] for v in ctor_vars ] )
    #            stub += "%s o( %s );\n\n" % ( decl.name.split(":")[-1], args )
    #        else:
    #            stub += "%s o;\n" % ( decl.name.split(":")[-1] )
    #
    #        stub += self.makeVarsArStub( decl, more_vars, style, False )
    #        stub += "construct( o );\n"
    #
    #    return stub

    def indent( self, s, i ):
        if not s: return s
        # Indent lines
        i = "\t" * i
        return "\n".join( [ i + l for l in s.splitlines() ] ) + "\n"

    def strip( self, s, tags ):
        # Strip tags at the beginning of string s
        old_s = ""
        while old_s != s:
            s = s.strip()
            old_s = s
            for t in tags:
                if s[ 0:len(t) ] == t:
                    s = s[ len(t): ]
        return s.strip()

    def insertHeaderStub( self, decl ):
        """
        Write the prototype that goes into the .hh file
        """
        stub = []

        if True :
        #if decl.has_vars or not decl.is_polymorphic:
            stub_body = self.indent(self.makeVarsArStub(decl, decl.membvars, 'save'), 2) if decl.inline else ""
            stub.append( "\ttemplate< class Archive > void save( Archive & arc ) const" )
            if decl.inline:
                stub.append( " {\n" + stub_body + "\t}\n" if stub_body else " {}\n" )
            else:
                stub.append( ";\n" )
            if decl.need_load_construct:

                stub_body = self.indent(self.makeVarsArStub(decl, decl.membvars, 'load_and_construct'), 2) if decl.inline else ""
                stub.append( "\ttemplate< class Archive > static void load_and_construct( Archive & arc, cereal::construct< %s > & construct )" \
                    % decl.name.split("::")[-1] )
                if decl.inline:
                    stub.append( " {\n" + stub_body + "\t}\n" if stub_body else " {}\n" )
                else:
                    stub.append( ";\n" )
            else:
                stub_body = self.indent(self.makeVarsArStub(decl, decl.membvars, 'load'), 2) if decl.inline else ""
                stub.append( "\ttemplate< class Archive > void load( Archive & arc )" )
                if decl.inline:
                    stub.append( " {\n" + stub_body + "\t}\n" if stub_body else " {}\n" )
                else:
                    stub.append( ";\n" )

        if not stub:
            return "";

        full_stub= [ "\n#ifdef    SERIALIZATION\n" ]
        if decl.decl == "class":
            if not decl.default_ctor and decl.constructors: # Add default c'tor prototype
                decl.requires_protected_default_ctor = True
                self.needs_access_fwd = True
                full_stub.append( "protected:\n\tfriend class cereal::access;\n\t%s();\n\n" % ( decl.name.split("::")[-1] ) )
                self.default_ctor.append( decl.name )
            elif decl.default_ctor and not decl.default_ctor.is_public:
                self.needs_access_fwd = True
                full_stub.append( "\tfriend class cereal::access;\n" )
            full_stub.append( "public:\n" )
        full_stub +=  stub
        full_stub.append( "#endif // SERIALIZATION\n" )

        self.hh.insertStub( decl.filerange[1] -1, "".join(full_stub), True ) # end of class/struct

    def insertCCStub(self, decl ):
        """
        Write the contents that go into the .cc file for a single class
        """
        if not decl.has_vars and not decl.is_polymorphic and not decl.base_class_names:
            print ("insertCCStub early exit")
            return

        stub = []
        stub.append("\n#ifdef    SERIALIZATION\n\n")
        if decl.requires_protected_default_ctor :
            stub.append( "/// @brief Default constructor required by cereal to deserialize this class\n" )
            stub.append( "%s::%s() {}\n\n" % ( decl.name, decl.name.split("::")[-1] ) )

        stub.append("/// @brief Automatically generated serialization method\n" )
        stub.append("template< class Archive >\n")
        stub.append("void\n")
        stub.append("%s::save( Archive &%s ) const {\n" % ( decl.name, " arc" if decl.has_vars or decl.base_class_names else "" ) )
        stub.append(self.indent( self.makeVarsArStub( decl, decl.membvars, "save" ), 1 ))
        stub.append("}\n\n")

        if decl.need_load_construct:
            stub.append("/// @brief Automatically generated deserialization method\n" )
            stub.append("template< class Archive >\n")
            stub.append("void\n")
            stub.append("%s::load_and_construct( Archive &%s, cereal::construct< %s > & construct ) {\n" % \
                        ( decl.name, " arc" if decl.has_vars or decl.base_class_names else "", decl.name ))
            stub.append(self.indent( self.makeVarsArStub( decl, decl.membvars, "load_and_construct" ), 1 ))
            stub.append("}\n")
        else:
            stub.append("/// @brief Automatically generated deserialization method\n" )
            stub.append("template< class Archive >\n")
            stub.append("void\n")
            stub.append("%s::load( Archive &%s ) {\n" % ( decl.name, " arc" if decl.has_vars or decl.base_class_names else "" ))
            stub.append(self.indent( self.makeVarsArStub( decl, decl.membvars, "load" ), 1 ))
            stub.append("}\n\n")

        if decl.need_load_construct:
            stub.append( "SAVE_AND_LOAD_AND_CONSTRUCT_SERIALIZABLE( %s );\n" % decl.name )
        else :
            stub.append( "SAVE_AND_LOAD_SERIALIZABLE( %s );\n" % decl.name );

        if decl.is_polymorphic :
            stub.append( "CEREAL_REGISTER_TYPE( %s )\n\n" % decl.name );
        stub.append( "#endif // SERIALIZATION\n" )
        stub = "".join(stub)

        #print("stub for", decl.name,"\n",stub,"inserted at",len(self.cc.contents))
        self.cc.insertStub( len(self.cc.contents), stub, False )

    def save(self):
        self.saveHH()
        #self.saveTmplHH()
        self.saveCC()

    def saveHH(self):
        # Make final modifications to the .hh file before outputting it

        # the #includes that are needed for the header
        stub = []
        if self.needs_access_fwd or self.contains_polymorphic_class :
            stub.append( "\n#ifdef    SERIALIZATION\n" )
            stub.append( "// Cereal headers\n" )
            if self.needs_access_fwd :
                stub.append( "#include <cereal/access.fwd.hpp>\n" )
            if self.contains_polymorphic_class :
                # need to #include polymorphic.fwd.hpp for the CEREAL_FORCE_DYNAMIC_INIT macro
                stub.append( "#include <cereal/types/polymorphic.fwd.hpp>\n" )
            stub.append( "#endif // SERIALIZATION\n\n" )

            # insert this #include before the first "namespace ... {"
            p = self.hh.contents.find("\nnamespace")
            self.hh.insertStub( p, "".join( stub ), False )

        # for dealing with dynamically linked libraries
        if self.contains_polymorphic_class :
            stub = ["\n"]
            stub.append( "#ifdef    SERIALIZATION\n" )
            stub.append( "CEREAL_FORCE_DYNAMIC_INIT( %s )\n" % self.hh.filename.partition("main/source/src/")[2].replace( "/", "_" )[:-3] )
            stub.append( "#endif // SERIALIZATION\n\n" )
            self.hh.insertStub( self.hh.contents.rfind("\n#endif"), "".join(stub), False )

        self.merge_serialization_blocks( self.hh )

        self.hh.save()

    def saveCC(self):
        """
        Make final modifications to the .cc file before outputting it
        """
        if not self.cc.contents:
            # Skip this file... should go to general .cc file
            print( "save CC early return, no contents" )
            return

        stub = [ "\n#ifdef SERIALIZATION\n" ]
        stub.append( "// Utility serialization headers\n" ) 
        stub.append( "#include <utility/serialization/serialization.hh>\n" )
        if self.uses_vector0 :
            stub.append( "#include <utility/vector0.srlz.hh>\n" )
        if self.uses_vector1 :
            stub.append( "#include <utility/vector1.srlz.hh>\n" )
        stub.append( "\n" )
        if self.uses_xyzVector :
            stub.append( "// Numeric serialization headers\n" )
            stub.append( "#include <numeric/xyz.serialization.hh>\n" )
            stub.append( "\n" )
        if self.cc_needs_cereal_base_class_hpp or self.stl_containers or self.needs_access_fwd:
            stub.append( "// Cereal headers\n" )
        includes = []
        if self.contains_polymorphic_class and "polymorphic" not in self.stl_containers :
            includes.append( "cereal/types/polymorphic.hpp" )        
        if self.needs_access_fwd :
            includes.append( "cereal/access.hpp" )
        if self.cc_needs_cereal_base_class_hpp :
            includes.append( "cereal/types/base_class.hpp" )
        for stltype in self.stl_containers :
            includes.append( "cereal/types/%s.hpp" % stltype )
        includes = sorted( includes )
        for include in includes :
            stub.append( "#include <%s>\n" % include )
        stub.append( "#endif // SERIALIZATION\n\n" )

        p = self.cc.contents.find("\nnamespace")
        self.cc.insertStub( p, "".join(stub), False ) # before first namespace

        # for dealing with dynamically linked libraries
        if self.contains_polymorphic_class :
            stub = []
            stub.append( "\n#ifdef    SERIALIZATION\n" )
            stub.append( "CEREAL_REGISTER_DYNAMIC_INIT( %s )\n"% self.hh.filename.partition("main/source/src/")[2].replace( "/", "_" )[:-3] )
            stub.append( "#endif // SERIALIZATION\n" )
            self.cc.insertStub( len(self.cc.contents), "".join(stub), False )
    
        self.merge_serialization_blocks( self.cc )

        self.cc.save()

    def saveTmplHH(self):
        """
        Create / update the .tmpl.hh file, if needed
        """
        if not self.templates:
            return

        guard_name = "INCLUDED_SERIALIZATION_" + self.tmpl.relative_filename.replace("/", "_").replace(".", "_")

        copyright_stub = """
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   %s
/// @brief  Auto-generated serialization templates
""" \
        % ( self.tmpl.relative_filename )

        stub = """
#ifdef SERIALIZATION
#ifndef %s
#define %s

#include <core/io/serialization/types.hh>

// Serialization templates
%s

#endif // %s
#endif // SERIALIZATION
""" \
        % ( guard_name, guard_name, "\n".join(self.templates), guard_name )

        if self.tmpl.contents:
            self.tmpl.contents += stub
        else:
            self.tmpl.contents = copyright_stub.lstrip() + stub

        self.tmpl.dirty = True # hack
        self.tmpl.save()

    def merge_serialization_blocks( self, buff ) :
        counter = 0
        while buff.contents.find( "#endif // SERIALIZATION", counter ) >= 0 :
            start = buff.contents.find( "\n#endif // SERIALIZATION", counter )
            end   = buff.contents.find( "\n#ifdef    SERIALIZATION", start )
            if end == -1 : break
            #print( "merge_serialization_blocks: Found start end pair:", start, end )

            # now lets check that there's nothing between start and end
            middle = buff.contents[ (start+len("#endif // SERIALIZATION")+1):end ]
            not_empty = False
            for i in range(len(middle)) :
                if middle[i] == " " : continue
                if middle[i] == "\t" : continue
                if middle[i] == "\n" : continue
                not_empty = True
                break
            if not_empty :
                #print( "not empty:", middle )
                counter = end
                continue
            #print( "EMPTY!", middle )

            # we're clear to delete the #endif..#ifdef block
            buff.deleteRange( start, end + len("\n#ifdef    SERIALIZATION"), False )
            counter = start

def processDef(def_filename, hhs):
    """
    Process and insert all stubs for the entire definition (many classes)
    """
    units = {}

    # Load definition data
    for line in open(def_filename).readlines():
        line = line.strip().split("\t")
        decl = line[0]

        if decl in ["class", "struct", "union"]:
            d = RecordDecl( line )
            if d.filename not in units:
                units[d.filename] = Unit( d.filename )
            units[ d.filename ].classes[ d.name ] = d
        elif decl in ["constructor"]:
            d = ConstructorDecl( line )
            if d.filename in units and d.classname in units[ d.filename ].classes:
                if d.filerange[1] - d.filerange[0] > 0:
                    units[ d.filename ].classes[ d.classname ].constructors.append(d)
        elif decl in ["field"]:
            d = FieldDecl( line )
            if d.filename in units:
                units[ d.filename ].classes[ d.classname ].fields.append(d)

    # Process all units
    for unit_name in hhs: #sorted(units.keys()) :
        if unit_name[0] != "/" :
            unit_name = os.getcwd()+ "/" + unit_name
        if unit_name not in units :
            print( "No classes in file", unit_name, "listed in the definition file; skipping")
            continue
        #print( "preprocessing", unit_name )
        units[ unit_name ].preprocess()

    # rewrite a portion of them
    for unit_name in hhs:
        if unit_name[0] != "/" :
            unit_name = os.getcwd()+ "/" + unit_name
        if unit_name not in units :
            print( "No classes in file", unit_name, "listed in the definition file; skipping")
            continue
        print ("Processing unit %s" % unit_name)
        unit = units[unit_name]
        unit.addSerializationRoutines()

if __name__ == '__main__':
    with blargs.Parser(locals()) as p :
        p.str( "definitions" ).shorthand( "d" ).required().described_as( "The file with the concatenated output of the class names and fields" )
        p.multiword( "hhs" ).required().cast( lambda x : x.split() ).described_as( "The list of header files to proccess" )
    processDef(definitions, hhs)