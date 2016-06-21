import willclang;
import sys,os;
from clang.cindex import Index,CursorKind,TypeKind;


class ClassDeclaration :
  def __init__( self ) :
    self.name = ""
    self.ns = ""
    self.parents = []
    self.file_declared_within = ""
    self.cursor_hash = None
    self.depth_declared_at = 0
    self.templated = False

class CodeReader :
  def __init__( self, ast ) :
    self.ast = ast
    self.classes = {} # dictionary from cursor-hash to ClassDeclarationObject
    self.last_class_processing = None

  def addclass( self, node ) :
    newclass = ClassDeclaration()
    newclass.name = node.cursor.spelling
    newclass.ns = node.namespace()
    newclass.cursor_hash = node.cursor.get_usr()
    newclass.file_declared_within = node.cursor.location.file
    newclass.depth_declared_at = node.depth
    self.classes[ newclass.cursor_hash ] = newclass
    self.last_class_processing = newclass

  def class_derives_from_refcount( self, cl ) :
    for parent in cl.parents :
      parentcl = self.classes[ parent ]
      if parentcl.name == "ReferenceCount" and parentcl.ns == "utility::pointer" :
        return True
      else :
        if self.class_derives_from_refcount( parentcl ) :
          return True
    return False


  def note_classes_deriving_from_reference_count( self ) :
    self.classes_deriving_from_refcount = {}
    for usr, cl in self.classes.items() :
      if self.class_derives_from_refcount( cl ) :
        self.classes_deriving_from_refcount[ usr ] = cl

  def find_class_declarations( self, node ) :
    #print "fcd: visiting", node.cursor.kind, node.cursor.location.file
    if self.last_class_processing :
      if node.depth <= self.last_class_processing.depth_declared_at :
        self.last_class_processing = None
    if node.cursor.kind == CursorKind.CLASS_DECL or node.cursor.kind == CursorKind.STRUCT_DECL :
      self.addclass( node )
      return "class_decl %s %s %s" % (node.cursor.kind, node.cursor.kind.is_declaration(), node.cursor.spelling)
    elif node.cursor.kind == CursorKind.CLASS_TEMPLATE :
      self.addclass( node )
      self.last_class_processing.templated = True
      return "class_temp_decl " + node.cursor.spelling
    elif node.cursor.kind == CursorKind.CXX_BASE_SPECIFIER :
      if not self.last_class_processing : return
      parent_usr = node.cursor.type.get_declaration().get_usr()
      # only keep parent classes that I have correctly captured as being classes
      if not parent_usr in self.classes : return
      self.last_class_processing.parents.append( parent_usr )
      #return "cxx_base_specifier %s %s %s %s spelling: %s namespace: %s" % (node.cursor.kind, node.cursor.kind.is_reference(), node.cursor.type.kind, node.cursor.type.get_declaration(), node.cursor.type.get_declaration().spelling, node.definition().namespace() )

  def visit( self, node ) :
    #print "visiting", node.cursor.get_usr()
    str = ""
    if node.cursor.kind == CursorKind.BINARY_OPERATOR :
      str = "binary_operator %s %s %s " % (node.cursor.kind, node.cursor.type.kind, binary_operator_token( node ) )
    elif node.cursor.kind == CursorKind.CALL_EXPR :
      str = "call expression: %s %s %s " % ( node.cursor.kind, node.cursor.type.kind, node.cursor.get_definition() )
    else :
      str =  "cursorkind= %s typekind=%s isdec=%s isref=%s" % (node.cursor.kind, node.cursor.type.kind, node.cursor.kind.is_declaration(), node.cursor.kind.is_reference())
      if node.definition() :
        str += " %s" % node.definition().cursor.kind
      str += " %s" % node.getstr(recursive=False)
    str += "\n"
    #for token in node.cursor.get_tokens() :
    #  str += " " * 2 * node.depth + "tok:(" + token.spelling + ")\n"
    return str

class CodeQualityChecker_FindRawPtrs_to_RefcountSubclasses :
  def __init__( self ) :
    self.problem_string = ""

  def has_problem( self ) :
    return self.problem_string != ""

  def problems( self ) :
    return self.problem_string

  def examine_ast( self, ast ) :
    cr = CodeReader(ast)
    ast.root.treemap( cr.find_class_declarations, allchild=True )
    cr.note_classes_deriving_from_reference_count()
    rfptrfinder = VarDecToRefCountSubclassPointerFinder( cr )
    ast.root.treemap( rfptrfinder.save_lines_w_rawptrs_assigned_opdata, allchild=False )
    if rfptrfinder.lines_w_rawpts_assigned_opdata :
      self.problem_string = "\n".join( rfptrfinder.lines_w_rawpts_assigned_opdata )

class VarDecToRefCountSubclassPointerFinder :
  def __init__( self, codereader ) :
    self.decs_to_refcount_subclass_pointers = {}
    self.code_reader = codereader
    self.lines_w_rawpts_assigned_opdata = []

  def print_lines_w_rawptrs_assigned_opdata( self, node ) :
    """
    Print summary information from nodes in the AST where the contents of
    an owning pointer is assigned to a raw pointer
    """
    if self.vardec_to_pointer_to_refcount_deriving_class( node ) or self.op_assignment_to_rawptr( node ) :
      print node.cursor.location.file, node.cursor.location.line
      print ">>> ",
      for token in node.cursor.get_tokens() :
        print token.spelling,
      print

  def save_lines_w_rawptrs_assigned_opdata( self, node ) :
    """
    Save a record of this node if it represents an owning pointer's contents being assigned to a raw pointer
    """
    if self.vardec_to_pointer_to_refcount_deriving_class( node ) or self.op_assignment_to_rawptr( node ) :
      saveout1 = "%s %s\n" % ( node.cursor.location.file, node.cursor.location.line )
      saveout2 = ">>> "
      for token in node.cursor.get_tokens() :
        saveout2 += token.spelling
      saveout2 += "\n"
      self.lines_w_rawpts_assigned_opdata.append( saveout1 )
      self.lines_w_rawpts_assigned_opdata.append( saveout2 )

  def vardec_to_pointer_to_refcount_deriving_class( self, node ) :
    """
    Does the given node represent the declaration of a raw pointer to a class derived from ReferenceCount?
    Requires that note_classes_deriving_from_reference_count() has been called
    """
    if node.cursor.kind != CursorKind.VAR_DECL : return False
    if node.cursor.type.kind != TypeKind.POINTER : return False
    return self.first_typeref_child_is_to_refcount_derived_class( node )

  def first_typeref_child_is_to_refcount_derived_class( self, node ) :
    """
    Look at the children of this node.  If the first non-namespace-ref child
    is a TYPE_REF to a class derived from ReferenceCount, return True
    """
    for child in node.allchild :
      #print child.cursor.kind,
      #if child.definition() :
      #  print child.definition().cursor.spelling
      #else :
      #  print
      if child.cursor.kind == CursorKind.NAMESPACE_REF :
        continue
      if child.cursor.kind == CursorKind.TYPE_REF :
        childdef = child.definition()
        if not childdef : return False
        defcursor = childdef.cursor
        return defcursor.get_usr() in self.code_reader.classes_deriving_from_refcount
      #print "ERROR: Non TYPE_REF, non NAMESPACE_REF cursor encountered in first_typeref_child_is_to_refcount_derived_class:", node.cursor.kind
      return False

  def op_assignment_to_rawptr( self, node ) :
    if self.vardec_to_pointer_to_refcount_deriving_class( node ) :
      if not node.allchild : return False
      return ptr_to_opheld_obj( node.allchild[-1] )
    if node.cursor.kind == CursorKind.BINARY_OPERATOR :
      bintok = binary_operator_token( node )
      if bintok == "=" :
        if not node.allchild : return False
        if len(node.allchild) != 2 : return False
        if not self.node_returns_ptr_to_refcount_subclass( node.allchild[0] ) : return False
        return ptr_to_opheld_obj( node.allchild[1] )
    return False

  def node_returns_ptr_to_refcount_subclass( self, node ) :
    """
    Recursive routine that returns True when a node returns a pointer to a class derived from refcount
    Handles variables that have previously been declared, functions that return pointers to classes
    derived from refcount, array accesses, and any of the previous three wrapped in parenthesis.
    """
    if node.cursor.kind == CursorKind.DECL_REF_EXPR :
      return self.decl_ref_expr_to_refcount_subclass( node )
    elif node.cursor.kind == CursorKind.CALL_EXPR :
      return self.function_returns_ptr_to_refcount_subclass( node )
    elif node.cursor.kind == CursorKind.ARRAY_SUBSCRIPT_EXPR :
      return self.array_index_returns_ptr_to_refcount_subclass( node )
    if node.cursor.kind == CursorKind.PAREN_EXPR :
      return self.node_returns_ptr_to_refcount_subclass( node.allchild[0] )
    return False

  def decl_ref_expr_to_refcount_subclass( self, node ) :
    """
    Returns true for nodes which represent variables that have been previously declared
    as pointers to classes deriving from ReferenceCount
    """
    if node.cursor.kind != CursorKind.DECL_REF_EXPR : return False
    nodedef = node.definition()
    if not nodedef : return False
    return self.vardec_to_pointer_to_refcount_deriving_class( nodedef )


  def function_returns_ptr_to_refcount_subclass( self, node ) :
    """
    Returns true for nodes which represent function calls where the function returns a
    pointer to a class derived from ReferenceCount
    """
    if node.cursor.kind == CursorKind.CALL_EXPR and node.cursor.type.kind == TypeKind.POINTER :
      nodedef = node.definition()
      if not nodedef : return False
      if not nodedef.allchild : return False
      nodedefchild0def = nodedef.allchild[0].definition()
      if not nodedefchild0def : return False
      return nodedefchild0def.cursor.get_usr() in self.code_reader.classes_deriving_from_refcount
    return False

  def array_index_returns_ptr_to_refcount_subclass( self, node ) :
    """
    Returns True for nodes which represent access to an array of pointers
    to a class derived from ReferenceCount
    """
    if node.cursor.kind != CursorKind.ARRAY_SUBSCRIPT_EXPR or node.cursor.type.kind != TypeKind.POINTER : return False
    if not node.allchild : return False
    child0def = node.allchild[0].definition()
    if not child0def : return False
    #print "child0def"
    #child0def.treeprint( lambda x : "%s %s" % ( x.cursor.kind, x.cursor.type.kind ) )
    if child0def.cursor.kind != CursorKind.VAR_DECL : return False
    if not child0def.allchild : return False
    child0defchild0 = child0def.allchild[0]
    if child0defchild0.cursor.kind != CursorKind.TYPE_REF or child0defchild0.cursor.type.kind != TypeKind.RECORD : return False
    child0defchild0def = child0defchild0.definition()
    if not child0defchild0def : return False
    #print "child0defchild0def"
    #child0defchild0def.treeprint( lambda x : "%s %s" % ( x.cursor.kind, x.cursor.type.kind ) )
    return child0defchild0def.cursor.get_usr() in self.code_reader.classes_deriving_from_refcount


  def initialization_source_is_owning_ptr( self, node ) :
    """
    Is the source data type used to initialize this variable an owning pointer that has been accessed?
    """

  def visit( self, node ) :
    if node.cursor.kind == CursorKind.VAR_DECL and node.cursor.type.kind == TypeKind.POINTER :
      print "visit: ", node.cursor.kind, node.cursor.type.kind, node.cursor.spelling
    if self.vardec_to_pointer_to_refcount_deriving_class( node ) :
      self.decs_to_refcount_subclass_pointers[ node.cursor.get_usr() ] = node

def owningptr_dereference( node ) :
  """
  Returns true when the given node represents a call to operator* of class owning_ptr<T>
  """
  if node.cursor.kind != CursorKind.CALL_EXPR : return False;
  if node.cursor.displayname != "operator*" : return False;
  if not node.allchild : return False
  firstchild = node.allchild[0]
  return node_is_owning_ptr( firstchild )

def node_is_owning_ptr( node ) :
  """
  Returns true if the type on an object is owning_ptr
  """
  if node.cursor.kind == CursorKind.MEMBER_REF_EXPR and node.allchild :
    return node_is_owning_ptr( node.allchild[0] )
  if node.cursor.type.kind != TypeKind.RECORD : return False
  childdec = node.cursor.type.get_declaration()
  if not childdec : return False
  if childdec.kind != CursorKind.CLASS_DECL : return False
  if childdec.spelling != "owning_ptr" : return False
  return True

def ptr_to_opheld_obj( node ) :
  """
  Recursive routine that returns true if the noode is the root of a tree
  that returns a raw pointer to the contents of an owning_ptr.
  Mutual recursion with ptr_to_opheld_obj.
  """
  if node.cursor.kind == CursorKind.CALL_EXPR :
    dspname = node.cursor.displayname
    #print dspname
    if dspname == "operator()" or dspname == "operator->" or dspname == "get" :
      #print dspname, len(node.allchild)
      #node.treeprint( lambda x : "%s %s" % (x.cursor.kind, x.cursor.type.kind) )
      if not node.allchild : return False
      return node_is_owning_ptr( node.allchild[0] )
  if node.cursor.kind == CursorKind.PAREN_EXPR :
    if not node.allchild : return False
    return ptr_to_opheld_obj( node.allchild[0] )
  if node.cursor.kind == CursorKind.UNARY_OPERATOR :
    if not node.allchild : return False
    toks = list(node.cursor.get_tokens())
    if not toks : return False
    #print "tok0", toks[0].spelling
    if toks[0].spelling == "&" :
      return ref_to_opheld_obj( node.allchild[ 0 ] )
  if node.cursor.kind == CursorKind.CXX_DYNAMIC_CAST_EXPR or node.cursor.kind == CursorKind.CXX_STATIC_CAST_EXPR :
    if node.cursor.type.kind == TypeKind.POINTER :
      if not node.allchild : return False
      return ptr_to_opheld_obj( node.allchild[-1] )
  return False

def ref_to_opheld_obj( node ) :
  """
  Recursive routine that returns true if the node is the root of a tree
  that returns a reference to the contents of an owning_ptr.
  Mutual recursion with ptr_to_opheld_obj.
  """
  if node.cursor.kind == CursorKind.CALL_EXPR :
    dspname = node.cursor.displayname
    if dspname == "operator*" :
      #print dspname
      if not node.allchild : return False
      return node_is_owning_ptr( node.allchild[0] )
  if node.cursor.kind == CursorKind.PAREN_EXPR :
    if not node.allchild : return False
    return ref_to_opheld_obj( node.allchild[0] )
  if node.cursor.kind == CursorKind.UNARY_OPERATOR :
    if not node.allchild : return False
    toks = list( node.cursor.get_tokens() )
    if not toks : return False
    #print "tok0", toks[0].spelling
    if toks[0].spelling == "*" :
      return ptr_to_opheld_obj( node.allchild[ 0 ] )
  if node.cursor.kind == CursorKind.CXX_DYNAMIC_CAST_EXPR or node.cursor.kind == CursorKind.CXX_STATIC_CAST_EXPR :
    if node.cursor.type.kind == TypeKind.RECORD :
      if not node.allchild : return False
      return ref_to_opheld_obj( node.allchild[-1] )
  return False

def stringify_token( token ) :
  return "(%s,%s,%s,%s,%s,%s,%s,%s)\n" % ( token.extent.start.file, token.extent.start.line, token.extent.start.column, token.extent.end.file, token.extent.end.line, token.extent.end.column, token.spelling, token.kind )

def binary_operator_token( node ) :
  assert( node.cursor.kind == CursorKind.BINARY_OPERATOR )
  if len(list(node.cursor.get_tokens())) == 0: return ""
  return list( node.allchild[ 0 ].cursor.get_tokens() )[ -1 ].spelling

def setup_cansampmover_ast() :
  #srcfile = "protocols/canonical_sampling/CanonicalSamplingMover.cc";
  srcfile = "protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMover.cc"
  #srcfile = "utility/pointer/ReferenceCount.hh"
  #srcfile = "stupid_test.cc"
  fn = os.path.abspath(srcfile);
  srcdir=os.path.abspath(".")
  rosetta_source_basedir="/Volumes/Macintosh HD/Users/andrew/scratch/decarboxy/GIT/newclone/rosetta/rosetta_source/"
  clangargs=["-I%s"%(srcdir),"-I%s"%(rosetta_source_basedir + "src/platform/macos"),"-I%s"%(rosetta_source_basedir + "external"), "-I%s"%(rosetta_source_basedir +"external/boost_1_46_1"), "-I%s"%(rosetta_source_basedir + "external/include")];

  print "defining source file"
  src = willclang.SourceFile(fn,clangargs);
  print "producing ast"
  ast = src.get_ast()

  #ast = willclang.get_doctest_ast( "test_apl.cc" )

  print "finding class declarations"
  cr = CodeReader(ast)
  ast.root.treemap( cr.find_class_declarations, allchild=True )
  cr.note_classes_deriving_from_reference_count()
  return src, ast, cr

class CodeQualityChecker_FindStackDeclaredObjectsReservedForTheHeap :
  def __init__( self ) :
    self.classes_that_should_only_live_on_the_heap = set( [
      "core::chemical::ResidueType"
    ])
    self.problem_string = ""

  def has_problem( self ) :
    return self.problem_string != ""

  def problems( self ) :
    return self.problem_string

  def examine_ast( self, ast ) :
    ast.root.treemap( self.find_vardecs_to_heaponly_classes, allchild=False )

  def find_vardecs_to_heaponly_classes( self, node ) :
    vardec_type = self.vardec_type_for_node( node )
    if not vardec_type : return
    scoped_class_name = vardec_type.namespace() + "::" + vardec_type.cursor.spelling;
    if scoped_class_name in self.classes_that_should_only_live_on_the_heap :
      problem = \
        "Variable declared in %s on line %d is of type %s which must only be declared on the heap.\n" % \
        ( node.cursor.location.file, node.cursor.location.line, scoped_class_name )
      problem += ">>> "
      for token in node.cursor.get_tokens() :
        problem += token.spelling
      problem += "\n"
      self.problem_string += problem

  def vardec_type_for_node( self, node ) :
    if node.cursor.kind == CursorKind.VAR_DECL and node.cursor.type.kind != TypeKind.LVALUEREFERENCE and node.cursor.type.kind != TypeKind.POINTER :
      for child in node.allchild :
        if child.cursor.kind == CursorKind.NAMESPACE_REF : continue
        if child.cursor.kind == CursorKind.TYPE_REF :
          if child.cursor.type.kind == TypeKind.RECORD :
            # direct reference to a type
            return child.definition()
          elif child.cursor.type.kind == TypeKind.TYPEDEF :
            canontype = child.cursor.type.get_canonical()
            if canontype.kind != TypeKind.RECORD : return None
            canontypedec = canontype.get_declaration()
            if canontypedec.get_usr() in node.ast.locmap :
              return node.ast.locmap[ canontypedec.get_usr() ]
    return None

#def vardec_to_residue_type( node ) :
#  def_for_node = vardec_type_for_node( node )
#  if not def_for_node : return False
#  if def_for_node.cursor.spelling == "ResidueType" and def_for_node.namespace() == "core::chemical" :
#    return True
#  return False
#
#def print_vardec_to_residue_type_nodes( node ) :
#  if vardec_to_residue_type( node ) :
#    print node.cursor.location.file, node.cursor.location.line
#    print ">>> ",
#    for token in node.cursor.get_tokens() :
#      print token.spelling,
#    print


if __name__ == "__main__" :

  src, ast, cr = setup_cansampmover_ast()

  print "found", len(cr.classes), "class declarations"
  print "found", len(cr.classes_deriving_from_refcount), "classes deriving from utility::pointer::ReferenceCount"
  for hashval, cl in cr.classes_deriving_from_refcount.items() :
    print "read class", cl.name, "in namespace", cl.ns, "in file", cl.file_declared_within
    for parentcl in cl.parents :
      print "  class", cl.name, "derives from", cr.classes[ parentcl ].name
  ast.root.treeprint( cr.visit, allchild=False )

  rfptrfinder = VarDecToRefCountSubclassPointerFinder( cr )
  ast.root.treemap( rfptrfinder.visit, allchild=False )
  print "found", len( rfptrfinder.decs_to_refcount_subclass_pointers ), "raw pointer variable declarations to classes derived from ReferenceCount"
  for hashval, node in rfptrfinder.decs_to_refcount_subclass_pointers.items() :
    print "node", node.cursor.spelling, "on line", node.cursor.location.line
