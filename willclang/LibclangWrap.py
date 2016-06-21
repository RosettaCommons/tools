import sys,os,re,math
from clang.cindex import Index,CursorKind,TypeKind
from willclang.util import *

ERRLEVEL = ("Ignored:","Note:","Warning:","Error:","Fatal:")
DECL_TYPES = (CursorKind.PARM_DECL,CursorKind.FIELD_DECL,CursorKind.VAR_DECL)

def get_doctest_file(fn):
	"""
	>>> get_doctest_file("types.cc")      #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
	'.../willclang/test/types.cc'
	"""
	fn = os.path.dirname(os.path.abspath(__file__))+"/test/"+fn
	if not os.path.exists(fn): raise IOError("can't find file: "+fn)
	return fn

def get_doctest_ast(fn):
	"""
	>>> print get_doctest_ast("types.cc")      #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
	TRANSLATION_UNIT /Users/sheffl..test/types.cc     'struct MyStru..\\tMSP msp;\\n}'
	| STRUCT_DECL MyStruct     'struct MyStruct {}'
	| TYPEDEF_DECL MSP     'typedef MyStruct MSP'
	| | TYPE_REF struct MyStruct     'MyStruct'
	| FUNCTION_DECL test_types     'void test_typ..\\tMSP msp;\\n}'
	| | COMPOUND_STMT      '{\\n\\tint i;\\n..\\tMSP msp;\\n}'
	| | | DECL_STMT      'int i;'
	| | | | VAR_DECL i     'int i'
	| | | DECL_STMT      'int *ip;'
	| | | | VAR_DECL ip     'int *ip'
	| | | DECL_STMT      'int & ir(i);'
	| | | | VAR_DECL ir     'int & ir(i'
	| | | | | DECL_REF_EXPR i     'i'
	| | | DECL_STMT      'MyStruct m;'
	| | | | VAR_DECL m     'MyStruct m'
	| | | | | TYPE_REF struct MyStruct     'MyStruct'
	| | | | | CALL_EXPR MyStruct     'm'
	| | | DECL_STMT      'MyStruct *mp;'
	| | | | VAR_DECL mp     'MyStruct *mp'
	| | | | | TYPE_REF struct MyStruct     'MyStruct'
	| | | DECL_STMT      'MyStruct &mr(m);'
	| | | | VAR_DECL mr     'MyStruct &mr(m'
	| | | | | TYPE_REF struct MyStruct     'MyStruct'
	| | | | | DECL_REF_EXPR m     'm'
	| | | DECL_STMT      'MSP msp;'
	| | | | VAR_DECL msp     'MSP msp'
	| | | | | TYPE_REF MSP     'MSP'
	| | | | | CALL_EXPR MyStruct     'msp'
	"""
	fn = get_doctest_file(fn)
	src = SourceFile(fn,clangargs=["-I"+os.path.dirname(fn)])
	return src.get_ast()

def print_raw_ast(cursor,depth=0):
	"""
	>>> ast = get_doctest_ast("test_type_simple.cc")
	>>> print_raw_ast(ast.tu.cursor)		           #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
	TRANSLATION_UNIT ...
	|  TYPEDEF_DECL __int128_t
	|  TYPEDEF_DECL __uint128_t
	|  STRUCT_DECL __va_list_tag
	|  |  FIELD_DECL gp_offset
	|  |  FIELD_DECL fp_offset
	|  |  FIELD_DECL overflow_arg_area
	|  |  FIELD_DECL reg_save_area
	|  TYPEDEF_DECL __va_list_tag
	|  |  STRUCT_DECL __va_list_tag
	|  |  |  FIELD_DECL gp_offset
	|  |  |  FIELD_DECL fp_offset
	|  |  |  FIELD_DECL overflow_arg_area
	|  |  |  FIELD_DECL reg_save_area
	|  TYPEDEF_DECL __builtin_va_list
	|  |  TYPE_REF __va_list_tag
	|  |  INTEGER_LITERAL
	|  VAR_DECL i
	|  VAR_DECL ir
	|  |  UNEXPOSED_EXPR i
	|  |  |  DECL_REF_EXPR i
	|  VAR_DECL ip
	|  |  UNARY_OPERATOR
	|  |  |  DECL_REF_EXPR i
	|  VAR_DECL ic
	|  |  UNEXPOSED_EXPR i
	|  |  |  DECL_REF_EXPR i
	"""
	print "|  "*depth+str(cursor.kind)[11:],
	print cursor.spelling if cursor.spelling else cursor.displayname
	for ch in cursor.get_children():
		print_raw_ast(ch,depth+1)

def hashcursor(c):
	"""
	get a unique string for a cursor

	>>> n = get_doctest_ast("test.cc").root.srcchild[1]
	>>> hashcursor(n.cursor)
	'c:@ir'
	"""
	return c.get_usr()


class ASTException(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)



class Include(object):
	"""docstring for IncludeTree"""
	def __init__(self,fname,namemap):
		super(Include, self).__init__()
		self.fname = fname
		self.provided_by = set()
		self.provided_by_transitive = set()
		self.namemap = namemap

	def __str__(self,depth=0,printed=None):
		s = self.fname+", provided by:"
		tmp = list(self.provided_by_transitive)
		tmp.sort()
		for c in tmp:
			s += "\n    "+c.fname
		return s

	def __cmp__(self,othr):
		return cmp(self.fname,othr.fname)



class Node(object):
	"""Wrapper of libclang ast 'cursor'"""
	def __init__(self, ast, cursor, parent=None):
		super(Node, self).__init__()
		self.ast = ast       #: AST this node belongs to
		self.cursor = cursor #: clang cursor
		self.parent = parent #: parent node
		self.depth = 1
		if parent : self.depth = parent.depth + 1
		self.code = []       #: Nchild+1 bits of code interdigitated between code owned by children
		childiter = self.cursor.get_children()
		if self.cursor.kind == CursorKind.TRANSLATION_UNIT:
			self.fname = self.ast.fname
			self.start = 0
			self.end = len(self.ast.code)
			# remove extra crap at beginning
			#assert childiter.next().spelling == "__int128_t"
			#assert childiter.next().spelling == "__uint128_t"
			#assert childiter.next().spelling == "__va_list_tag"
			#assert childiter.next().spelling == "__va_list_tag"
			#assert childiter.next().spelling == "__builtin_va_list"
		else:
			self.fname = str(self.cursor.location.file)  #: filename of sourcefile this node came from
			self.start = self.cursor.extent.start.offset #: start position in sourcefile
			self.end   = self.cursor.extent.end.offset   #: end position in sourcefile
		if self.cursor.get_definition():
		 	self.ast.locmap[hashcursor(self.cursor)] = self
		self.allchild = [Node(ast,c,self) for c in childiter] #: all children
		self.update_srcchild()

	def update_srcchild(self):
		"""
		updates list of same-srcfile children
		"""
		self.srcchild = []                                    #: children from same srcfile
		for c in self.allchild:
			if self.ast.fname != str(c.cursor.location.file): continue
			#if c.type().startswith("UNEXPOSED"): continue
			self.srcchild.append(c)
		for c in self.srcchild: c.update_srcchild()

	def includes_used(self):
		"""
		returns [sorted list of includes] {map of nodes which think they need a particular include}
		attempt to find the headers which are actually needed for the core in the sourcefile
		in the following example, many headers are included, but only something in tdef3.hh is actually used.

		WORK IN PROGRESS!

		>>> ast = get_doctest_ast("test_inc.cc")

		print all headers:

		>>> for h in sorted(ast.all_headers): print h            #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		/.../willclang/test/include/class1.hh
		/.../willclang/test/include/class2.hh
		/.../willclang/test/include/class3.hh
		/.../willclang/test/include/class4.hh
		/.../willclang/test/include/def1.hh
		/.../willclang/test/include/def2.hh
		/.../willclang/test/include/def3.hh
		/.../willclang/test/include/def4.hh
		/.../willclang/test/include/func1.hh
		/.../willclang/test/include/func2.hh
		/.../willclang/test/include/func3.hh
		/.../willclang/test/include/func4.hh
		/.../willclang/test/include/ns1.hh
		...

		but these aren't all needed, the ast it just this:

		>>> print ast.root.getstr()                               #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		TRANSLATION_UNIT ...
		| VAR_DECL i     'tdef3 i'
		| | TYPE_REF tdef3     'tdef3'

		print only headers "referred" to in the actual source:

		>>> inc,inc_why = ast.root.includes_used()
		>>> for h in sorted(inc): print h                          #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		/.../willclang/test/include/tdef3.hh

		The node(s) which think they need a particular header are reported in the map 'inc_why' (2nd return):

		>>> for h in inc:                                          #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		...    print h
		...    for n in inc_why[h]:
		...        print "    needed by",n.getstr()
		/.../willclang/test/include/tdef3.hh
		    needed by TYPE_REF tdef3     'tdef3'
		"""
		inc = set()
		incwhy = dict()
		if self.cursor:
			td = self.cursor
			if td.type.get_canonical().get_declaration().location.file:
				h = str(td.type.get_canonical().get_declaration().location.file)
				if h != self.ast.fname:
					inc.add(h)
					incwhy.setdefault(h,[]).append(self)
			if td.kind == CursorKind.TYPE_REF:
				while td.kind == CursorKind.TYPE_REF:
					td = td.get_definition()
					if td is None: break
					chld = list(td.get_children())
					if td.location.file:
						fn = str(td.location.file)
						if self.parent.is_decl():
							if self.parent.is_ptr_or_ref():	fn = self.ast.switch_to_fwd(fn)
							else:                           fn = self.ast.switch_to_hh(fn)
							if fn != self.ast.fname:
								inc.add(fn)
								incwhy.setdefault(fn,[]).append(self)
					if len(chld)==1: td = chld[0]
		for c in self.srcchild:
		 	cinc,cincwhy = c.includes_used()
			inc |= cinc
			for k in incwhy: cincwhy.setdefault(k,[]).extend(incwhy[k])
			incwhy = cincwhy
		return inc,incwhy

	def search(self,name):
		rslt = []
		if self.cursor and (str(self.cursor.displayname).endswith(name) or str(self.cursor.spelling).endswith(name)):
			rslt.append(self)
		for c in self.srcchild: rslt.extend(c.search(name))
		return rslt

	def splitcode(self):
		if self.fname == self.ast.fname:
			if self.srcchild:
				self.code.append(self.ast.code[self.start:self.srcchild[0].start])
				for i in range(1,len(self.srcchild)):
					#assert self.srcchild[i-1].end <= self.srcchild[i].start
					self.code.append(self.ast.code[self.srcchild[i-1].end:self.srcchild[i].start])
				self.code.append(self.ast.code[self.srcchild[-1].end:self.end])
			else:
				self.code.append(self.ast.code[self.start:self.end])
			for c in self.srcchild: c.splitcode()

	def isdefinition(self):
		"""
		>>> def visit(node): print "IS_DEF:",str(node.isdefinition()).ljust(5),node.getstr(recursive=False)
		>>> ast = get_doctest_ast("types.cc")
		>>> ast.root.treemap(visit) #doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
		IS_DEF: False TRANSLATION_UNIT /Users/sheffl..test/types.cc     'struct MyStru..\\tMSP msp;\\n}'
		IS_DEF: True  STRUCT_DECL MyStruct     'struct MyStruct {}'
		IS_DEF: True  TYPEDEF_DECL MSP     'typedef MyStruct MSP'
		IS_DEF: False TYPE_REF struct MyStruct     'MyStruct'
		IS_DEF: True  FUNCTION_DECL test_types     'void test_typ..\\tMSP msp;\\n}'
		IS_DEF: False COMPOUND_STMT      '{\\n\\tint i;\\n..\\tMSP msp;\\n}'
		IS_DEF: False DECL_STMT      'int i;'
		IS_DEF: True  VAR_DECL i     'int i'
		IS_DEF: False DECL_STMT      'int *ip;'
		IS_DEF: True  VAR_DECL ip     'int *ip'
		...
		"""
		return self.cursor.is_definition()

	def definition(self):
		"""
		>>> ast = get_doctest_ast("test_ref.cc")
		>>> print ast.root.getstr(recurseall=True)                     #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		TRANSLATION_UNIT /.../test_ref.cc     '#include <tes..sfunc(ms);\\n}'
		| STRUCT_DECL MyStruct
		| STRUCT_DECL MyStruct
		| | FIELD_DECL i
		| FUNCTION_DECL msfunc
		| | PARM_DECL ms
		| | | TYPE_REF struct MyStruct
		| | COMPOUND_STMT
		| | | RETURN_STMT
		| | | | UNEXPOSED_EXPR i
		| | | | | MEMBER_REF_EXPR i
		| | | | | | DECL_REF_EXPR ms
		| FUNCTION_DECL func1     'void func1(My..sfunc(ms);\\n}'
		| | PARM_DECL ms     'MyStruct & ms'
		| | | TYPE_REF struct MyStruct     'MyStruct'
		| | COMPOUND_STMT      '{\\n\\tint i = msfunc(ms);\\n}'
		| | | DECL_STMT      'int i = msfunc(ms);'
		| | | | VAR_DECL i     'int i = msfunc(ms)'
		| | | | | CALL_EXPR msfunc     'msfunc(ms)'
		| | | | | | UNEXPOSED_EXPR msfunc     'msfunc'
		| | | | | | | DECL_REF_EXPR msfunc     'msfunc'
		| | | | | | DECL_REF_EXPR ms     'ms'
		>>> ast.root.treeprint(lambda n: n.definition(),allchild=True) #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		None
		| STRUCT_DECL MyStruct
		| STRUCT_DECL MyStruct
		| | FIELD_DECL i
		| FUNCTION_DECL msfunc
		| | PARM_DECL ms
		| | | STRUCT_DECL MyStruct
		| | None
		| | | None
		| | | | FIELD_DECL i
		| | | | | FIELD_DECL i
		| | | | | | PARM_DECL ms
		| FUNCTION_DECL func1     'void func1(My..sfunc(ms);\\n}'
		| | PARM_DECL ms     'MyStruct & ms'
		| | | STRUCT_DECL MyStruct
		| | None
		| | | None
		| | | | VAR_DECL i     'int i = msfunc(ms)'
		| | | | | FUNCTION_DECL msfunc
		| | | | | | FUNCTION_DECL msfunc
		| | | | | | | FUNCTION_DECL msfunc
		| | | | | | PARM_DECL ms     'MyStruct & ms'
		"""
		d = self.cursor.get_definition()
		if not d: return None
		try: return self.ast.locmap[hashcursor(d)]
		except KeyError: return None

	def type(self):
		return str(self.cursor.kind)[11:]

	def children_of_type(self,t,recurseall=False):
		"""
		>>> ast = get_doctest_ast("types.cc")
		>>> for decl in ast.root.children_of_type("DECL_STMT"):
		...    print decl
		DECL_STMT      'int i;'
		DECL_STMT      'int *ip;'
		DECL_STMT      'int & ir(i);'
		DECL_STMT      'MyStruct m;'
		DECL_STMT      'MyStruct *mp;'
		DECL_STMT      'MyStruct &mr(m);'
		DECL_STMT      'MSP msp;'
		"""
		n = [self] if self.type() == t else list()
		for c in (self.allchild if recurseall else self.srcchild):
			n.extend(c.children_of_type(t))
		return n

	def namespace(self):
		"""
		get the full namespace of any-old-thing. you *can't* do this with grep!

		.. todo::
			find another way to get namespace of arbitrary things, because this the parent-reference approach
			is really the only reason each node needs a wrapper. this is expensive

		>>> ast = get_doctest_ast("test_ns.cc")
		>>> visit = lambda x: x.codeline() + (" IN " if x.namespace() else "") + x.namespace()
		>>> ast.root.treeprint(visit)     #doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
		#...
		| | | ns1::val1                     IN ns1                     ...
		| | | ns1::ns2::val2                IN ns1::ns2                ...
		| | | ns1::ns2::ns3::val3           IN ns1::ns2::ns3           ...
		| | | ns1::ns2::ns3::ns4::val4      IN ns1::ns2::ns3::ns4      ...
		| | | ns1::ns2::ns3::ns4::ns5::val5 IN ns1::ns2::ns3::ns4::ns5 ...
		| | | using namespace ns1; ...
		| | | using namespace ns2; ...
		| | | using namespace ns3; ...
		| | | val1           IN ns1
		| | | val2           IN ns1::ns2
		| | | val3           IN ns1::ns2::ns3
		| | | ns4::val4      IN ns1::ns2::ns3::ns4      ...
		| | | ns4::ns5::val5 IN ns1::ns2::ns3::ns4::ns5 ...
		| | | using namespace ns4; ...
		| | | using namespace ns5; ...
		| | | val1 IN ns1
		| | | val2 IN ns1::ns2
		| | | val3 IN ns1::ns2::ns3
		| | | val4 IN ns1::ns2::ns3::ns4
		| | | val5 IN ns1::ns2::ns3::ns4::ns5
		"""
		if not self.cursor: return None
		ns = ""
		p = self.definition()
		if not p: return ""
		while p and p.cursor:
			if p.cursor.kind == CursorKind.NAMESPACE:
				ns = p.cursor.spelling + "::" + ns
			p = p.parent
		return ns.rstrip(":")

	def codeline(self):
		"""
		>>> n = get_doctest_ast("types.cc").root.search("test_types")[0]
		>>> print n.codeline()
		void test_types() {i ...| ct &mr(m);MSP msp;}
		"""
		assert self.fname == self.ast.fname
		code = self.ast.code[self.start:self.end]
		l = (code if len(code) < 50 else code[:22]+" ...| "+code[-22:])
		return l.replace("\n","").replace("\t","")

	def treemap(self,func,allchild=False):
		"""
		>>> def visit(node): print node.cursor.location
		>>> ast = get_doctest_ast("types.cc")
		>>> ast.root.treemap(visit) #doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
		<SourceLocation file None, line 0, column 0>
		<SourceLocation file '.../willclang/test/types.cc', line 1, column 8>
		<SourceLocation file '.../willclang/test/types.cc', line 3, column 18>
		<SourceLocation file '.../willclang/test/types.cc', line 3, column 9>
		<SourceLocation file '.../willclang/test/types.cc', line 5, column 6>
		<SourceLocation file '.../willclang/test/types.cc', line 5, column 19>
		<SourceLocation file '.../willclang/test/types.cc', line 6, column 2>
		<SourceLocation file '.../willclang/test/types.cc', line 6, column 6>
		<SourceLocation file '.../willclang/test/types.cc', line 7, column 2>
		...
		"""
		func(self)
		children = self.allchild if allchild else self.srcchild
		for c in children:
			c.treemap(func,allchild)

	def treeprint(self,func,allchild=False,depth=0):
		"""
		>>> ast = get_doctest_ast("test_ref.cc")
		>>> ast.root.treeprint(lambda n: "FOO "+str(n))
		FOO TRANSLATION_UNIT /Users/sheffl..t/test_ref.cc     '#include <tes..sfunc(ms);\\n}'
		| FOO FUNCTION_DECL func1     'void func1(My..sfunc(ms);\\n}'
		| | FOO PARM_DECL ms     'MyStruct & ms'
		| | | FOO TYPE_REF struct MyStruct     'MyStruct'
		| | FOO COMPOUND_STMT      '{\\n\\tint i = msfunc(ms);\\n}'
		| | | FOO DECL_STMT      'int i = msfunc(ms);'
		| | | | FOO VAR_DECL i     'int i = msfunc(ms)'
		| | | | | FOO CALL_EXPR msfunc     'msfunc(ms)'
		| | | | | | FOO UNEXPOSED_EXPR msfunc     'msfunc'
		| | | | | | | FOO DECL_REF_EXPR msfunc     'msfunc'
		| | | | | | FOO DECL_REF_EXPR ms     'ms'
		"""
		print "| "*depth+str(func(self))
		children = self.allchild if allchild else self.srcchild
		for c in children:
			c.treeprint(func,allchild,depth+1)

	def is_decl(self):
		"""
		>>> def visit(node): print "IS_DECL:",node.is_decl(),node.getstr(recursive=False)
		>>> ast = get_doctest_ast("types.cc")
		>>> ast.root.treemap(visit) #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		IS_DECL: False TRANSLATION_UNIT ...test/types.cc     ...
		IS_DECL: False STRUCT_DECL MyStruct       'struct MyStruct {}'
		...
		IS_DECL: False DECL_STMT                  'int i;'
		IS_DECL: True  VAR_DECL i                 'int i'
		IS_DECL: False DECL_STMT                  'int *ip;'
		IS_DECL: True  VAR_DECL ip                'int *ip'
		IS_DECL: False DECL_STMT                  'int & ir(i);'
		IS_DECL: True  VAR_DECL ir                'int & ir(i'
		...
		"""
		return self.cursor.kind in DECL_TYPES

	def is_ptr_or_ref(self):
		"""
		returns True if this node is a reference or pointer type;
		returns False if this node is a non-reference and non-pointer type;
		should return None if not a declaration

		WORK IN PROGRESS

		>>> def visit(node): print "IS_PtrRef:",node.is_ptr_or_ref(),node.getstr(recursive=False)
		>>> ast = get_doctest_ast("types.cc")
		>>> ast.root.treemap(visit)                            #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		IS_PtrRef: None  TRANSLATION_UNIT ...test/types.cc ...
		IS_PtrRef: None  STRUCT_DECL MyStruct         'struct MyStruct {}'
		...
		IS_PtrRef: False VAR_DECL i                   'int i'
		IS_PtrRef: None  DECL_STMT                    'int *ip;'
		IS_PtrRef: True  VAR_DECL ip                  'int *ip'
		IS_PtrRef: None  DECL_STMT                    'int & ir(i);'
		IS_PtrRef: True  VAR_DECL ir                  'int & ir(i'
		IS_PtrRef: None  DECL_REF_EXPR i              'i'
		IS_PtrRef: None  DECL_STMT                    'MyStruct m;'
		IS_PtrRef: False VAR_DECL m                   'MyStruct m'
		...
		IS_PtrRef: True  VAR_DECL mp                  'MyStruct *mp'
		IS_PtrRef: None  TYPE_REF struct MyStruct     'MyStruct'
		IS_PtrRef: None  DECL_STMT                    'MyStruct &mr(m);'
		IS_PtrRef: True  VAR_DECL mr                  'MyStruct &mr(m'
		...
		IS_PtrRef: False VAR_DECL msp                 'MSP msp'
		...
		"""
		if not self.is_decl(): return None
		for c in self.cursor.get_children():
			if c.kind==CursorKind.TEMPLATE_REF and c.displayname=="owning_ptr":	return True
			if c.kind==CursorKind.TYPE_REF:
				if c.displayname.endswith("OP"): return True
				if c.displayname.endswith("OPs"): return True
				if c.displayname.endswith("CAP"): return True
		return self.cursor.type.get_canonical().get_pointee().kind!=TypeKind.INVALID

	def getstr(self,recursive=True,recurseall=False,depth=0):
		"""
		make a printable string from a node

		say raw source is this (test2.hh contains ns1::i):

		>>> ast = get_doctest_ast("test2.cc")
		>>> print ast.code                                          #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		#include <test2.hh>
		using namespace ns1;
		int j = ns1::i;

		raw ast looks like this:

		>>> print_raw_ast(ast.tu.cursor)                             #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		TRANSLATION_UNIT ...
		|  NAMESPACE ns1
		|  |  VAR_DECL i
		|  USING_DIRECTIVE
		|  |  NAMESPACE_REF ns1
		|  VAR_DECL j
		|  |  UNEXPOSED_EXPR i
		|  |  |  DECL_REF_EXPR i
		|  |  |  |  NAMESPACE_REF ns1

		by default prints only nodes from parsed src file:

		>>> print ast.root.getstr(recursive=True,recurseall=False)   #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		TRANSLATION_UNIT ...
		| USING_DIRECTIVE      'using namespace ns1'
		| | NAMESPACE_REF ns1     'ns1'
		| VAR_DECL j     'int j = ns1::i'
		| | UNEXPOSED_EXPR i     'ns1::i'
		| | | DECL_REF_EXPR i     'ns1::i'
		| | | | NAMESPACE_REF ns1     'ns1'

		with recursive=False prints only this node

		>>> print ast.root.getstr(recursive=False)                   #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		TRANSLATION_UNIT ...

		with recurseall, you get lines from included headers (NAMESPACE,VAR_DECL):

		>>> print ast.root.getstr(recursive=True,recurseall=True)    #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		TRANSLATION_UNIT ...
		| NAMESPACE ns1
		| | VAR_DECL i
		| USING_DIRECTIVE      'using namespace ns1'
		| | NAMESPACE_REF ns1     'ns1'
		| VAR_DECL j     'int j = ns1::i'
		| | UNEXPOSED_EXPR i     'ns1::i'
		| | | DECL_REF_EXPR i     'ns1::i'
		| | | | NAMESPACE_REF ns1     'ns1'
		"""
		s = "| "*depth
		subcode = ""#SOURCE UNREAD FOR: "+self.fname
		if self.fname==self.ast.fname:
			subcode = self.ast.code[self.start:self.end]
			subcode = subcode.replace("\n","\\n").replace("\t","\\t")
		if len(subcode) > 30: subcode = subcode[:13]+".."+subcode[-13:]
		name = self.cursor.spelling if self.cursor.spelling else self.cursor.displayname
		if len(name) > 30: name = name[:13]+".."+name[-13:]
		s += self.type()+" "+name
		if subcode: s += "     '"+subcode+"'"
		#s += '"'+self.codeline()+'" '
		#s += str(self.code)+" "
		#s += os.path.basename(str(self.cursor.location.file))+":"+str(self.cursor.extent.start.line)
		if recursive:
			chld = self.allchild if recurseall else self.srcchild
			for c in chld:
				s += "\n"
				s += c.getstr(recursive,recurseall,depth+1)
		return s

	def gencode(self,SP=""):
		s = ""
		if self.allchild:
			if self.code: s += self.code[0]
			for i in range(len(self.srcchild)):
				s += self.srcchild[i].gencode(SP)
				if self.code: s += self.code[i+1]
		else:
			if self.code: s += self.code[0]
		return SP+s+SP

	def __str__(self):
		return self.getstr(0,False)

	def check(self,recursive=True):
		"""
		>>> node = get_doctest_ast("test.cc").root
		>>> node.check()
		>>> test = node
		>>> test.allchild[1].allchild[0].fname = "NONSENSE"
		>>> test.check(recursive=False)
		>>> test.check()
		Traceback (most recent call last):
		    assert self.fname == str(self.cursor.location.file)
		AssertionError
		>>> test = node
		>>> test.allchild[1].allchild = []
		>>> test.check()
		Traceback (most recent call last):
		    for c in self.srcchild: assert c in self.allchild
		AssertionError
		>>> test.fname = str(node.cursor.location.file)
		>>> test.check()
		Traceback (most recent call last):
		    assert self.fname == self.ast.fname
		AssertionError
		"""
		if self.type()=="TRANSLATION_UNIT":
			assert self.fname == self.ast.fname
		else:
			assert self.fname == str(self.cursor.location.file)
		assert len(self.code) == len(self.srcchild)+1 or (len(self.code) is 0 and len(self.srcchild) is 0)
		for c in self.srcchild: assert c in self.allchild
		if self.parent:
			if not self.parent.fname == self.fname: assert self in self.parent.allchild
			else: assert self in self.parent.srcchild  and self in self.parent.allchild
		if not recursive: return
		for c in self.allchild: c.check()

	def checksrcoverlap(self):
		fail = False
		for c in self.srcchild:
			if self.start > c.start or self.end < c.end:
				fail = True
		for i in range(1,len(self.srcchild)):
			if self.srcchild[i-1].end > self.srcchild[i].start:
				fail = True
		if fail:
			# print "I",self.code[0]
			# for i in range(len(self.srcchild)):
			# 	print "C",self.srcchild[i].gencode()
			# 	print "I",self.code[i+1]
			# print
			if self.type() in CUSTOM_CGEN_TYPES:
				return
			CUSTOM_CGEN_TYPES.append( self.type() )
			print self.type()
			# print str(self.start).rjust(5)+"-"+str(self.end).ljust(5),self
			# for c in self.srcchild:
			# 	print str(c.start).rjust(5)+"-"+str(c.end).ljust(5),"   ",c
			# print
		for c in self.srcchild: c.checksrcoverlap()




class AST(object):
	"""
	Wrapper of libclang ast

	>>> ast = get_doctest_ast("test.cc")
	"""
	def __init__(self,src,clangargs=None,clangwd=None):
		"""
		recursively create Node wrapper around chang ast cursors
		"""
		super(AST, self).__init__()
		if not clangargs: clangargs = list()
		self.src = src
		self.fname = src.fname
		pform = None
		self.locmap = {}
		with open(src.fname) as o: self.code = o.read()
		clang = Index.create()
		with open(self.src.fname) as o:
			self.code = o.read()
		wd = os.getcwd()
		if clangwd:
			os.chdir(clangwd)
		ClangOpt_None = 0x0
		ClangOpt_DetailedPreprocessingRecord = 0x01
		ClangOpt_Incomplete = 0x02
		ClangOpt_PrecompiledPreamble = 0x04
		ClangOpt_CacheCompletionResults = 0x08
		ClangOpt_CXXPrecompiledPreamble = 0x10
		ClangOpt_CXXChainedPCH = 0x20
		ClangOpt_NestedMacroExpansions = 0x40
		#print "libclang parsing"
		self.tu = clang.parse(src.fname,args=clangargs,options=ClangOpt_None)
		os.chdir(wd)
		#print "parse done"
		self.check_tu()
		#print "building wrapper tree"
		self.root = Node(self,self.tu.cursor,None)
		# self.root.allchild = []
		# self.root.fname = self.fname
		# for x in self.tu.cursor.get_children():
		# 	# fn = str(x.location.file)
		# 	# if fn == self.fname or fn[0] != "/":
		# 			self.root.allchild.append( Node(self,x,None) )
		# print self.root.allchild[5]
		#print "update srcchild"
		self.root.update_srcchild()
		#print "dividing code"
		self.root.splitcode()
		#print "consistency check"
		#self.root.check()
		self.all_headers = self.get_all_headers()

	def get_all_headers(self):
		"""
		>>> ast = get_doctest_ast("test_fwd.cc")
		>>> for i in sorted(ast.get_all_headers()): print i     #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		/.../willclang/test/include/class1.fwd.hh
		/.../willclang/test/include/class1.hh
		/.../willclang/test/include/class2.hh
		/.../willclang/test/include/class3.hh
		"""
		return set(str(x.include) for x in self.tu.get_includes())

	def check_tu(self):
		"""
		check for parse errors

		>>> ast = get_doctest_ast("bad.cc")               #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		Traceback (most recent call last):
		ASTException: "Parse failed from Errors:...Error: no member named 'not_func' in namespace 'foo::bar' AT .../test/bad.cc line 13 col 38"
		"""
		fail = False
		msg = ""
		for err in list(self.tu.diagnostics):
			msg += "\n"+str(ERRLEVEL[err.severity])+" "+str(err.spelling)+" AT "+str(err.location.file)+" line "+str(err.location.line)+" col "+str(err.location.column)
			sys.stdout.flush()
			if err.severity >= err.Error: fail = True
		if fail:
			raise ASTException("Parse failed from Errors: "+msg)

	def check(self):
		"""
		consistency check for ast

		>>> ast = get_doctest_ast("bad.cc")               #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		Traceback (most recent call last):
		ASTException: "Parse failed from Errors:...Error: no member named 'not_func' in namespace 'foo::bar' AT .../test/bad.cc line 13 col 38"
		>>> ast = get_doctest_ast("test.cc")
		>>> ast.root.allchild[2].fname = "Nonsense"
		>>> ast.check()
		Traceback (most recent call last):
		    assert self.fname == str(self.cursor.location.file)
		AssertionError
		"""
		self.check_tu()
		self.root.check()

	def get_necessary_headers(self):
		"""
		wrapper for self.root.includes_used()
		"""
		inc,incwhy = self.root.includes_used()
		return inc,incwhy

	def transitive_include_gragh(self,inc_subset,inc_why=None):
		"""
		return {map of headers H to sets of all other headers which provide them}
		find all headers which (possibly transitively) include a header H

		.. todo::
			replace the lame Include class with a map or something, no reason for it to exist right now

		>>> ast = get_doctest_ast("test_inc.cc")
		>>> necessary_headers = [str(x.include) for x in ast.tu.get_includes()]
		>>> inodes,inc_subset = ast.transitive_include_gragh(necessary_headers)
		>>> for k in sorted(inodes.keys()):          #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		...    print os.path.basename(k)
		...    for i in sorted(inodes[k]): print "    provided by:",os.path.basename(i)
		class1.hh
		class2.hh
		class3.hh
		...
		tdefchain1.hh
		tdefchain2.hh
		    provided by: tdefchain1.hh
		tdefchain3.hh
		    provided by: tdefchain1.hh
		    provided by: tdefchain2.hh
		tdefchain4.hh
		    provided by: tdefchain1.hh
		    provided by: tdefchain2.hh
		    provided by: tdefchain3.hh
		...
		>>> ast.transitive_include_gragh(["imaginary_header.hh"])    #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		Traceback (most recent call last):
		ASTException:
		'Bad inc_subset! Include: imaginary_header.hh not referenced in clang parse of: /.../willclang/test/test_inc.cc!'
		"""
		inc_all = list(self.tu.get_includes())
		abs2relpath = dict()
		root = None
		# inc_subset from ast parse may be leaf in include graph
		# so abs2relpath mapping not in all includes
		for cur in inc_all:
			dn = str(cur.include)
			for i in inc_subset:
				if i.endswith(dn):
					abs2relpath[i] = dn
		# make name mapping
		for cur in inc_all:
			fn = str(cur.location.file)
			incfn = fn
			for cur2 in inc_all:
				dn = str(cur2.include)
		 		if fn.endswith(dn):
					if len(dn) < len(incfn) and incfn[-len(dn)-1]=="/":
						incfn = dn
			abs2relpath[fn] = incfn
	 	# for k in sorted(abs2relpath.keys()):
	 	# 	 	print k.ljust(100),abs2relpath[k]
		try:
			inc_subset = [abs2relpath[x] for x in inc_subset]
		except KeyError as e:
			msg = "Bad inc_subset! Include: "+e.message+" not referenced in clang parse of: "+self.fname+"!"
			if inc_why and e.message in inc_why:
				msg += "\ninclude deemed necessary because:"
				for n in inc_why[e.message]:
					msg += "\n    "+n.fname+":"+str(n.cursor.location.line)+" "+n.getstr(recursive=False)
			raise ASTException(msg)
		t = dict()
		inodes = dict()
		# incrs = set()
		for cur in inc_all:
			incer = abs2relpath[str(cur.location.file)]
			incee = str(cur.include)
			if incer not in inodes: inodes[incer] = Include(incer,abs2relpath)
			if incee not in inodes: inodes[incee] = Include(incee,abs2relpath)
			#print "IOAERSN",incer,incee
			inodes[incee].provided_by.add(inodes[incer])
		for fn in inodes:
			inodes[fn].provided_by_transitive |= inodes[fn].provided_by
		for i in range(int(math.log(len(inodes))+1)): # logN times -- like matrix mult to get X^N
			for fn in inodes:
				tmp = set()
				for inode2 in inodes[fn].provided_by_transitive:
				 	tmp |= inode2.provided_by_transitive
				inodes[fn].provided_by_transitive |= tmp
		inodes_nes = dict()
		for k,v in inodes.items():
			if k in inc_subset:
				inodes_nes[k] = set(x.fname for x in v.provided_by_transitive if x.fname in inc_subset)
		return inodes_nes,inc_subset

	def switch_to_fwd(self,fn):
		"""
		>>> ast = get_doctest_ast("test_fwd.cc")
		>>> for h in sorted(ast.all_headers): print h                       #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		/.../willclang/test/include/class1.fwd.hh
		/.../willclang/test/include/class1.hh
		/.../willclang/test/include/class2.hh
		/.../willclang/test/include/class3.hh
		>>> for h in sorted(ast.all_headers): print ast.switch_to_fwd(h)	  #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		/.../willclang/test/include/class1.fwd.hh
		/.../willclang/test/include/class1.fwd.hh
		/.../willclang/test/include/class2.hh
		/.../willclang/test/include/class3.hh
		"""
		if not fn.endswith(".fwd.hh") and fn.endswith(".hh"):
			tmp = fn[:-3]+".fwd.hh"
			if tmp in self.all_headers: return tmp
		return fn

	def switch_to_hh(self,fn):
		"""
		>>> ast = get_doctest_ast("test_fwd.cc")
		>>> for h in sorted(ast.all_headers): print h		                    #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		/.../willclang/test/include/class1.fwd.hh
		/.../willclang/test/include/class1.hh
		/.../willclang/test/include/class2.hh
		/.../willclang/test/include/class3.hh
		>>> for h in sorted(ast.all_headers): print ast.switch_to_hh(h)	   #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
		/.../willclang/test/include/class1.hh
		/.../willclang/test/include/class1.hh
		/.../willclang/test/include/class2.hh
		/.../willclang/test/include/class3.hh
		"""
		if fn.endswith(".fwd.hh"):
			tmp = fn[:-7]+".hh"
			if tmp in self.all_headers: return tmp
		return fn

	def __str__(self):
		return self.root.getstr()


class SourceFile(object):
	"""
	SourceFile holds an individual C++ sourcefile
	this should eventually include non-ast locic to do things fast, if possible
	resort to ast only if necessary
	"""
	def __init__(self,fname,clangargs=None,clangwd=None):
		super(SourceFile, self).__init__()
		self.fname = os.path.abspath(fname)
		# self.code = None
		self.ast = None
		if not clangargs: clangargs = list()
		self.clangargs = clangargs
		self.clangwd = clangwd

	# def get_code(self):
		# if not self.code:
		# 	with open(fname) as o:
		# 		self.code = o.read()

	def get_ast(self):
		if not self.ast:
			self.ast = AST(self,self.clangargs,self.clangwd)
		return self.ast



class RosettaSourceFile(SourceFile):
	"""Same as SourceFile, except sepcifies rosetta includes for clang"""
	def __init__(self, fname):
		clangargs = ["-Isrc","-Iexternal","-Iexternal/include","-Iexternal/dbio","-Isrc/platform/"+rosetta_platform()]
		clangwd="/Users/sheffler/svn/rosetta/rosetta_source/"
		if not os.path.exists(clangwd):
			clangwd = None
		super(RosettaSourceFile,self).__init__(fname,clangargs,clangwd)



def setup_extra_doctests():
	TEST = dict()
	TEST["test_inc_all"] = """
	>>> ast = get_doctest_ast("test_inc.cc")
	>>> print ast.root.getstr(recurseall=True)         #doctest: +ELLIPSIS, +REPORT_NDIFF, +NORMALIZE_WHITESPACE
	TRANSLATION_UNIT ...
	| STRUCT_DECL class1
	| | FIELD_DECL i
	| STRUCT_DECL class2
	| | FIELD_DECL i
	| STRUCT_DECL class3
	| | FIELD_DECL i
	| STRUCT_DECL class4
	| | FIELD_DECL i
	| VAR_DECL def1
	| | INTEGER_LITERAL
	| VAR_DECL def2
	| VAR_DECL def3
	| | INTEGER_LITERAL
	| VAR_DECL def4
	| | INTEGER_LITERAL
	| FUNCTION_DECL func1
	| | PARM_DECL t
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | UNEXPOSED_EXPR t
	| | | | | DECL_REF_EXPR t
	| FUNCTION_DECL func2
	| | PARM_DECL t
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | UNEXPOSED_EXPR t
	| | | | | DECL_REF_EXPR t
	| FUNCTION_DECL func3
	| | PARM_DECL t
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | UNEXPOSED_EXPR t
	| | | | | DECL_REF_EXPR t
	| FUNCTION_DECL func4
	| | PARM_DECL t
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | UNEXPOSED_EXPR t
	| | | | | DECL_REF_EXPR t
	| NAMESPACE ns1
	| | VAR_DECL i
	| NAMESPACE ns2
	| | VAR_DECL i
	| NAMESPACE ns3
	| | VAR_DECL i
	| NAMESPACE ns4
	| | VAR_DECL i
	| CLASS_TEMPLATE tclass1
	| | TEMPLATE_TYPE_PARAMETER T
	| | FIELD_DECL i
	| | | TYPE_REF T
	| CLASS_TEMPLATE tclass2
	| | TEMPLATE_TYPE_PARAMETER T
	| | FIELD_DECL i
	| | | TYPE_REF T
	| CLASS_TEMPLATE tclass3
	| | TEMPLATE_TYPE_PARAMETER T
	| | FIELD_DECL i
	| | | TYPE_REF T
	| CLASS_TEMPLATE tclass4
	| | TEMPLATE_TYPE_PARAMETER T
	| | FIELD_DECL i
	| | | TYPE_REF T
	| TYPEDEF_DECL tdef1
	| TYPEDEF_DECL tdef2
	| TYPEDEF_DECL tdef3
	| TYPEDEF_DECL tdef4
	| TYPEDEF_DECL tdefchain4
	| TYPEDEF_DECL tdefnochain4
	| TYPEDEF_DECL tdefchain3
	| | TYPE_REF tdefchain4
	| TYPEDEF_DECL tdefnochain3
	| TYPEDEF_DECL tdefchain2
	| | TYPE_REF tdefchain3
	| TYPEDEF_DECL tdefnochain2
	| TYPEDEF_DECL tdefchain1
	| | TYPE_REF tdefchain2
	| TYPEDEF_DECL tdefnochain1
	| TYPEDEF_DECL tdefchain4
	| TYPEDEF_DECL tdefnochain4
	| TYPEDEF_DECL tdefchain3
	| | TYPE_REF tdefchain4
	| TYPEDEF_DECL tdefnochain3
	| TYPEDEF_DECL tdefchain2
	| | TYPE_REF tdefchain3
	| TYPEDEF_DECL tdefnochain2
	| TYPEDEF_DECL tdefchain4
	| TYPEDEF_DECL tdefnochain4
	| TYPEDEF_DECL tdefchain3
	| | TYPE_REF tdefchain4
	| TYPEDEF_DECL tdefnochain3
	| TYPEDEF_DECL tdefchain4
	| TYPEDEF_DECL tdefnochain4
	| FUNCTION_TEMPLATE tfunc1
	| | TEMPLATE_TYPE_PARAMETER T
	| | PARM_DECL t
	| | | TYPE_REF T
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | DECL_REF_EXPR t
	| FUNCTION_TEMPLATE tfunc2
	| | TEMPLATE_TYPE_PARAMETER T
	| | PARM_DECL t
	| | | TYPE_REF T
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | DECL_REF_EXPR t
	| FUNCTION_TEMPLATE tfunc3
	| | TEMPLATE_TYPE_PARAMETER T
	| | PARM_DECL t
	| | | TYPE_REF T
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | DECL_REF_EXPR t
	| FUNCTION_TEMPLATE tfunc4
	| | TEMPLATE_TYPE_PARAMETER T
	| | PARM_DECL t
	| | | TYPE_REF T
	| | COMPOUND_STMT
	| | | RETURN_STMT
	| | | | DECL_REF_EXPR t
	"""
	return TEST

__TEST__ = setup_extra_doctests()

if __name__ == '__main__':
	# ast = get_doctest_ast("template.cc")
	# print_raw_ast(ast.root.cursor)
	# inc_nes = ast.src.get_necessary_headers()
	# for i in inc_nes: print "PREPRUNE",i
#	print ast.root.getstr(recursive=True)
#	print ast.root.srcchild[-1].is_ptr_or_ref()
	import doctest
	tr = doctest.testmod()
	print "tests passed:",tr.attempted-tr.failed
	print "tests failed:",tr.failed
