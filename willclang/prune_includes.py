import sys,os,willclang
from clang.cindex import Index,CursorKind,TypeKind

def cursor_is_includ(c):
	return c.kind==CursorKind.INCLUSION_DIRECTIVE

if __name__ == '__main__':
	if len(sys.argv)==1:
		srcdir = os.path.dirname(willclang.__file__)+"/test"
		srcfile = srcdir+"/types.cc"
		clangargs = clangargs=["-I%s"%srcdir]
		print 'parsing',srcfile,'with args',clangargs
		src = willclang.SourceFile(srcfile,clangargs)
 		ast = src.get_ast()
		print ast.root.getstr()
	else:
		for fn in sys.argv[1:]:
			fn = os.path.abspath(fn)
			src = willclang.RosettaSourceFile(fn)
			ast = src.get_ast()
			# for x in src.get_ast().root.search("NamedStubID"):
			# 	print x.cursor.get_definition()
			# sys.exit()
#			willclang.util.print_raw_ast(src.get_ast().tu.cursor)
			inc_nes,inc_why = ast.get_necessary_headers()			
			for i in inc_nes: print "PREPRUNE",i
			# for i in inc_why:
			# 	print i
			# 	for j in inc_why[i]:
			# 		print "   ",j.getstr(recursive=False),j.cursor.location
			#inc_all = list(src.ast.tu.get_includes())
			inc_provided,inc_nes = ast.transitive_include_gragh(inc_nes,inc_why)
			# for inc in inc_provided:
			# 	print inc,"provided by","the following:" if inc_provided[inc] else "NONE"
			# 	for x in inc_provided[inc]:
			# 		print "   ",x
			for ip in sorted(inc_provided.keys()):
				print ip
				for ipv in sorted(inc_provided[ip]):
					print "   ",ipv
			while True:
				torm = None
				for h in inc_nes:
					if len(inc_provided[h]):
						torm = h
						break
				if not torm: break
				inc_nes.remove(torm)
				for v in inc_provided.values(): 
					if torm in v: v.remove(torm)
			for i in inc_nes:
				if i.startswith("src/"): i = i[4:]
				print "#include <%s>"

		
	
	
