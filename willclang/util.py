import sys,os,re,platform

ERRLEVEL = ("Ignored:","Note:","Warning:","Error:","Fatal:")

def rosetta_platform():
	if   platform.platform().startswith("Darwin"): pform = "macos"
	elif platform.platform().startswith("Linux" ): pform = "linux"
	else: raise ASTException("AST doesn't recognize platform: %s"%platform.platform())
	return pform


def rospath(fname,checkfs=True):
	"""
	get path to rosetta_source
	
	>>> fname = "/SOME_DIRECTORY/SOME_OTHER_DIR/rosetta_source/src/core/scoring"
	>>> rospath(fname,checkfs=False)
	'/SOME_DIRECTORY/SOME_OTHER_DIR/rosetta_source'
	>>> fname == rospath(fname,checkfs=False)+"/"+rosbase(fname,checkfs=False)
	True
	"""
	if checkfs: assert os.path.exists(fname)
	if checkfs: fname = os.path.abspath(fname)
	fname = fname.rstrip("/")
	mark = "rosetta_source/src"
	assert fname.find(mark) > 0
	r = fname[:fname.find(mark)+len(mark)-4]	
	return r


def rosbase(fname,checkfs=True):
	"""
	get the rosetta-relative path
	
	>>> fname = "/SOME_DIRECTORY/SOME_OTHER_DIR/rosetta_source/src/core/scoring/methods"
	>>> rosbase(fname,checkfs=False)
	'src/core/scoring/methods'
	>>> fname == rospath(fname,checkfs=False)+"/"+rosbase(fname,checkfs=False)
	True
	"""
	if checkfs: assert os.path.exists(fname)
	if checkfs: fname = os.path.abspath(fname)
	mark = "rosetta_source/src"
	assert fname.find(mark) > 0
	return fname[fname.find(mark)+15:]


def dir2ns(fname,checkfs=True):
	"""
	get namespace corresponding to dir
	
	>>> fname = "/SOME_DIRECTORY/SOME_OTHER_DIR/rosetta_source/src/core/pack/task"
	>>> dir2ns(fname,checkfs=False)
	'core::pack::task'
	"""
	if checkfs: assert os.path.isdir(fname)
	if checkfs: fname = os.path.abspath(fname)
	fname = fname.rstrip("/")
	mark = "rosetta_source/src"
	assert fname.find(mark) > 0
	fname = fname[fname.find(mark)+len(mark)+1:]
	return fname.replace("/","::")


def swapfwdhh(fn):
	"""
	>>> print swapfwdhh("test.hh")
	test.fwd.hh
	>>> print swapfwdhh("t/e/s/t.fwd.hh")
	t/e/s/t.hh
	>>> print swapfwdhh("test.cc")
	test.cc
	"""
	if fn.endswith(".fwd.hh"): return fn.replace(".fwd.hh",".hh")
	elif fn.endswith(".hh"):   return fn.replace(".hh",".fwd.hh")
	else: return fn

if __name__ == '__main__':
	import doctest
	tr = doctest.testmod()	
	print "tests passed:",tr.attempted-tr.failed
	print "tests failed:",tr.failed
