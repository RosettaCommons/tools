#!/usr/bin/env python

import os
import re
import subprocess
import sys

class Clang:

	def __init__(self, clang_path, include_paths):
	
		self.clang_path = clang_path
		self.include_paths = include_paths
	
	def ast_xml(self, filename):
	
		args = [self.clang_path, "-cc1", filename]
		args += ["-I" + path for path in self.include_paths]
		args.append("-ast-print-xml")
		args.append("-fexceptions")
		
		return subprocess.call(args)

class ClangXML:

	def __init__(self, clang_xml_filename):
		
		clang_xml_file = open(clang_xml_filename)
		self.lines = clang_xml_file.readlines()
		self.access_dict = None
		clang_xml_file.close()

	def ids_for_file(self, filepath):
	
		filere = re.compile('      <File id="([^"]+)" name="' + filepath + '"\/>')
		
		return [filere.sub("\\1", line.rstrip()) for line in self.lines if filere.match(line)]
	
	def ids_for_var(self, var_name, file_id):
	
		varre = re.compile('.+<Var id="([^"]+)" file="' + file_id + '".+name="' + var_name +'".+/>')
		
		return [varre.sub("\\1", line.rstrip()) for line in self.lines if varre.match(line)]
	
	def ids_for_method(self, var_name, file_id):
	
		methodre = re.compile('.+<CXXMethod id="([^"]+)" file="' + file_id + '".+name="' + var_name +'".+/>')
		
		return [methodre.sub("\\1", line.rstrip()) for line in self.lines if methodre.match(line)]

	def linenums_for_ref(self, ref_id):
	
		refre = re.compile('ref="' + ref_id +'".+>')
		
		return [i for i in xrange(len(self.lines)) if refre.search(self.lines[i])]

	def linenum_for_id(self, id):
	
		idre = re.compile('id="' + id +'".+/*>')
		
		for i in xrange(len(self.lines)):
			if idre.search(self.lines[i]):
				return i

	def parent_linenum(self, linenum, parent_types):
	
		linespaces = re.sub("^( +).+$", "\\1", self.lines[linenum].rstrip())
		parentre = re.compile("^ {0," + str(len(linespaces)-1) + "}<(" + "|".join(parent_types) + ")")
		
		for i in xrange(linenum-1, -1, -1):
			if parentre.match(self.lines[i]):
				return i

	def line_tag(self, linenum):
		
		return re.sub(" +<([^ ]+).+", "\\1", self.lines[linenum].rstrip())
	
	def line_attributes(self, linenum):
	
		attrs = re.sub(" +<[^ ]+ (.+)/*>", "\\1", self.lines[linenum].rstrip())
		
		attr_dict = {}
		attrre = re.compile("([^=]+)=\"(.+)\"")
		#print self.lines[linenum],
		for attr in attrs.split(" "):
			#print attr
			attr_match = attrre.match(attr)
			attr_dict[attr_match.group(1)] = attr_match.group(2)

		return attr_dict

	def line_parents(self, linenum):
	
		spacere = re.compile("^( +)")
		numspaces = len(spacere.match(self.lines[linenum]).group(1))
		parents = []
		
		last_line = linenum
		for i in xrange(linenum-1, -1, -1):
			line_numspaces = len(spacere.match(self.lines[i]).group(1))
			if line_numspaces < numspaces:
				if self.line_tag(i) == "Contexts":
					last_context = self.line_attributes(last_line)["context"]
					new_context_linenum = self.linenum_for_id(last_context)
					parents.append(self.line_attributes(new_context_linenum)["name"])
					parents.reverse()
					return self.line_parents(new_context_linenum) + parents
				#print self.lines[i],
				parents.append(self.line_attributes(i)["name"])
				last_line = i;
				numspaces = line_numspaces
			if numspaces <= 4:
				break
		
		parents.reverse()
		
		return parents

	def range_for_block(self, linenum):
	
		if self.lines[linenum].rstrip().endswith("/>"):
			return (linenum, linenum+1)
		
		spacere = re.compile("^( +)")
		numspaces = len(spacere.match(self.lines[linenum]).group(1))
		
		if self.lines[linenum].lstrip().startswith("</"):
			for i in xrange(linenum-1, -1, -1):
				line_numspaces = len(spacere.match(self.lines[i]).group(1))
				if numspaces == line_numspaces:
					return (i, linenum+1)
		else:
			for i in xrange(linenum+1, len(self.lines)):
				line_numspaces = len(spacere.match(self.lines[i]).group(1))
				if numspaces == line_numspaces:
					return (linenum, i+1)

	def ref_ids_for_range(self, start_linenum, end_linenum):
		
		refre = re.compile(' ref="([^"]+)"')
		
		refs = []
		
		for line in self.lines[start_linenum:end_linenum]:
		
			refmatch = refre.search(line)
			if refmatch:
				refs.append(refmatch.group(1))

		return refs

	def access_for_id(self, id):
	
		if self.access_dict is None:
			self.access_dict = {}
			accessre = re.compile('id="([^"]+)".+access="([^"]+)"')
			for line in self.lines:
				accessmatch = accessre.search(line)
				if accessmatch:
					self.access_dict[accessmatch.group(1)] = accessmatch.group(2)

		return self.access_dict.get(id)

	def refs_for_id(self, id):
	
		ref_list = []
		linenum_list = []
		count_list = []
	
		ref_linenums = self.linenums_for_ref(id)
	
		for linenum in ref_linenums:
		
			parent_linenum =  self.parent_linenum(linenum, ["CXXMethod", "CXXConstructor", "Function"])
			#print "p: " + self.lines[parent_linenum],
			try:
				index = linenum_list.index(parent_linenum)
				count_list[index] += 1
				continue
			except:
				pass
			parent_attrs = self.line_attributes(parent_linenum)
			parent_file_linenum = self.linenum_for_id(parent_attrs["file"])
			parent_file_attrs = self.line_attributes(parent_file_linenum)
			
			context_line = self.linenum_for_id(parent_attrs["context"])
			#print "c: " + self.lines[context_line],
			if self.line_tag(context_line) == "TranslationUnit":
				context_parents = []
			else:
				context_parents = self.line_parents(context_line)
				context_parents.append(self.line_attributes(context_line)["name"])
			context_parents.append(parent_attrs["name"])
			
			parent_file_name = parent_file_attrs["name"]
			if parent_file_name.startswith("./"):
				parent_file_name = parent_file_name[2:]
			ref_line = parent_file_name + ":" + parent_attrs["line"] + " - " + "::".join(context_parents) + " (" + self.line_tag(parent_linenum)

			ref_list.append(ref_line)
			linenum_list.append(parent_linenum)
			count_list.append(1)
		
		for i in xrange(len(ref_list)):
			parent_range = self.range_for_block(linenum_list[i])
			parent_refs = self.ref_ids_for_range(parent_range[0], parent_range[1])
			
			public_refs = set()
			protected_refs = set()
			private_refs = set()
				
			for ref in parent_refs:
				access_type = self.access_for_id(ref)
				if access_type == "public":
					public_refs.add(ref)
				if access_type == "protected":
					protected_refs.add(ref)
				if access_type == "private":
					private_refs.add(ref)
			
			if False:
				print "Public:"
				for ref in public_refs:
					print self.lines[self.linenum_for_id(ref)]
				print "Protected:"
				for ref in protected_refs:
					print self.lines[self.linenum_for_id(ref)]
				print "Private:"
				for ref in private_refs:
					print self.lines[self.linenum_for_id(ref)]
			
			ref_list[i] += ", option: %i, public: %i, protected: %i, private: %i)" % (count_list[i], len(public_refs), len(protected_refs), len(private_refs))
		
		return ref_list

def makepaths(files):

	paths = []
	
	sortedkeys = files.keys()
	sortedkeys.sort()
	for key in sortedkeys:
		value = files[key]
		if type(value) == type([]):
			newpaths = [valueitem.lstrip("/") + ".cc" for valueitem in value]
		elif type(value) == type({}):
			newpaths = makepaths(value)
		else:
			print type(value)
			continue
		
		paths += [os.path.join(key, path) for path in newpaths]

	return paths

if __name__ == "__main__":

	clang = Clang("clang", [".", "platform/linux", "../external/include", "../external/boost_1_38_0"])

	refs_list = []
	refs_set = set()

	filepaths = []
	for library in ["core", "protocols", "devel"]:
		files = {}
		execfile(library + ".src.settings", files)
		filepaths += [os.path.join(library, path) for path in makepaths(files["sources"])]
	
	#print "\n".join(filepaths)

	#filepaths = ["core/chemical/ResidueType.cc"]

	for filepath in filepaths:
	
		print filepath
	
		xmlpath = os.path.splitext(filepath)[0] + ".xml"
		
		if True or not os.path.isfile(xmlpath):
			retcode = clang.ast_xml(filepath)
			if retcode != 0:
				continue
		clang_xml = ClangXML(xmlpath)
		os.remove(xmlpath)
	
		option_file_ids = clang_xml.ids_for_file("./core/options/option.hh")
		if len(option_file_ids) == 0:
			continue
		option_var_ids = clang_xml.ids_for_var("option", option_file_ids[0])
		
		new_refs_list = clang_xml.refs_for_id(option_var_ids[0])
		
		for ref in new_refs_list:
			if not ref in refs_set:
				#print ref
				refs_list.append(ref)
				refs_set.add(ref)
		
	print ""
	print "Unique core::options::option[] accessors:"
	print "\n".join(refs_list)
