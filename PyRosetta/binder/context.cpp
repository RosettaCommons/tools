// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/context.cpp
/// @brief  Data structures to represent root context and modules
/// @author Sergey Lyskov


#include <context.hpp>

#include <sstream>
#include <fstream>


using namespace llvm;
using namespace clang;
using std::string;
using std::vector;
using std::unordered_map;


namespace binder {

vector<string> split(string const &buffer, string const & separator="\n")
{
	string line;
	vector<string> lines;

	for(uint i=0; i<buffer.size(); ++i) {
		if( buffer.compare(i, separator.size(), separator) ) line.push_back( buffer[i] );
		else {
			lines.push_back(line);
			line.resize(0);
		}
	}

	if( line.size() ) lines.push_back(line);

	return lines;
}


string indent(string const &code, string const &indentation)
{
	auto lines = split(code);
	string r;

	for(auto & l : lines) r += indentation + l + '\n';

	return r;
}



Item::Item(NamedDecl *decl)
{
	string qn = decl->getQualifiedNameAsString();
	string n  = decl->getNameAsString();

	name = decl->getQualifiedNameAsString();
	path = decl->getQualifiedNameAsString().substr(0, qn.size() - n.size() );
}



llvm::raw_ostream & operator << (llvm::raw_ostream & os, Item const &i)
{
	return os << "I{name=" << i.name << ", path=" << i.path << ", include=" << i.include << ", code=\n" << i.code << "\n}\n";
}



void Context::add(Item const &I)
{
	modules[""].items.push_back(I);
}


const char * module_header = R"_(#include <pybind11/pybind11.h>
#include <self_test.hpp>

PYBIND11_PLUGIN(example) {
	pybind11::module m("example", "example module");

)_";

void Context::generate()
{
	std::stringstream s;

	//s << "#include <pybind11/pybind11.h>\nPYBIND11_PLUGIN(example) {\n\tpybind11::module m(\"example\", \"example module\");\n";

	s << module_header;

	vector<Item> & items( modules[""].items );

	for(auto &i : items) s << indent(i.code, "\t");

	s << "\n\treturn m.ptr();\n}\n";

	errs() << s.str() << "\n";

	std::ofstream f("../tools/clang/tools/extra/binder/_example_/example.cpp");
	f << s.str();
}


} // namespace binder
