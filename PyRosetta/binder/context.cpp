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

#include <fmt/format.h>


#include <sstream>
#include <fstream>

using namespace llvm;
using namespace clang;
using std::string;
using std::vector;
using std::unordered_map;

using namespace fmt::literals;
using fmt::format;

namespace binder {

const std::string _root_module_variable_name_{"m"};


/// Split string using given separator
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


/// Replace all occurrences of string
void replace(string &s, string const & from, string const &to)
{
	size_t i = s.size();
	while( ( i = s.rfind(from, i) ) != string::npos) {
		s.replace(i, from.size(), to);
		--i;
	}
}


string indent(string const &code, string const &indentation)
{
	auto lines = split(code);
	string r;

	for(auto & l : lines) r += indentation + l + '\n';

	return r;
}


string namespace_from_named_decl(NamedDecl *decl)
{
	string qn = decl->getQualifiedNameAsString();
	string n  = decl->getNameAsString();
	//name = decl->getQualifiedNameAsString();

	int namespace_len = qn.size() - n.size();

	string path = decl->getQualifiedNameAsString().substr(0, namespace_len > 1 ? namespace_len-2 : namespace_len );  // removing trailing '::'

	return path;
}


llvm::raw_ostream & operator << (llvm::raw_ostream & os, Binder const &b)
{
	clang::NamedDecl *decl = b.get_named_decl();

	string name  = decl->getNameAsString();
	string qualified_name = decl->getQualifiedNameAsString();
	string path = decl->getQualifiedNameAsString().substr(0, qualified_name.size() - name.size() );

	return os << "B{name=" << name << ", path=" << path << ", include= , code=\n" << b("module") << "\n}\n";
}


/// Calculate base (upper) namespace for given one: core::pose::motif --> core::pose
string base_namespace(string const & ns)
{
	size_t f = ns.rfind("::");
	if( f == string::npos ) return "";
	else return ns.substr(0, f);
}


/// Calculate last namespace for given one: core::pose::motif --> motif
string last_namespace(string const & ns)
{
	size_t f = ns.rfind("::");
	if( f == string::npos ) return ns;
	else return ns.substr(f+2, ns.size()-f-2);
}

// Generate C++ name of pybind11 module variable
string module_variable(string const & namespace_)
{
	if( namespace_.size() ) {
		string m{namespace_};
		replace(m, "::", "_"); // std::replace(m.begin(), m.end(), ':', '_');
		return _root_module_variable_name_+"_" + m;
	}
	else return _root_module_variable_name_;
}

// Generate prototype/call for function binding particular module
string module_binder_function(string const & namespace_, bool prototype, string const &extra="", string const arg=_root_module_variable_name_)
{
	string r;

	if( namespace_.size() ) {
		string m{namespace_};
		replace(m, "::", "_");
		r = "bind_" + _root_module_variable_name_ + "_" + m;
	}
	else r = "bind_" + _root_module_variable_name_;

	if( extra.size() ) r += "_" + extra;

	r += "({})"_format( prototype ? "pybind11::module &" + _root_module_variable_name_ : arg );

	r = prototype ? "void " + r : '\t' + r;

	return r;
}


void Context::add(BinderOP &b)
{
	string ns = namespace_from_named_decl( b->get_named_decl() );
	modules[ns].push_back(b);
	//outs() << "Adding to: " << ns << "\n";
}


/// recursivly create all nested namespaces in case some is missings
void Context::create_all_nested_namespaces()
{
	vector<string> n;

	for(auto & p : modules) {
		string ns = p.first;

		while( ns.size() ) {
			ns = base_namespace(ns);
			n.push_back(ns);
		}
	}

	for(auto & i : n) modules[i];
	modules[""];
}

/// create vector of all namespaces and sort it
vector<string> Context::sorted_namespaces()
{
	vector<string> n;
	for(auto & p : modules) n.push_back(p.first);

	std::sort( std::begin(n), std::end(n) );

	return n;
}


/// generate bindings for namespace and split it in chunks so maximum size is ~ max_code_size
vector<string> Context::bind_namespaces(string const &namespace_, size_t max_code_size)
{
	vector<BinderOP> & binders( modules[namespace_] );

	string prototypes;

	vector<string> R;
	string r;

	for(auto &b : binders) {
		if( r.size() == 0 ) {
			r += module_binder_function(namespace_, true, std::to_string( R.size() ) );
			prototypes += r + ";\n";
			r += " {\n";
		}

		r += (*b)(_root_module_variable_name_, "\t");

		if( r.size() >= max_code_size ) { r+= "}\n"; R.push_back(r); r=""; }
	}

	if( r.size() ) { r+= "}\n"; R.push_back(r); }

	string module_main = "\n\n" + prototypes + "\n" + module_binder_function(namespace_, true) + " {\n";
	for(size_t i=0; i<R.size(); ++i) {
		module_main += module_binder_function(namespace_, false, std::to_string(i) ) + ";\n"; // calls to binder functions
	}
	module_main += "}\n";

	if( R.size() ) R.back() += module_main;
	else R.push_back(module_main);

	return R;
}



const char * module_header = R"_(#include <pybind11/pybind11.h>
#include <self_test.hpp>

{}
PYBIND11_PLUGIN(example) {{
	pybind11::module {}("{}", "{}");

)_";

void Context::generate()
{

	//s << "#include <pybind11/pybind11.h>\nPYBIND11_PLUGIN(example) {\n\tpybind11::module m(\"example\", \"example module\");\n";

	create_all_nested_namespaces();
	vector<string> namespaces = sorted_namespaces();

	string module_var_defs, module_function_defs, module_function_calls;

	for(auto & n : namespaces) {
		if( n.size() ) {
			module_var_defs += "\tpybind11::module {} = {}.def_submodule(\"{}\", \"Bindings for {} namespace\");\n"_format( module_variable(n), module_variable( base_namespace(n) ), last_namespace(n), n);
		}

		module_function_defs += module_binder_function(n, true) + ";\n";
		module_function_calls += module_binder_function(n, false, "", module_variable(n)) + ";\n";
	}

	std::stringstream s;
	s << format(module_header, module_function_defs, _root_module_variable_name_, "example", "example module");

	s << module_var_defs << "\n" << module_function_calls;
	s << "\treturn m.ptr();\n}\n";

	//for(auto & p : modules) outs() << "Ns: " << p.first << "\n";

	// vector<BinderOP> & binders( modules[""] );
	// //for(auto &b : binders) s << indent(b("m"), "\t");
	// for(auto &b : binders) s << (*b)("m", "\t");
	// s << "\treturn m.ptr();\n}\n";

	for(auto &ns : modules) {
		vector<string> files = bind_namespaces(ns.first, 1024*256);
		for(auto &f : files) s << f;
	}

	errs() << s.str() << "\n";

	std::ofstream f("../tools/clang/tools/extra/binder/_example_/example.cpp");
	f << s.str();
}


} // namespace binder
