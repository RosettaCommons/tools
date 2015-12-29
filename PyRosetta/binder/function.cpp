// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/function.cpp
/// @brief  Binding generation for static and member functions
/// @author Sergey Lyskov

#include <function.hpp>

#include <fmt/format.h>

#include <clang/AST/DeclCXX.h>

#include <vector>


using namespace clang;

using std::string;
using std::vector;
using std::unordered_map;

using namespace fmt::literals;

namespace binder {


// Generate function argument list separate by comma: int, bool, std::sting
std::string function_arguments(clang::FunctionDecl *record)
{
	string r;

	for(uint i=0; i<record->getNumParams(); ++i) {
		r += record->getParamDecl(i)->getOriginalType().getAsString();
		if( i != record->getNumParams()-1 ) r += ", ";
	}

	return r;
}


// Generate function pointer type string for given function: void (*)(int, doule)_ or  void (ClassName::*)(int, doule)_ for memeber function
string function_pointer_type(FunctionDecl *F)
{
	string r;
	string prefix { F->isCXXClassMember() ? cast<CXXRecordDecl>( F->getParent() )->getQualifiedNameAsString() + "::" : "" };

	r += F->getReturnType().getAsString();  r+= " ({}*)("_format(prefix);

	r += function_arguments(F);

	r += ") ";

	return r;
}

// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
string bind_function(FunctionDecl *F)
{
	string function_name { F->getNameAsString() };
	string function_qualified_name { F->getQualifiedNameAsString() };

	string r = R"(.def("{}", ({}) &{}, "doc")"_format(function_name, function_pointer_type(F), function_qualified_name);


	for(uint i=0; i<F->getNumParams(); ++i) r += ", pybind11::arg(\"{}\")"_format( string( F->getParamDecl(i)->getName() ) );

	r += ')';

	return r;
}


Item bind_function(string const &module, FunctionDecl *F)
{
	Item I{ F };

	string &c(I.code);

	c += module + bind_function(F) + ';';

	return I;
}



} // namespace binder
