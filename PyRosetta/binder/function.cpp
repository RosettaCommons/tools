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

#include <vector>


using namespace clang;

using std::string;
using std::vector;
using std::unordered_map;

using namespace fmt::literals;

namespace binder {


// Generate function argument list separate by comma
std::string function_arguments(clang::FunctionDecl *record)
{
	string r;

	for(uint i=0; i<record->getNumParams(); ++i) {
		r += record->getParamDecl(i)->getOriginalType().getAsString();
		if( i != record->getNumParams()-1 ) r += ", ";
	}

	return r;
}


// Generate function pointer type string for given function. Example void (*)(int, doule)_ or  void (ClassName::*)(int, doule)_ for memeber function
string function_pointer_type(FunctionDecl *record)
{
	string r;
	r += record->getReturnType().getAsString();  r+= " (*)(";

	r += function_arguments(record);

	r += ") ";

	return r;
}


Item bind_function(FunctionDecl *R)
{
	Item I{ R };

	string &c(I.code);

	string function_name { R->getNameAsString() };
	string function_qualified_name { R->getQualifiedNameAsString() };

	//c += _module_variable_name_ + ".def(\"" + function_name + "\", " + function_pointer_type(R) + "&" + function_name + ", \"doc\");\n";
	c += _module_variable_name_ + R"(.def("{}", ({}) &{}, "doc");)"_format(function_name, function_pointer_type(R), function_qualified_name);

	return I;
}



} // namespace binder
