// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/class.cpp
/// @brief  Binding generation for C++ struct and class objects
/// @author Sergey Lyskov

#include <class.hpp>
#include <function.hpp>

#include <fmt/format.h>

//#include <vector>


using namespace llvm;
using namespace clang;

using std::string;
//using std::vector;
//using std::unordered_map;

using namespace fmt::literals;

namespace binder {

Item bind_class(clang::CXXRecordDecl *C)
{
	Item I{ C };

	string &c(I.code);

	string class_name { C->getNameAsString() };
	string class_qualified_name { C->getQualifiedNameAsString() };

	//class_<A>(module_a, "A")
	c += R"(pybind11::class_<{}>({}, "{}"))"_format(class_qualified_name, _module_variable_name_, class_name) + '\n';

	if( C->ctor_begin() == C->ctor_end() ) {  // No constructors defined, adding default constructor
		c+= "\t.def(pybind11::init<>())\n\n";
	}
	else {
		for(auto t = C->ctor_begin(); t != C->ctor_end(); ++t) {
			if( t->getAccess() == AS_public ) c+= "\t.def(pybind11::init<{}>())\n"_format( function_arguments(*t) );
		}
		c += '\n';
	}

	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( m->getAccess() == AS_public  and   !isa<CXXConstructorDecl>(*m) ) {
			c += '\t' + bind_function(*m) + '\n';
		}
	}

	c += ";\n\n";


	return I;
}



} // namespace binder
