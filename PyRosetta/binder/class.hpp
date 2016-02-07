// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/class.hpp
/// @brief  Binding generation for C++ struct and class objects
/// @author Sergey Lyskov

#ifndef _INCLUDED_class_hpp_
#define _INCLUDED_class_hpp_

#include <context.hpp>

#include <clang/AST/DeclCXX.h>

#include <string>


namespace binder {

class ClassBinder : public Binder
{
public:
	ClassBinder(clang::CXXRecordDecl *c) : C(c) {}


	/// check if generator can create binding
	//bool is_bindable() const override;


	/// generate binding code
	string operator()(string const &module_variable_name, string const &indentation="\t") const override;


	clang::NamedDecl * get_named_decl() const override { return C; };

private:
	clang::CXXRecordDecl *C;
};


// generate class name that could be used in bindings code indcluding template specialization if any
std::string class_name(clang::CXXRecordDecl *C);


// generate qualified class name that could be used in bindings code indcluding template specialization if any
std::string class_qualified_name(clang::CXXRecordDecl *C);


/// check if generator can create binding
bool is_bindable(clang::CXXRecordDecl *C);


//Item bind_class(clang::CXXRecordDecl *R);



} // namespace binder

#endif // _INCLUDED_class_hpp_
