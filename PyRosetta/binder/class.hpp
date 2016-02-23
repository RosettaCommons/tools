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
#include <vector>


namespace binder {

class ClassBinder : public Binder
{
public:
	ClassBinder(clang::CXXRecordDecl *c) : C(c) {}


	/// check if generator can create binding
    bool bindable() const override;

	/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
	string id() const override;

	/// generate binding code for this object and all its dependencies
	void bind(Context &) override;

	clang::NamedDecl * named_decl() const override { return C; };

	// check if bindings for object should be skipped
	bool is_skipping_requested(std::vector<std::string> const & namespaces_to_skip) const override;


private:
	clang::CXXRecordDecl *C;
};


// generate classtemplate specialization for ClassTemplateSpecializationDecl or empty string otherwise
std::string template_specialization(clang::CXXRecordDecl const *C);


// generate class name that could be used in bindings code indcluding template specialization if any
std::string class_name(clang::CXXRecordDecl *C);


// generate qualified class name that could be used in bindings code indcluding template specialization if any
std::string class_qualified_name(clang::CXXRecordDecl *C);


/// check if generator can create binding
bool is_bindable(clang::CXXRecordDecl const *C);


/// check if bindings for particular object was requested
bool is_binding_requested(clang::CXXRecordDecl const *C, std::vector<std::string> const & namespaces_to_bind);

// check if bindings for object should be skipped
bool is_skipping_requested(clang::CXXRecordDecl const *C, std::vector<std::string> const & namespaces_to_skip);



} // namespace binder

#endif // _INCLUDED_class_hpp_
