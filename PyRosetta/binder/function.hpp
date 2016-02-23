// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/function.hpp
/// @brief  Binding generation for static and member functions
/// @author Sergey Lyskov

#ifndef _INCLUDED_function_hpp_
#define _INCLUDED_function_hpp_

#include <context.hpp>

#include <clang/AST/Decl.h>

#include <string>


namespace binder {


// Generate function argument list separate by comma
std::string function_arguments(clang::FunctionDecl *record);


// Generate function pointer type string for given function. Example void (*)(int, doule)_ or  void (ClassName::*)(int, doule)_ for memeber function
std::string function_pointer_type(clang::FunctionDecl *record);


// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
std::string bind_function(clang::FunctionDecl *F);


//Item bind_function(std::string const &module, clang::FunctionDecl *F);


class FunctionBinder : public Binder
{
public:
	FunctionBinder(clang::FunctionDecl *f) : F(f) {}

	/// check if generator can create binding
	bool bindable() const override;

	/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
	string id() const override;

	/// generate binding code for this object and all its dependencies
	void bind(Context &) override;

    clang::NamedDecl * named_decl() const override { return F; };

	// check if bindings for object should be skipped
	bool is_skipping_requested(std::vector<std::string> const & namespaces_to_skip) const override { return false; }


private:
	clang::FunctionDecl *F;
};


/// check if generator can create binding
bool is_bindable(clang::FunctionDecl const *F);



} // namespace binder

#endif // _INCLUDED_function_hpp_
