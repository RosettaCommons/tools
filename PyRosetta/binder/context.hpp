// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/context.hpp
/// @brief  Data structures to represent root context and modules
/// @author Sergey Lyskov

#ifndef _INCLUDED_context_hpp_
#define _INCLUDED_context_hpp_

#include <clang/AST/Decl.h>

#include <llvm/Support/raw_ostream.h>

#include <string>
#include <unordered_map>


namespace binder {

/// Bindings Generator - represent object that can generate binding info for function, class, enum or data variable
class Binder
{
public:
	typedef std::string string;

	virtual ~Binder() {}

	/// check if generator can create binding
	virtual bool is_bindable() const = 0;

	/// generate binding code
	virtual string operator()(string const &module_variable_name, string const &indentation="\t") const = 0;

	virtual clang::NamedDecl * get_named_decl() const = 0;
};


typedef std::shared_ptr< Binder > BinderOP;

typedef std::vector<BinderOP> Binders;


// struct Item_
// {
// 	typedef std::string string;

// 	Item_(clang::NamedDecl *decl);

// 	string name; // full C++ name of this element including namespace
// 	string path; // path for Python module for this namespace


// 	//vector<string> includes;
// 	string include;  // name of header file that needed to be include in order to compile this bindings

// 	std::vector<string> dependencies;  // list of C++ names that need to be defined/binded before this item could be compiled

// 	string code;  // C++ code that create bindins for this item
// };


llvm::raw_ostream & operator << (llvm::raw_ostream & os, Binder const &b);

/// Module - represent bindings of individual Python module
// struct Module
// {
// 	std::string path;
// 	std::vector<BinderOP> binders;
// };


/// Context - root, hold bindings info for whole TranslationUnit
class Context
{
	typedef std::string string;

public:

	void add(BinderOP &);

	void generate(std::string const &root_module, std::string const &prefix, int maximum_file_length);

private:
	std::unordered_map<string, Binders> modules;

	std::unordered_map<string, BinderOP> system_binders;

	void create_all_nested_namespaces();

	/// create vector of all namespaces and sort it
	std::vector<string> sorted_namespaces();
	std::vector<string> bind_namespaces(string const &namespace_, size_t max_code_size);
};




} // namespace binder

#endif // _INCLUDED_context_hpp_
