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

//const string _module_variable_name_{"~!@#module_variable_name#@!~"};
const std::string _module_variable_name_{"m"};

/// Item - replresent binding info for function, class, enum or variable
struct Item
{
	typedef std::string string;

	Item(clang::NamedDecl *decl);

	string name; // full C++ name of this element including namespace
	string path; // path for Python module for this namespace


	//vector<string> includes;
	string include;  // name of header file that needed to be include in order to compile this bindings

	std::vector<string> dependencies;  // list of C++ names that need to be defined/binded before this item could be compiled

	string code;  // C++ code that create bindins for this item
};


llvm::raw_ostream & operator << (llvm::raw_ostream & os, Item const &i);

/// Module - represent bindings of individual Python module
struct Module
{
	std::string path;

	std::vector<Item> items;
};


/// Context - root, hold bindings info for whole TranslationUnit
struct Context
{
	typedef std::string string;


	std::unordered_map<string, Module> modules;

	std::unordered_map<string, Item> system_items;

	void add(Item const &);

	void generate();
};



} // namespace binder

#endif // _INCLUDED_context_hpp_
