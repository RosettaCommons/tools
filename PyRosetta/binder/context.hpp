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
#include <vector>
#include <set>
#include <unordered_map>


namespace binder {

class Context;


/// Bindings Generator - represent object that can generate binding info for function, class, enum or data variable
class Binder
{
public:
	typedef std::string string;

	virtual ~Binder() {}

	/// check if generator can create binding
	virtual bool bindable() const = 0;

	/// generate binding code for this object and all its dependencies
	virtual void bind(Context &) = 0;

	virtual clang::NamedDecl * named_decl() const = 0;

	/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
	virtual string id() const = 0;

	// check if bindings for object should be skipped
	virtual bool is_skipping_requested(std::vector<std::string> const & namespaces_to_skip) const = 0;

	// return true if code was already generate for this object
	bool is_binded() const { return code_.size(); }

	// return binding code
	string & code() { return code_; }

	/// return true if object declared in system header
	bool is_in_system_header();

	explicit operator std::string() const { return id(); /*named_decl()->getQualifiedNameAsString();*/ }

private:
	string code_;
};


typedef std::shared_ptr< Binder > BinderOP;

typedef std::vector<BinderOP> Binders;



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

	void generate(string const &root_module, std::vector<string> const & namespaces_to_bind, std::vector<string> const & namespaces_to_skip, std::string const &prefix, uint maximum_file_length);

	void bind(std::string const &object);

	/// generate C++ expression for module variable for namespace_
	string module_variable_name(string const & namespace_);

	//void request_bindings(std::vector<string> const & namespaces);
	// is binding requested

private:

	/// array of all binderes from translation unit
	std::vector<BinderOP> binders;

	/// types → binder
	std::unordered_map<string, BinderOP> types;

	std::set<string> ids;

	// binder.id() → binder
	//std::unordered_map<string, Binders> binders_map;

	std::unordered_map<string, Binders> modules;
	std::unordered_map<string, BinderOP> system_binders;


	/// bind all objects residing in namespaces and it dependency
	void bind(std::vector<string> const & namespaces_to_bind, std::vector<string> const & namespaces_to_skip);

	std::set<string> create_all_nested_namespaces();

	/// create vector of all namespaces and sort it
	//std::vector<string> sorted_namespaces();
	std::vector<string> bind_namespaces(string const &namespace_, size_t max_code_size);
};




} // namespace binder

#endif // _INCLUDED_context_hpp_
