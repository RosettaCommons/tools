// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/util.hpp
/// @brief  Various helper functions
/// @author Sergey Lyskov


#ifndef _INCLUDED_util_hpp_
#define _INCLUDED_util_hpp_

#include <clang/AST/Decl.h>
#include <clang/AST/DeclTemplate.h>

#include <string>
#include <vector>

namespace binder {


/// Split string using given separator
std::vector<std::string> split(std::string const &buffer, std::string const & separator="\n");

/// Replace all occurrences of string
std::string replace(std::string const &s, std::string const & from, std::string const &to);


/// check if string begins with given prefix
bool begins_wtih(std::string const &source, std::string const &prefix);


/// indent given code
std::string indent(std::string const &code, std::string const &indentation);


/// calculate namespace path from given NamedDecl, like: std, core::pose
std::string namespace_from_named_decl(clang::NamedDecl const *decl);


std::string typename_from_type_decl(clang::TypeDecl *decl);


/// Calculate base (upper) namespace for given one: core::pose::motif --> core::pose
std::string base_namespace(std::string const & ns);


/// Calculate last namespace for given one: core::pose::motif --> motif
std::string last_namespace(std::string const & ns);



/// replace all _Bool types with bool
void fix_boolean_types(std::string &type);

/// Generate string representation of given expression
std::string expresion_to_string(clang::Expr *e);

/// Generate string representation of given TemplateArgument
std::string template_argument_to_string(clang::TemplateArgument const &);



/// calcualte line in source file for NamedDecl
std::string line_number(clang::NamedDecl *decl);


// extract include path needed for declaration itself (without template dependency if any), return empty string if no include could be found
std::string relevant_include(clang::NamedDecl const *decl);

// extract include path needed for declaration itself (without template dependency if any), do nothing if include could not be found (ie for build-in's)
void add_relevant_include(clang::NamedDecl const *decl, std::vector<std::string> &includes);

/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(clang::QualType qt, /*clang::ASTContext const &context,*/ std::vector<std::string> &includes);



/// extract include needed for declaration and add it to includes
//bool add_relevant_includes(clang::NamedDecl *decl, std::vector<std::string> &includes);

/// Try to read exisitng file and if content does not match to code - write a new version. Also create nested dirs starting from prefix if nessesary.
void update_source_file(std::string const &prefix, std::string const &file_name, std::string const &code);



} // namespace binder

#endif // _INCLUDED_util_hpp_
