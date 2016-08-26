// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_find_field_decl_HH
#define INCLUDED_matchers_find_field_decl_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Find instances where get_self_ptr() or get_self_weak_ptr() is used in a c'tor.
	This is illegal because the weak self-pointer isn't set yet, and will
	result in bad_weak_ptr exception at runtime.
*/

class FieldDeclaration {
public:
	std::string full_name_;
	std::string var_name_;
	std::string type_;
	std::string type_desugared_;
	std::string cls_;
	std::string loc_;
	clang::CharSourceRange range_;
	unsigned int start_;
	unsigned int end_;
};

class FieldDeclFinder : public ReplaceMatchCallback {

public:
	FieldDeclFinder(clang::tooling::Replacements *Replace);
	FieldDeclFinder(clang::tooling::Replacements *Replace, bool verbose, bool match_target_file_only );

	// Main callback for all matches
	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

	std::list< std::string >
	all_classes_with_fields() const;

	std::list< FieldDeclaration >
	fields_for_class( std::string const & classname ) const;

private:
	bool verbose_;
	bool match_target_file_only_;
	std::map< std::string, std::list< FieldDeclaration > > class_fields_;

};

clang::ast_matchers::DeclarationMatcher
match_to_field_decl();

void
add_field_decl_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
