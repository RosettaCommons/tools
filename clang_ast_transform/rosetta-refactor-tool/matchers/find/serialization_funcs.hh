// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_matchers_find_serialization_funcs_HH
#define INCLUDED_matchers_find_serialization_funcs_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

#include <set>
#include <string>
#include <utility>


class SerializationFunctionFinder : public ReplaceMatchCallback {
private:
	typedef std::set< std::string > class_names;

	class_names  classes_w_serialization_funcs_;
	bool verbose_;
public:

	SerializationFunctionFinder( clang::tooling::Replacements * replace, bool verbose );
	virtual ~SerializationFunctionFinder();

	// Main callback for all matches
	virtual void run( clang::ast_matchers::MatchFinder::MatchResult const & result);

	class_names const & classes_w_serialization_funcs() const;

};

clang::ast_matchers::DeclarationMatcher
match_to_serialization_method();


class SerializedMemberFinder : public ReplaceMatchCallback {
private:
	typedef std::set< std::pair< std::string, std::string > > data_members;
	typedef std::set< std::string > class_names;

	class_names  classes_w_serialization_funcs_;
	data_members save_variables_;
	data_members load_variables_;
	bool verbose_;
public:
	SerializedMemberFinder( clang::tooling::Replacements * replace, bool verbose = false );
	virtual ~SerializedMemberFinder();

	// Main callback for all matches
	virtual void run( clang::ast_matchers::MatchFinder::MatchResult const & result);

	class_names const & classes_w_serialization_funcs() const;
	data_members const & members_saved() const;
	data_members const & members_loaded() const;

};

clang::ast_matchers::StatementMatcher
match_to_saved_data_members();

clang::ast_matchers::StatementMatcher
match_to_loaded_data_members();


void
add_serialization_func_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );


#endif
