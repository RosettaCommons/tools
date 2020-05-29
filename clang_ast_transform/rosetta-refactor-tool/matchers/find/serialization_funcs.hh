// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_find_serialization_funcs_HH
#define INCLUDED_matchers_find_serialization_funcs_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

#include <set>
#include <string>
#include <utility>


class SerializationFunctionFinder : public ReplaceMatchCallback {
public:
	typedef std::set< std::string > class_names;
	typedef std::set< std::pair< std::string, std::string > > data_members;

	SerializationFunctionFinder( clang::tooling::Replacements * replace, bool verbose );
	virtual ~SerializationFunctionFinder();

	// Main callback for all matches
	virtual void run( clang::ast_matchers::MatchFinder::MatchResult const & result);

	class_names const & classes_w_serialization_funcs() const;

	data_members const & exempted_members_from_save() const;
	data_members const & exempted_members_from_save_w_opts() const;
	data_members const & exempted_members_from_load() const;
	data_members const & exempted_members_from_load_w_opts() const;

private:

	class_names  classes_w_serialization_funcs_;
	data_members save_members_exempted_;
	data_members save_w_opts_members_exempted_;
	data_members load_members_exempted_;
	data_members load_w_opts_members_exempted_;
	bool verbose_;
};

clang::ast_matchers::DeclarationMatcher
match_to_serialization_method_definition();


class SerializedMemberFinder : public ReplaceMatchCallback {
public:
	typedef std::set< std::pair< std::string, std::string > > data_members;
	typedef std::set< std::string > class_names;

	SerializedMemberFinder( clang::tooling::Replacements * replace, bool verbose = false );
	virtual ~SerializedMemberFinder();

	// Main callback for all matches
	virtual void run( clang::ast_matchers::MatchFinder::MatchResult const & result);

	class_names const & classes_w_serialization_funcs() const;
	data_members const & members_saved() const;
	data_members const & members_saved_w_opts() const;
	data_members const & members_loaded() const;
	data_members const & members_loaded_w_opts() const;

private:
	class_names  classes_w_serialization_funcs_;
	data_members save_variables_;
	data_members save_w_opts_variables_;
	data_members load_variables_;
	data_members load_w_opts_variables_;
	bool verbose_;
};

clang::ast_matchers::StatementMatcher
match_to_serialized_data_members();

clang::ast_matchers::StatementMatcher
match_to_externally_serialized_data_members();

void
add_serialization_func_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );


#endif
