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


class SerializationFuncFinder : public ReplaceMatchCallback {
private:
	typedef std::set< std::pair< std::string, std::string > > data_members;

	data_members save_variables_;
	data_members load_variables_;

public:
	SerializationFuncFinder( clang::tooling::Replacements * replace );
	virtual ~SerializationFuncFinder();

	// Main callback for all matches
	virtual void run( clang::ast_matchers::MatchFinder::MatchResult const & result);

};

void
add_serialization_func_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );


#endif
