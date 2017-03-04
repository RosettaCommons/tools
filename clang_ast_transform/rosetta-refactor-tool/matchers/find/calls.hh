// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/*
	Find all function calls in class methods
	To get a parsable list, run:
	./run.sh test-cases.cc -matchers=find_calls -colors=false|grep ' => '|sed 's/ => /\t/'
*/
#ifndef INCLUDED_matchers_find_calls_HH
#define INCLUDED_matchers_find_calls_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

class CallsFinder : public ReplaceMatchCallback {
public:
	CallsFinder(
		clang::tooling::Replacements *Replace,
		const char *tag = "CallsFinder");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

};

void
add_calls_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
