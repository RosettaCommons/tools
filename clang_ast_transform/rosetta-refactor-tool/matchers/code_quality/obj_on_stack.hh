// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_code_quality_obj_on_stack_HH
#define INCLUDED_matchers_code_quality_obj_on_stack_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Code quality checker finder:
	- Find locations where an object that uses enable_shared_from_this<>
	  is declared on the stack; see matchers at end of this file
	  for a list of classes as of Sept 2014.

	Example:

class X : public std::enable_shared_from_this< X > {
public:
	X() {}
	~X() {}
};

void foo() {
	X x;          // BAD
	X & x_r = x;  // OK
}

*/

class ObjOnStackFinder : public ReplaceMatchCallback {
public:
	ObjOnStackFinder(
		clang::tooling::Replacements *Replace,
		const char *tag = "ObjOnStackFinder");

	ObjOnStackFinder(
		clang::tooling::Replacements *Replace,
		bool verbose,
		const char *tag = "ObjOnStackFinder");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult & Result);

private:
	bool verbose_;
};

void
add_obj_on_stack_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
