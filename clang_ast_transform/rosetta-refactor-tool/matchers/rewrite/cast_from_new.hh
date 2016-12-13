// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_cast_from_new_HH
#define INCLUDED_matchers_rewrite_cast_from_new_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace implicit casts in "new" instantiation:

		ClassAOP new_aap() { return new ClassA; }
		AtomCOP foo() { return this; }
		foo(AtomCOP a); foo(0);

	variables_[ varname ] = new VariableExpression( varname );
*/

class RewriteCastFromNew : public ReplaceMatchCallback {
public:
	RewriteCastFromNew(
		clang::tooling::Replacements * Replace,
		const char *tag = "CastFromNew");

	RewriteCastFromNew(
		clang::tooling::Replacements * Replace,
		bool debug,
		const char *tag = "CastFromNew");

	void run(const clang::ast_matchers::MatchFinder::MatchResult & Result);
};

void
add_cast_from_new_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif

