// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_cast_from_new_vardecl_HH
#define INCLUDED_matchers_rewrite_cast_from_new_vardecl_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace implicit casts in "new" instantiation in variable declarations:

	FooOP foo( new Foo );
	FooOP foo = new Foo;
*/

class RewriteCastFromNewVarDecl : public ReplaceMatchCallback {
public:
	RewriteCastFromNewVarDecl(
		clang::tooling::Replacements *Replace,
		const char *tag = "CastFromNewVarDecl");

	RewriteCastFromNewVarDecl(
		clang::tooling::Replacements *Replace,
		bool debug,
		const char *tag = "CastFromNewVarDecl");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);
};

void
add_cast_from_new_vardecl_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
