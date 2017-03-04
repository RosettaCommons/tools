// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_typedef_HH
#define INCLUDED_matchers_rewrite_typedef_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

class RewriteDecl : public ReplaceMatchCallback {

public:
	RewriteDecl(
		clang::tooling::Replacements *Replace,
		const char *tag ="Decl",
		const char *delim ="");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

private:
	const char *delim;
};

void
add_typedef_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
