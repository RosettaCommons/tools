// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_real_comparison_HH
#define INCLUDED_matchers_rewrite_real_comparison_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace comparison of Reals (doubles or floats):
	double a, b;
	if(a == b) {}
	if(a != b) {}
*/

class RewriteRealComparisons : public ReplaceMatchCallback {
public:
	RewriteRealComparisons(
		clang::tooling::Replacements *Replace,
		const char *tag,
		const char *fn);

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);
private:
	std::string fn_;
};

void
add_real_comparison_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
