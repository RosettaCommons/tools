// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_code_quality_naked_ptr_op_casts_HH
#define INCLUDED_matchers_code_quality_naked_ptr_op_casts_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Find naked pointer to OP/AP casts:
		SomeOP op;
		Some *x = op;
		Some &x = *op;
		SomeOP op = x;
		SomeOP op = &x;
*/

class CastFinder : public ReplaceMatchCallback {
public:
	CastFinder(
		clang::tooling::Replacements *Replace,
		const char *tag = "CastFinder");

	CastFinder(
		clang::tooling::Replacements *Replace,
		bool debug,
		const char *tag = "CastFinder");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);
};

void
add_naked_ptr_op_casts_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif

