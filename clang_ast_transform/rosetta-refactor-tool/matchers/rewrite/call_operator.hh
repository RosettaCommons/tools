// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_call_operator_HH
#define INCLUDED_matchers_rewrite_call_operator_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace calls to operator() on owning_ptrs which return a named pointer to the object

		SomeOP a;
		foo( a() );
		foo( *a() );

	TODO:
	Replace with a.get() on OPs and a.lock().get() instead? Do we want this?
*/

class RewriteCallOperator : public ReplaceMatchCallback {
public:
	RewriteCallOperator(
		clang::tooling::Replacements *Replace,
		const char *tag = "CallOperator");

	RewriteCallOperator(
		clang::tooling::Replacements *Replace,
		bool debug,
		const char *tag = "CallOperator");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult & Result);
};

void
add_call_operator_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
