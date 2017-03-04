// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_pointer_name_HH
#define INCLUDED_matchers_rewrite_pointer_name_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace owning_ptr/access_ptr name in templates, etc.
*/

class RewritePointerName : public ReplaceMatchCallback {
	
public:
	RewritePointerName(
		clang::tooling::Replacements *Replace,
		const char *tag ="PointerName");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);
};

void
add_pointer_name_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif



