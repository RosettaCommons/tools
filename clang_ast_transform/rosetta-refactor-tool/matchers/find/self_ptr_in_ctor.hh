// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_find_self_ptr_in_ctor_HH
#define INCLUDED_matchers_find_self_ptr_in_ctor_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Find instances where get_self_ptr() or get_self_weak_ptr() is used in a c'tor.
	This is illegal because the weak self-pointer isn't set yet, and will
	result in bad_weak_ptr exception at runtime.
*/

class SelfPtrInCtorFinder : public ReplaceMatchCallback {
public:
	SelfPtrInCtorFinder(
		clang::tooling::Replacements *Replace,
		const char *tag = "SelfPtrInCtorFinder");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);
};

void
add_self_ptr_in_ctor_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
