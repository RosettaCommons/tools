// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_ctor_initializer_HH
#define INCLUDED_matchers_rewrite_ctor_initializer_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace 0 and NULL in constructor initializers:

	class ClassA {
		ClassA() : op1_( 0 ), ap1_( 0 ), op2_( NULL ), ap_( NULL ) { }
		ClassAOP op1_, op2_;
		ClassAAP ap1_, ap2_;
	}

	Note:
	std::weak_ptr and std::shard_ptr default c'tor initializes the object to null.
	While std::shard_ptr can be initialized with 0 or NULL as written,
	this will not work for std::weak_ptr because it can only by initialized with a shared_ptr.
*/

////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace 0 in constructor initializers

class RewriteCtorInitializer : public ReplaceMatchCallback {

public:
	RewriteCtorInitializer(
		clang::tooling::Replacements *Replace,
		const char *tag ="CtorInitializer");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);
};

class RewriteCtorInitializerHacky : public ReplaceMatchCallback {

public:
	RewriteCtorInitializerHacky(
		clang::tooling::Replacements *Replace,
		const char *tag ="RewriteCtorInitializerHacky");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);
};


void
add_ctor_initializer_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
