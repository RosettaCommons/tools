// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_replace_member_calls_HH
#define INCLUDED_matchers_replace_member_calls_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace implicit casts in class member calls:
	This may not be safe to do!

	std::vector<ClassAOP> as_vector_;
		void set_a_vector1(ClassA *a) {
			as_vector_.push_back(a);
		}
		void set_aref_vector1(ClassA & a) {
			as_vector_.push_back(&a);
		}
*/

class RewriteClassMemberCalls : public ReplaceMatchCallback {
public:
	RewriteClassMemberCalls(
		clang::tooling::Replacements *Replace,
		const char *tag = "ClassMemberCalls",
		bool debug = false
	);


	void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

};


void
add_member_calls_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements, bool dangerous_rewrites );

#endif

