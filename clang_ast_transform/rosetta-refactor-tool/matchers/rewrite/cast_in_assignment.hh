// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_cast_in_assignment_HH
#define INCLUDED_matchers_rewrite_cast_in_assignment_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace implicit casts of named pointers into owning/access pointers:

	Explicit cast needed to indicate that the shared_ptr is in charge
	of controlling the lifetime of the object; probably not a good
	idea to automagically rewrite these casts after all!

		ClassAOP a_;
		ClassA *a;
		a_ = a;

		utility::vector1<ClassAOP> as_;
		void set_a_vector1(ClassA *a) {
			as_[0] = a;
		}

		ClassAOP a_;
		ClassA & a;
		a_ = &a;

		utility::vector1<ClassAOP> as_;
		void set_aref_vector1(ClassA & a) {
			as_[0] = &a;
		}

		std::map<std::string,ClassAOP> as_map_;
		as_map_["some"] = a;
*/

////////////////////////////////////////////////////////////////////////////////////////////////////

class RewriteAssignmentsOper : public ReplaceMatchCallback {
public:
	RewriteAssignmentsOper(
		clang::tooling::Replacements *Replace,
		const char *tag = "AssignmentsOper",
		bool debug = false
	);

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

};


void
add_cast_in_assignment_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements, bool dangerous_rewrites );

#endif
