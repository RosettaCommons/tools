// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_pose_dynamic_cast_HH
#define INCLUDED_matchers_rewrite_pose_dynamic_cast_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

/*
	Replace dynamic casts in pose/conformation special cases:

	From:
		ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &( pose.conformation() ) );

	To:
		ProtectedConformationCOP conf = utility::pointer::dynamic_pointer_cast< ProtectedConformation const >( pose.conformation_ptr() );
*/

class RewritePoseDynamicCast : public ReplaceMatchCallback {
public:
	RewritePoseDynamicCast(
		clang::tooling::Replacements *Replace,
		const char *replacementCastCode,
		const char *tag = "PoseDynamicCast");

	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

private:
	std::string replacementCastCode;
};

void
add_pose_dynamic_cast_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif

