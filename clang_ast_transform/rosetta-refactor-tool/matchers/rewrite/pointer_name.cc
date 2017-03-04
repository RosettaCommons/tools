// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/*
	Find instances where get_self_ptr() or get_self_weak_ptr() is used in a c'tor.
	This is illegal because the weak self-pointer isn't set yet, and will
	result in bad_weak_ptr exception at runtime.
*/

#include "pointer_name.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

/*
	Replace owning_ptr/access_ptr name in templates, etc.
*/

RewritePointerName::RewritePointerName(
	tooling::Replacements *Replace,
	const char *tag ) :
	ReplaceMatchCallback(Replace, tag) {}

void RewritePointerName::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const Stmt *node = Result.Nodes.getStmtAs<Stmt>("stmt");

	if(!rewriteThisFile(node, sm))
		return;
	
	std::string origCode( getText(sm, node) );
	std::string newCode( origCode );

	replace(newCode, "owning_ptr", "shared_ptr");
	replace(newCode, "access_ptr", "weak_ptr");

	doRewrite(sm, node, origCode, newCode);
}

void
add_pointer_name_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;
	// CXXUnresolvedConstructExpr in templates
	finder.addMatcher(
		unresolvedConstructExpr( isUtilityPointer() ).bind("stmt"),
		new RewritePointerName(replacements));
}
