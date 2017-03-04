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

#include "self_ptr_in_ctor.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include <string>

using namespace clang;

SelfPtrInCtorFinder::SelfPtrInCtorFinder(
	tooling::Replacements *Replace,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag) {}

void SelfPtrInCtorFinder::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const CXXConstructorDecl *ctor = Result.Nodes.getStmtAs<CXXConstructorDecl>("ctor");

	if(!rewriteThisFile(ctor, sm))
		return;

	if(!ctor->hasBody())
		return;

	// TODO: retrieve enclosing namespace of caller and call node -- how?
	// DeclRefExpr::getFoundDecl()->getQualifiedNameAsString()

	const std::string locStr( ctor->getSourceRange().getBegin().printToString(sm) );
	const std::string code = getText(sm, ctor);

	if(
		code.find("get_self_ptr()") == std::string::npos &&
		code.find("get_self_weak_ptr()") == std::string::npos
	)
		return;

	llvm::outs()
		<< "@ " << locStr << /* color("cyan") << " (" << tag << ")" << color("") << */ "\n"
		;
}

void
add_self_ptr_in_ctor_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;
	finder.addMatcher(
		constructorDecl(
			// Doesn't work... so ugly "solution" above
			//forEachDescendant(
			//	memberCallExpr(
			//		hasName("get_self_ptr")
			//	)
			//)
		).bind("ctor"),
		new SelfPtrInCtorFinder(replacements));
}
