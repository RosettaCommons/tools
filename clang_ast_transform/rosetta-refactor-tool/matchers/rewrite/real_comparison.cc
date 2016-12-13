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

#include "real_comparison.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

/*
	Replace comparison of Reals (doubles or floats):
	double a, b;
	if(a == b) {}
	if(a != b) {}
*/

RewriteRealComparisons::RewriteRealComparisons(
	tooling::Replacements *Replace,
	const char *tag,
	const char *fn) :
	ReplaceMatchCallback(Replace, tag), fn_(fn) {}

void RewriteRealComparisons::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const Expr *a = Result.Nodes.getStmtAs<Expr>("a");
	const Expr *b = Result.Nodes.getStmtAs<Expr>("b");
	const BinaryOperator *op = Result.Nodes.getStmtAs<BinaryOperator>("op");

	if(!rewriteThisFile(op, sm))
		return;

	const std::string origCode = getText(sm, op);
	const std::string aCode = getText(sm, a);
	const std::string bCode = getText(sm, b);

	if(origCode.empty())
		return;

	std::string newCode = fn_ + "( " + aCode + ", " + bCode + " )";
	doRewrite(sm, op, origCode, newCode);
}

void
add_real_comparison_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;

	// a == b
	finder.addMatcher(
		binaryOperator(
			hasOperatorName("=="),
			hasLHS(ignoringParenImpCasts(
				anyOf(
					declRefExpr(to(varDecl(hasType(isReal())))).bind("a"),
					//declRefExpr(to(varDecl(hasType(isInteger())))).bind("a"),
					floatLiteral().bind("a"),
					integerLiteral().bind("a")
				)
			)),
			hasRHS(ignoringParenImpCasts(
				anyOf(
					declRefExpr(to(varDecl(hasType(isReal())))).bind("b"),
					declRefExpr(to(varDecl(hasType(isInteger())))).bind("b"),
					floatLiteral().bind("b"),
					integerLiteral().bind("b")
				)
			))
		).bind("op"),
		new RewriteRealComparisons(replacements, "RewriteRealComparisons:==", "numeric::equal_by_epsilon"));

	// a != b
	finder.addMatcher(
		binaryOperator(
			hasOperatorName("!="),
			hasLHS(ignoringParenImpCasts(
				anyOf(
					declRefExpr(to(varDecl(hasType(isReal())))).bind("a"),
					floatLiteral().bind("a"),
					integerLiteral().bind("a")
				)
			)),
			hasRHS(ignoringParenImpCasts(
				anyOf(
					declRefExpr(to(varDecl(hasType(isReal())))).bind("b"),
					floatLiteral().bind("b"),
					integerLiteral().bind("b")
				)
			))
		).bind("op"),
		new RewriteRealComparisons(replacements, "RewriteRealComparisons:!=", "!numeric::equal_by_epsilon"));
}
