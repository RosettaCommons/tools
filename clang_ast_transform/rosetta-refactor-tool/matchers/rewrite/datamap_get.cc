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

#include "datamap_get.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

/*
	Replace:
		data.get<SomeClass*>(x, y)
	with:
		data.get_ptr<SomeClass>(x, y)

	Example:
	scorefxn_ = data.get<core::scoring::ScoreFunction*>("scorefxns", scorefxn_name);

  | | `-ExprWithCleanups 0xa6ddfd8 <line:132:3, col:81> 'class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction>' lvalue
  | |   `-CXXOperatorCallExpr 0xa6ddf90 <col:3, col:81> 'class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction>' lvalue
  | |     |-ImplicitCastExpr 0xa6ddf78 <col:13> 'class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction> &(*)(pointer)' <FunctionToPointerDecay>
  | |     | `-DeclRefExpr 0xa6ddef0 <col:13> 'class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction> &(pointer)' lvalue CXXMethod 0x469aa30 'operator=' 'class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction> &(pointer)'
  | |     |-MemberExpr 0xa6dd130 <col:3> 'core::scoring::ScoreFunctionCOP':'class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction>' lvalue ->scorefxn_ 0x46b3ef0
  | |     | `-CXXThisExpr 0xa6dd118 <col:3> 'class protocols::features::InterfaceFeatures *' this
  | |     `-ImplicitCastExpr 0xa6dded8 <col:15, col:81> 'pointer':'const class core::scoring::ScoreFunction *' <NoOp>
  | |       `-CXXMemberCallExpr 0xa6dd880 <col:15, col:81> 'class core::scoring::ScoreFunction *':'class core::scoring::ScoreFunction *'
  | |         |-MemberExpr 0xa6dd7f0 <col:15, col:53> '<bound member function type>' .get 0xa6dd6b0
  | |         | `-ImplicitCastExpr 0xa6dd8b8 <col:15> 'const class basic::datacache::DataMap' lvalue <NoOp>
  | |         |   `-DeclRefExpr 0xa6dd160 <col:15> 'basic::datacache::DataMap':'class basic::datacache::DataMap' lvalue ParmVar 0xa6da450 'data' 'basic::datacache::DataMap &'
  | |         |-CXXBindTemporaryExpr 0xa6dda98 <col:55> 'std::string':'class std::basic_string<char>' (CXXTemporary 0xa6dda90)
  | |         | `-CXXConstructExpr 0xa6dda58 <col:55> 'std::string':'class std::basic_string<char>' 'void (class std::basic_string<char> &&) noexcept' elidable
  | |         |   `-MaterializeTemporaryExpr 0xa6dda38 <col:55> 'class std::basic_string<char>' xvalue
  | |         |     `-CXXBindTemporaryExpr 0xa6dda18 <col:55> 'std::string':'class std::basic_string<char>' (CXXTemporary 0xa6dda10)
  | |         |       `-ImplicitCastExpr 0xa6dd9f8 <col:55> 'std::string':'class std::basic_string<char>' <ConstructorConversion>
  | |         |         `-CXXConstructExpr 0xa6dd9b8 <col:55> 'std::string':'class std::basic_string<char>' 'void (const char *, const class std::allocator<char> &)'
  | |         |           |-ImplicitCastExpr 0xa6dd8d0 <col:55> 'const char *' <ArrayToPointerDecay>
  | |         |           | `-StringLiteral 0xa6dd2f0 <col:55> 'const char [10]' lvalue "scorefxns"
  | |         |           `-CXXDefaultArgExpr 0xa6dd990 <<invalid sloc>> 'const class std::allocator<char>':'const class std::allocator<char>' lvalue
  | |         `-CXXBindTemporaryExpr 0xa6ddaf8 <col:68> 'std::string':'class std::basic_string<char>' (CXXTemporary 0xa6ddaf0)
  | |           `-CXXConstructExpr 0xa6ddab8 <col:68> 'std::string':'class std::basic_string<char>' 'void (const class std::basic_string<char> &)'
  | |             `-DeclRefExpr 0xa6dd328 <col:68> 'const string':'const class std::basic_string<char>' lvalue Var 0xa6dc8c0 'scorefxn_name' 'const string':'const class std::basic_string<char>'
*/

RewriteDataMapGet::RewriteDataMapGet(
	tooling::Replacements *Replace,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag)
{}

void RewriteDataMapGet::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;

	const Expr *declrefexpr = Result.Nodes.getStmtAs<Expr>("declrefexpr");
	// const Expr *membercallexpr = Result.Nodes.getStmtAs<Expr>("membercallexpr");
	const Expr *operatorcallexpr = Result.Nodes.getStmtAs<Expr>("operatorcallexpr");
	const MemberExpr *memberexpr = Result.Nodes.getStmtAs<MemberExpr>("memberexpr");

	if(!rewriteThisFile(operatorcallexpr, sm))
		return;

	// Get original code and cast type
	const std::string origCode = getText(sm, memberexpr);
	const std::string memberType = QualType::getAsString( declrefexpr->getType().split() );
	const std::string memberName = memberexpr->getMemberNameInfo().getName().getAsString();

	if(origCode.empty() || memberType != "basic::datacache::DataMap" || memberName != "get")
		return;

	size_t p = origCode.find("<");
	std::string newCode =
		trim(std::string(origCode, 0, p)) + "_ptr<" +
		trim(std::string(origCode, p), "<>* ") +
		">";

	doRewrite(sm, memberexpr, origCode, newCode);
}

void
add_datamap_get_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;

	finder.addMatcher(
		operatorCallExpr(
			isUtilityPointer(),
			has(
				memberCallExpr(
					has(
						memberExpr(
							has(
								declRefExpr().bind("declrefexpr")
							)
						).bind("memberexpr")
					)
				).bind("membercallexpr")
			)
		).bind("operatorcallexpr"),
		new RewriteDataMapGet(replacements));
}
