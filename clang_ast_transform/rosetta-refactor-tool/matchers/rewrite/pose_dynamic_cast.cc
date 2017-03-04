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

#include "pose_dynamic_cast.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

/*
	Replace dynamic casts in pose/conformation special cases:

	From:
		ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &( pose.conformation() ) );

	To:
		ProtectedConformationCOP conf = utility::pointer::dynamic_pointer_cast< ProtectedConformation const >( pose.conformation_ptr() );
*/

RewritePoseDynamicCast::RewritePoseDynamicCast(
	tooling::Replacements *Replace,
	const char *replacementCastCode,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag),
	replacementCastCode(replacementCastCode) {}

void RewritePoseDynamicCast::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;

	const Expr *dyncastexpr = Result.Nodes.getStmtAs<Expr>("dyncastexpr");
	const Expr *declrefexpr = Result.Nodes.getStmtAs<Expr>("declrefexpr");
	const Stmt *unaryoperator = Result.Nodes.getStmtAs<Stmt>("unaryoperator");
	const Expr *membercallexpr = Result.Nodes.getStmtAs<Expr>("membercallexpr");

	if(!rewriteThisFile(dyncastexpr, sm))
		return;

	// Get original code and cast type
	const std::string origCode = getText(sm, dyncastexpr);
	const std::string origCastType = QualType::getAsString( dyncastexpr->getType().split() );
	const std::string declRefType = QualType::getAsString( declrefexpr->getType().split() );
	const std::string origCastCode = getText(sm, unaryoperator);
	const std::string origCastCodeInside = getText(sm, membercallexpr);

	if(origCode.empty() || origCastCode.empty() || origCastType.empty())
		return;

	// We only are about Pose.conformation() here
	if(!endsWith(declRefType, "core::pose::Pose") || !endsWith(origCastCodeInside, ".conformation"))
		return;

	std::string castType = stripQualifiers(origCastType);
	if(endsWith(castType, "*"))
		castType = trim(std::string(castType, 0, castType.length() -1));
	if(origCastType.find("const ") != std::string::npos)
		castType += " const";

	std::string newCastCode = origCastCodeInside;
	replace(newCastCode, ".conformation", ".conformation_ptr");

	std::string newCode =
		replacementCastCode + "< " + castType + " > "
			+ "( " + newCastCode + "() )";

	doRewrite(sm, dyncastexpr, origCode, newCode);
}

void
add_pose_dynamic_cast_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;
	/*
	  ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &( pose.conformation() ) );

	  |-DeclStmt 0xa2675a0 <line:61:3, col:107>
	  | `-VarDecl 0xa2643b0 <col:3, col:106> col:28 conf 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>'
	  |   `-ExprWithCleanups 0xa267588 <col:28, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>'
	  |     `-CXXConstructExpr 0xa267550 <col:28, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' 'void (const class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation> &)' elidable
	  |       `-MaterializeTemporaryExpr 0xa267530 <col:35, col:106> 'const class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' lvalue
	  |         `-ImplicitCastExpr 0xa267518 <col:35, col:106> 'const class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' <NoOp>
	  |           `-CXXBindTemporaryExpr 0xa2671f8 <col:35, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' (CXXTemporary 0xa2671f0)
	  |             `-ImplicitCastExpr 0xa2671d8 <col:35, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' <ConstructorConversion>
	  |               `-CXXConstructExpr 0xa2671a0 <col:35, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' 'void (pointer)'
	  |                 `-CXXDynamicCastExpr 0xa264578 <col:35, col:106> 'const class protocols::environment::ProtectedConformation *' dynamic_cast<const class protocols::environment::ProtectedConformation *> <Dynamic>
	  |                   `-UnaryOperator 0xa264548 <col:81, col:104> 'const Conformation *' prefix '&'
	  |                     `-ParenExpr 0xa2644f8 <col:82, col:104> 'const Conformation':'const class core::conformation::Conformation' lvalue
	  |                       `-CXXMemberCallExpr 0xa2644d0 <col:84, col:102> 'const Conformation':'const class core::conformation::Conformation' lvalue
	  |                         `-MemberExpr 0xa2644a0 <col:84, col:89> '<bound member function type>' .conformation 0x97efd70
	  |                           `-DeclRefExpr 0xa264408 <col:84> 'const core::pose::Pose':'const class core::pose::Pose' lvalue ParmVar 0xa263b80 'pose' 'const core::pose::Pose &'
	*/

	finder.addMatcher(
		dynamicCastExpr(
			hasParent(
				constructExpr( isUtilityPointer() ).bind("constructexpr")
			),
			has(
				unaryOperator(
					has(
						memberCallExpr(
							has(
								memberExpr(
									has(
										declRefExpr().bind("declrefexpr")
									)
								).bind("membercallexpr")
							)
						)
					)
				).bind("unaryoperator")
			)
		).bind("dyncastexpr"),
		new RewritePoseDynamicCast(replacements, "utility::pointer::dynamic_pointer_cast", "DynamicCast:constructExpr"));
}
