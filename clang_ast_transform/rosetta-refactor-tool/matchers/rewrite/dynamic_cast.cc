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

#include "dynamic_cast.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

/*
	Replace dynamic casts with smart pointers, i.e.

	From:
		WrappedRealOP val = dynamic_cast< WrappedReal * > ( data()() );
		WrappedRealOP val = static_cast< WrappedReal * > ( data()() );
		CacheableStringFloatMapCOP data = dynamic_cast< CacheableStringFloatMap const * > ( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
	To:
		WrappedRealOP val = utility::pointer::dynamic_pointer_cast< WrappedReal > ( data() );
		WrappedRealOP val = utility::pointer::static_pointer_cast< WrappedReal > ( data() );
		CacheableStringFloatMapCOP data = utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap const > ( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
*/

RewriteDynamicCast::RewriteDynamicCast(
	tooling::Replacements *Replace,
	const char *replacementCastCode,
	const char *tag ) :
	ReplaceMatchCallback(Replace, tag),
	replacementCastCode(replacementCastCode) {}

void RewriteDynamicCast::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;

	const Expr *dyncastexpr = Result.Nodes.getStmtAs<Expr>("dyncastexpr");
	const Expr *operatorexpr = Result.Nodes.getStmtAs<Expr>("operatorexpr");
	const Expr *membercallexpr = Result.Nodes.getStmtAs<Expr>("membercallexpr");

	if(!rewriteThisFile(dyncastexpr, sm))
		return;

	// Get original code and cast type
	const std::string origCode = getText(sm, dyncastexpr);
	const std::string origCastType = QualType::getAsString( dyncastexpr->getType().split() );
	const std::string origCastCode = getText(sm, operatorexpr ? operatorexpr : membercallexpr);

	if(origCode.empty() || origCastCode.empty() || origCastType.empty())
		return;

	std::string castType = stripQualifiers(origCastType);
	if(endsWith(castType, "*"))
		castType = trim(std::string(castType, 0, castType.length() -1));
	if(origCastType.find("const ") != std::string::npos)
		castType += " const";

	std::string newCastCode = origCastCode;

	// origCastCode should end with () or .get, so strip that and hope for the best :)
	if(endsWith(origCastCode, ".get()"))
		newCastCode = std::string(origCastCode, 0, origCastCode.length() -6);
	else if(endsWith(origCastCode, "()"))
		newCastCode = std::string(origCastCode, 0, origCastCode.length() -2);

	std::string newCode =
		replacementCastCode + "< " + castType + " > "
			+ "( " + newCastCode + " )";

	doRewrite(sm, dyncastexpr, origCode, newCode);
}

void
add_dynamic_cast_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{

	using namespace clang::ast_matchers;
	/*
	 * Dynamic cast with call operator() to return naked pointer

	WrappedRealOP val = dynamic_cast< WrappedReal * > ( data()() );

	CXXDynamicCastExpr 0x7f5caef5af28 <col:22, col:63> 'WrappedReal *' dynamic_cast<WrappedReal *> <Dynamic>
	`-CXXOperatorCallExpr 0x7f5caef59de0 <col:54, col:61> 'pointer':'class utility::pointer::ReferenceCount *'
	  |-ImplicitCastExpr 0x7f5caef59dc8 <col:60, col:61> 'pointer (*)(void) const' <FunctionToPointerDecay>
	  | `-DeclRefExpr 0x7f5caef59d48 <col:60, col:61> 'pointer (void) const' lvalue CXXMethod 0x7f5caf1f97e0 'operator()' 'pointer (void) const'
	  `-ImplicitCastExpr 0x7f5caef59e20 <col:54, col:59> 'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' <NoOp>
	    `-CXXBindTemporaryExpr 0x7f5caef59d28 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' (CXXTemporary 0x7f5caef59d20)
	      `-CXXMemberCallExpr 0x7f5caef59cf0 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>'
	        `-MemberExpr 0x7f5caef59cc0 <col:54> '<bound member function type>' ->data 0x7f5caef8a990

	CXXDynamicCastExpr 0x54731c0 <col:33, col:82> 'const class numeric::expression_parser::FunctionToken *' dynamic_cast<const class numeric::expression_parser::FunctionToken *> <Dynamic>
	`-CXXOperatorCallExpr 0x5473158 <col:72, col:81> 'pointer':'const class numeric::expression_parser::Token *'
	  |-ImplicitCastExpr 0x5473140 <col:80, col:81> 'pointer (*)(void) const' <FunctionToPointerDecay>
	  | `-DeclRefExpr 0x5473118 <col:80, col:81> 'pointer (void) const' lvalue CXXMethod 0x5358c50 'operator()' 'pointer (void) const'
	  `-ImplicitCastExpr 0x5473198 <col:72> 'const class utility::pointer::owning_ptr<const class numeric::expression_parser::Token>' lvalue <NoOp>
	    `-DeclRefExpr 0x54730f0 <col:72> 'TokenCOP':'class utility::pointer::owning_ptr<const class numeric::expression_parser::Token>' lvalue Var 0x5472ea0 'toptoken' 'TokenCOP':'class utility::pointer::owning_ptr<const class numeric::expression_parser::Token>'
	*/

	finder.addMatcher(
		dynamicCastExpr(
			has(
				operatorCallExpr(
					allOf(
						// CHILD EXPR: operator() for owning_ptr::operator()
						has(
							declRefExpr( isCallOperator() )
						),
						// CHILD EXPR: castFrom
						anyOf(
							has(
								memberExpr( isUtilityPointer() )
							),
							has(
								declRefExpr( isUtilityPointer() )
							),
							has(
								bindTemporaryExpr( isUtilityPointer() )
							)
						)
					)
				).bind("operatorexpr")
			)
		).bind("dyncastexpr"),
		new RewriteDynamicCast(replacements, "utility::pointer::dynamic_pointer_cast", "DynamicCast:operatorCallExpr"));

	/*
	 * Same as above, but for static_cast< >, i.e. staticCastExpr
	 */

	finder.addMatcher(
		staticCastExpr(
			has(
				operatorCallExpr(
					allOf(
						// CHILD EXPR: operator() for owning_ptr::operator()
						has(
							declRefExpr( isCallOperator() )
						),
						// CHILD EXPR: castFrom
						anyOf(
							has(
								memberExpr( isUtilityPointer() )
							),
							has(
								declRefExpr( isUtilityPointer() )
							),
							has(
								bindTemporaryExpr( isUtilityPointer() )
							)
						)
					)
				).bind("operatorexpr")
			)
		).bind("dyncastexpr"),
		new RewriteDynamicCast(replacements, "utility::pointer::static_pointer_cast", "StaticCast"));

	/*
	 * Dynamic cast without call operator() to return naked pointer

	CacheableStringFloatMapCOP data = dynamic_cast< CacheableStringFloatMap const * > ( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );

	`-CXXDynamicCastExpr 0xa9f8160 <line:527:5, line:528:102> 'const class basic::datacache::CacheableStringFloatMap *' dynamic_cast<const class basic::datacache::CacheableStringFloatMap *> <Dynamic>
	  `-CXXMemberCallExpr 0xa9f80e8 <col:6, col:100> 'const class basic::datacache::CacheableData *'
	    |-MemberExpr 0xa9f80b8 <col:6, col:18> '<bound member function type>' .get_raw_const_ptr 0x98f5530
	    | `-ImplicitCastExpr 0xa9f8118 <col:6, col:16> 'const class basic::datacache::DataCache<class basic::datacache::CacheableData>' lvalue <UncheckedDerivedToBase (DataCache)>
	    |   `-CXXMemberCallExpr 0xa9f7f70 <col:6, col:16> 'const BasicDataCache':'const class basic::datacache::BasicDataCache' lvalue
	    |     `-MemberExpr 0xa9f7f40 <col:6, col:11> '<bound member function type>' .data 0x8a23230
	    |       `-DeclRefExpr 0xa9f7ea8 <col:6> 'const core::pose::Pose':'const class core::pose::Pose' lvalue ParmVar 0xa9f76a0 'pose' 'const core::pose::Pose &'
	    `-ImplicitCastExpr 0xa9f8138 <col:37, col:79> 'std::size_t':'unsigned long' <IntegralCast>
	      `-DeclRefExpr 0xa9f8070 <col:37, col:79> 'enum core::pose::datacache::CacheableDataType::Enum' EnumConstant 0x90dca30 'ARBITRARY_FLOAT_DATA' 'enum core::pose::datacache::CacheableDataType::Enum'
	*/

	finder.addMatcher(
		dynamicCastExpr(
			has(
				memberCallExpr().bind("membercallexpr")
			),
			hasParent(
				constructExpr( isUtilityPointer() )
			)
		).bind("dyncastexpr"),
		new RewriteDynamicCast(replacements, "utility::pointer::dynamic_pointer_cast", "DynamicCast:memberCallExpr"));

	/*
	 * Same as above, but for static_cast< >, i.e. staticCastExpr
	 */

	finder.addMatcher(
		staticCastExpr(
			has(
				memberCallExpr().bind("membercallexpr")
			),
			hasParent(
				constructExpr( isUtilityPointer() )
			)
		).bind("dyncastexpr"),
		new RewriteDynamicCast(replacements, "utility::pointer::static_pointer_cast", "StaticCast:memberCallExpr"));
}
