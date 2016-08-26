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

#include "ctor_initializer.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

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

RewriteCtorInitializer::RewriteCtorInitializer(
	tooling::Replacements *Replace,
	const char *tag ) :
	ReplaceMatchCallback(Replace, tag) {}

void RewriteCtorInitializer::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const Expr *node = Result.Nodes.getStmtAs<Expr>("literal");

	if(!rewriteThisFile(node, sm))
		return;

	std::string origCode( getText(sm, node) );
	if(origCode != "0" && origCode != "NULL")
		return;

	std::string newCode = "/* " + origCode + " */";
	doRewrite(sm, node, origCode, newCode);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace NULL in constructor initializers

/*
 * NOTE: This is a separate callback because NULL is actually __null in stddef.h,
 * and that's in a file that we don't rewrite here.
 */

RewriteCtorInitializerHacky::RewriteCtorInitializerHacky(
		tooling::Replacements *Replace,
		const char *tag
) :
	ReplaceMatchCallback(Replace, tag)
{}

void RewriteCtorInitializerHacky::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const Expr *node = Result.Nodes.getStmtAs<Expr>("literal");
	const Expr *construct = Result.Nodes.getStmtAs<Expr>("construct");

	std::string matchedCode( getText(sm, node) );
	if(matchedCode != "__null")
		return;

	std::string origCode( getText(sm, construct) );
	std::string newCode = origCode;
	replace(newCode, "NULL", "/* NULL */");
	doRewrite(sm, construct, origCode, newCode);
}


void
add_ctor_initializer_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{

	using namespace clang::ast_matchers;
	/*
	CXXCtorInitializer Field 0x4f43710 'this_weak_ptr_' 'AtomTreeCAP':'class utility::pointer::access_ptr<const class core::kinematics::AtomTree>'
	|-CXXConstructExpr 0x6779c38 <line:62:2, col:18> 'AtomTreeCAP':'class utility::pointer::access_ptr<const class core::kinematics::AtomTree>' 'void (pointer)'
	| `-ImplicitCastExpr 0x6779c20 <col:17> 'pointer':'const class core::kinematics::AtomTree *' <NullToPointer>
	|   `-IntegerLiteral 0x6779b68 <col:17> 'int' 0
	*/
	
	finder.addMatcher(
		constructorDecl(
			forEachConstructorInitializer(
				ctorInitializer(
					withInitializer(
						constructExpr(
							has(
								integerLiteral().bind("literal")
							)
						).bind("construct")
					)
				)
			)
		),
		new RewriteCtorInitializer(replacements, "CtorInitializer:0"));



	/*
	|-CXXCtorInitializer Field 0x4ea45f0 'pep_prev_' 'NodeAP':'class utility::pointer::access_ptr<class core::environment::FoldTreeSketch::Node>'
	| |-CXXConstructExpr 0x5b0dcd8 <line:325:3, col:19> 'NodeAP':'class utility::pointer::access_ptr<class core::environment::FoldTreeSketch::Node>' 'void (pointer)'
	| | `-ImplicitCastExpr 0x5b0dcc0 </local/luki/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> 'pointer':'class core::environment::FoldTreeSketch::Node *' <NullToPointer>
	| |   `-GNUNullExpr 0x5b0dbe8 <col:18> 'long'
	*/

	finder.addMatcher(
		constructorDecl(
			forEachConstructorInitializer(
				ctorInitializer(
					withInitializer(
						constructExpr(
							has(
								expr(isNullExpr()).bind("literal")
							)
						).bind("construct")
					)
				)
			)
		),
		new RewriteCtorInitializerHacky(replacements, "CtorInitializer:nullExpr"));
}
