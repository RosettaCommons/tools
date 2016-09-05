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

#include "member_calls.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

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

RewriteClassMemberCalls::RewriteClassMemberCalls(
	tooling::Replacements *Replace,
	const char *tag,
	bool debug
) :
	ReplaceMatchCallback(Replace, tag, debug ) {}

void RewriteClassMemberCalls::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
	const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
	const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");

	if(!rewriteThisFile(expr, sm))
		return;

	// Get castFrom and castTo variable types
	const std::string castFromType(
		stripQualifiers(
			castFrom ? QualType::getAsString( castFrom->getType().split() ) : ""
		)
	);
	const std::string castToType(
		stripQualifiers(
			castTo ? QualType::getAsString( castTo->getType().split() ) : ""
		)
	);

	// Desugared types
	const std::string castFromTypeD(
		stripQualifiers(
			castFrom ? QualType::getAsString( castFrom->getType().getSplitDesugaredType() ) : ""
		)
	);
	const std::string castToTypeD(
		stripQualifiers(
			castTo ? QualType::getAsString( castTo->getType().getSplitDesugaredType() ) : ""
		)
	);

	const std::string origCode = getText(sm, expr);

	if(debug()) {
		const std::string locStr( expr->getSourceRange().getBegin().printToString(sm) );

		llvm::errs() << locStr << "\n";
		llvm::errs() << color("red") << origCode << color("") << "\n";

		llvm::errs() << "castFrom: " << color("purple") << castFromType << color("");
		if(castFromType != castFromTypeD)
			llvm::errs() << "          " << color("purple") << castFromTypeD << color("");
		llvm::errs() << "\n";

		llvm::errs() << "castTo:   " << color("brown") << castToType << color("");
		if(castToType != castToTypeD)
			llvm::errs() << "          " << color("brown") << castToTypeD << color("");
		llvm::errs() << "\n\n";
	}

	// Get original code
	if(origCode.empty())
		return;

	// Same thing, do nothing
	if(castFromTypeD == castToTypeD)
		return;
	if(castFromTypeD == extractContainedType(castToTypeD))
		return;

	// Determine cast type
	std::string type(castToType);

	// This cast to is a container, then cast to the contained type
	std::string contained_type = extractContainedType(castToType);
	if(!contained_type.empty())
		type = contained_type;
	else {
		// ... using full desugared definition
		contained_type = extractContainedType(castToTypeD);
		if(!contained_type.empty())
			type = contained_type;
	}

	// Full type definition not yet rewritten in original code, so do it here
	replace(type, "owning_ptr", "shared_ptr");
	replace(type, "access_ptr", "weak_ptr");

	std::string newCode;
	if(origCode == "0" || origCode == "NULL")
		newCode = type + "()";
	else
		newCode = type + "( " + origCode + " )";

	doRewrite(sm, expr, origCode, newCode);
}

void
add_member_calls_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements, bool dangerous_rewrites )
{
	using namespace clang::ast_matchers;

	if ( ! dangerous_rewrites ) return;
	/*
		std::vector<ClassAOP> as_vector_;
		void set_a_vector1(ClassA *a) {
			as_vector_.push_back(a);
		}

	  | `-CXXMemberCallExpr 0x2ef07d0 <col:3, col:25> 'void'
	  |   |-MemberExpr 0x2ef07a0 <col:3, col:14> '<bound member function type>' .push_back 0x2e616a0
	  |   | `-MemberExpr 0x2ef0670 <col:3> 'std::vector<ClassAOP>':'class std::vector<class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue ->as_vector_ 0x2e78910
	  |   |   `-CXXThisExpr 0x2ef0658 <col:3> 'class ClassB *' this
	  |   `-MaterializeTemporaryExpr 0x2ef08a8 <col:24> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' xvalue
	  |     `-CXXBindTemporaryExpr 0x2ef0888 <col:24> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x2ef0880)
	  |       `-CXXConstructExpr 0x2ef0848 <col:24> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
	  |         `-ImplicitCastExpr 0x2ef0830 <col:24> 'class ClassA *' <LValueToRValue>
	  |           `-DeclRefExpr 0x2ef0710 <col:24> 'class ClassA *' lvalue ParmVar 0x2ee9c70 'a' 'class ClassA *'
	*/


	finder.addMatcher(
		memberCallExpr(
			allOf(
				has(
					bindTemporaryExpr(
						has(
							constructExpr(
								has(
									declRefExpr().bind("castFrom")
								)
							)
						)
					).bind("expr")
				),
				has(
					memberExpr(
						has(
							memberExpr( containsUtilityPointer() ).bind("castTo")
						)
					)
				)
			)
		),
		new RewriteClassMemberCalls(replacements, "ClassMemberCalls"));


	/*
		std::vector<ClassAOP> as_vector_;
		void set_aref_vector1(ClassA & a) {
			as_vector_.push_back(&a);
		}

	  | `-CXXMemberCallExpr 0x3e22090 <col:3, col:26> 'void'
	  |   |-MemberExpr 0x3e22060 <col:3, col:14> '<bound member function type>' .push_back 0x3d7e850
	  |   | `-MemberExpr 0x3e21f18 <col:3> 'std::vector<ClassAOP>':'class std::vector<class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue ->as_vector_ 0x3d95ac0
	  |   |   `-CXXThisExpr 0x3e21f00 <col:3> 'class ClassB *' this
	  |   `-MaterializeTemporaryExpr 0x3e22158 <col:24, col:25> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' xvalue
	  |     `-CXXBindTemporaryExpr 0x3e22138 <col:24, col:25> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x3e22130)
	  |       `-CXXConstructExpr 0x3e220f0 <col:24, col:25> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
	  |         `-UnaryOperator 0x3e21fe0 <col:24, col:25> 'class ClassA *' prefix '&'
	  |           `-DeclRefExpr 0x3e21fb8 <col:25> 'class ClassA' lvalue ParmVar 0x3e071a0 'a' 'class ClassA &'
	*/

	finder.addMatcher(
		memberCallExpr(
			allOf(
				has(
					bindTemporaryExpr(
						has(
							constructExpr(
								has(
									unaryOperator(
										has(
											declRefExpr().bind("castFrom")
										)
									)
								)
							)
						)
					).bind("expr")
				),
				has(
					memberExpr(
						has(
							memberExpr( containsUtilityPointer() ).bind("castTo")
						)
					)
				)
			)
		),
		new RewriteClassMemberCalls(replacements, "ClassMemberCalls:unaryOperator"));

}
