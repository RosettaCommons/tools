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

#include "call_operator.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

/*
	Replace calls to operator() on owning_ptrs which return a named pointer to the object

		SomeOP a;
		foo( a() );
		foo( *a() );

	TODO:
	Replace with a.get() on OPs and a.lock().get() instead? Do we want this?
*/

RewriteCallOperator::RewriteCallOperator(
	tooling::Replacements *Replace,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag)
{}

RewriteCallOperator::RewriteCallOperator(
	tooling::Replacements *Replace,
	bool debug,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag, debug)
{}

void RewriteCallOperator::run(const ast_matchers::MatchFinder::MatchResult &Result)
{
	SourceManager &sm = *Result.SourceManager;

	const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
	const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
	const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");

	if(!rewriteThisFile(expr, sm))
		return;

	// Get original code
	const std::string origCode = getText(sm, expr);
	if(origCode.empty())
		return;

	if(origCode.find('*') != std::string::npos) {
		// Hack: try to avoid *e1() etc.
		return;
	}

	const std::string castFromType(
		castFrom ? QualType::getAsString( castFrom->getType().getSplitDesugaredType() ) : ""
	);
	const std::string castToType(
		castTo ? QualType::getAsString( castTo  ->getType().getSplitDesugaredType() ) : ""
	);

	if(debug()) {
		const std::string locStr( expr->getSourceRange().getBegin().printToString(sm) );

		llvm::errs() << locStr << "\n";
		llvm::errs() << color("red") << origCode << color("") << "\n";
		llvm::errs() << "castFrom: " << color("purple") << castFromType << color("");
		llvm::errs() << "\n";
		llvm::errs() << "castTo:   " << color("brown") << castToType << color("");
		llvm::errs() << "\n\n";
	}

	// origCode should end with (), so strip that call operator
	std::string newCode = std::string(origCode, 0, origCode.length() -2);
	doRewrite(sm, expr, origCode, newCode);
}

void
add_call_operator_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{

	using namespace clang::ast_matchers;
// Round 1
// (anyOf handles max 5 statements)

// The below matcher set covers these AST's

/*
CXXConstructExpr 0x7fcd48e1a460 <col:10, col:26> 'AtomOP':'class utility::pointer::owning_ptr<class core::kinematics::tree::Atom>' 'void (pointer)'
`-CXXOperatorCallExpr 0x7fcd48e1a3a0 <col:18, col:26> 'pointer':'class core::kinematics::tree::Atom *'
  |-ImplicitCastExpr 0x7fcd48e1a388 <col:25, col:26> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x7fcd48e1a360 <col:25, col:26> 'pointer (void) const' lvalue CXXMethod 0x7fcd48defcc0 'operator()' 'pointer (void) const'
  `-ImplicitCastExpr 0x7fcd48e1a3e0 <col:18> 'const class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' lvalue <NoOp>
    `-MemberExpr 0x7fcd48e1a330 <col:18> 'AtomAP':'class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' lvalue ->parent_ 0x7fcd48df0380

CXXConstructExpr 0x7fcd48e19a58 <col:10, col:27> 'AtomCOP':'class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' 'void (pointer)'
`-ImplicitCastExpr 0x7fcd48e19a40 <col:19, col:27> 'pointer':'const class core::kinematics::tree::Atom *' <NoOp>
  `-CXXOperatorCallExpr 0x7fcd48e19940 <col:19, col:27> 'pointer':'class core::kinematics::tree::Atom *'
    |-ImplicitCastExpr 0x7fcd48e19928 <col:26, col:27> 'pointer (*)(void) const' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x7fcd48e198a0 <col:26, col:27> 'pointer (void) const' lvalue CXXMethod 0x7fcd48defcc0 'operator()' 'pointer (void) const'
    `-MemberExpr 0x7fcd48e19870 <col:19> 'const AtomAP':'const class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' lvalue ->parent_ 0x7fcd48df0380
*/

finder.addMatcher(
	operatorCallExpr(
		allOf(
			// CHILD EXPR: operator()
			has(
				declRefExpr( isCallOperator() )
			),
			// CHILD EXPR: castFrom
			anyOf(
				has(
					memberExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					bindTemporaryExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					declRefExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					operatorCallExpr( isUtilityPointer() ).bind("castFrom")
				)
			),
			// PARENT Expr/Stmt: castTo
			// could have used hasAncestor() here to cover both cases
			// (with and without implicit cast) but this is more strict/specific
			anyOf(
				hasParent(
					constructExpr( isUtilityPointer() ).bind("castTo")
				),
				hasParent(
					implicitCastExpr(
						hasParent(
							constructExpr( isUtilityPointer() ).bind("castTo")
						)
					)
				)
			)
		)
	).bind("expr"),
	new RewriteCallOperator(replacements, "CallOperator:constructExpr"));


// Round 2

/*
ConditionalOperator 0x48c49e0 <col:8, col:31> 'const char *'
`-ImplicitCastExpr 0x48c4998 <col:8, col:18> '_Bool' <PointerToBoolean>
  `-CXXOperatorCallExpr 0x48c48d0 <col:8, col:18> 'pointer':'class utility::pointer::
    |-ImplicitCastExpr 0x48c48b8 <col:17, col:18> 'pointer (*)(void) const' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x48c4838 <col:17, col:18> 'pointer (void) const' lvalue CXXMethod 0x48990f0 'operator()' 'pointer (void) const'
    `-MemberExpr 0x48c4808 <col:8, col:11> 'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>':'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' lvalue ->second 0x489db20
      `-CXXOperatorCallExpr 0x48c47c8 <col:8> 'pointer':'const struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > *'
        |-ImplicitCastExpr 0x48c47b0 <col:9> 'pointer (*)(void) const' <FunctionToPointerDecay>
        | `-DeclRefExpr 0x48c4788 <col:9> 'pointer (void) const' lvalue CXXMethod 0x489afc0 'operator->' 'pointer (void) const'
        `-ImplicitCastExpr 0x48c4770 <col:8> 'const struct std::_Rb_tree_const_iterator<struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > >' lvalue <NoOp>
          `-DeclRefExpr 0x48c4748 <col:8> 'class ResourceManager::ResourcesMap::const_iterator':'struct std::_Rb_tree_const_iterator<struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > >' lvalue Var 0x48c0dc0 'r' 'class ResourceManager::ResourcesMap::const_iterator':'struct std::_Rb_tree_const_iterator<struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > >'

BinaryOperator 0x1ee6ec0 <col:12, /data/rosetta/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> '_Bool' '=='
`-CXXOperatorCallExpr 0x1ee6e50 </data/rosetta/main/source/src/utility/signals/Link.hh:108:12, col:18> 'pointer':'struct utility::signals::LinkUnit *'
  |-ImplicitCastExpr 0x1ee6e38 <col:17, col:18> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x1ee6e10 <col:17, col:18> 'pointer (void) const' lvalue CXXMethod 0x1ee38d0 'operator()' 'pointer (void) const'
  `-MemberExpr 0x1ee6de0 <col:12> 'const LinkUnitOP':'const class utility::pointer::owning_ptr<struct utility::signals::LinkUnit>' lvalue ->unit_ 0x1ee5be0
    `-CXXThisExpr 0x1ee6dc8 <col:12> 'const class utility::signals::Link *' this
*/

finder.addMatcher(
	operatorCallExpr(
		allOf(
			// CHILD EXPR: operator()
			has(
				declRefExpr( isCallOperator() )
			),
			// CHILD EXPR: castFrom
			anyOf(
				has(
					memberExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					declRefExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					bindTemporaryExpr( isUtilityPointer() ).bind("castFrom")
				)
			),
			// PARENT Expr/Stmt: castTo
			anyOf(
				hasParent(
					implicitCastExpr(
						hasParent(
							conditionalOperator()
						)
					).bind("castTo")
				),
				hasParent(
					binaryOperator().bind("castTo")
				),
				hasParent(
					unaryOperator().bind("castTo")
				)
			)
		)
	).bind("expr"),
	new RewriteCallOperator(replacements, "CallOperator:misc"));
}
