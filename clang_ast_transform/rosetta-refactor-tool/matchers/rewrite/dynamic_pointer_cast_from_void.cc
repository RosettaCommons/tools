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

#include "field_decl.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include <string>

using namespace clang;

/*
	Find all function calls in class methods
	To get a parsable list, run:
	./run.sh test-cases.cc -matchers=find_calls -colors=false|grep ' => '|sed 's/ => /\t/'
*/

DynamicPointerCastFromVoidFinder::DynamicPointerCastFromVoidFinder(
	tooling::Replacements *Replace,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag)
{}

void DynamicPointerCastFromVoidFinder::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	DeclRefExpr const * dpc = Result.Nodes.getStmtAs<DeclRefExpr>("dpc_node");
	//DeclRefExpr const * vptr = Result.Nodes.getStmtAs<DeclRefExpr>("ptr_to_void");

	if(!rewriteThisFile(dpc, sm))
		return;

	//std::cout << "Found a match!" << std::endl;
	//
	//ASTTemplateKWAndArgsInfo const * tmpl_args = dpc->getTemplateKWAndArgsInfo();
	//std::cout << " num template args: " << tmpl_args->NumTemplateArgs << std::endl;
	//TemplateArgumentLoc const & arg1 = (*tmpl_args)[0];
	//std::cout << " arg1: " << arg1.getTypeSourceInfo()->getType().getAsString() << std::endl;
	//
	//std::cout << " vptr: " << vptr->getFoundDecl()->getQualifiedNameAsString() << std::endl;

	std::string origCode( getText(sm,dpc) );
	doRewrite( sm, dpc, origCode, "dynamic_void_pointer_cast" );

}

void
add_dynamic_pointer_cast_from_void_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;

	//  |   |             | `-DeclRefExpr 0xb20d180 <col:41, col:109> 'shared_ptr<class core::conformation::Residue> (const shared_ptr<void> &)' lvalue Function 0xb20d080 'dynamic_pointer_cast' 'shared_ptr<class core::conformation::Residue> (const shared_ptr<void> &)' (FunctionTemplate 0x3f34df0 'dynamic_pointer_cast')
	//  |   |             `-ImplicitCastExpr 0xb20d2b0 <col:113> 'const shared_ptr<void>':'const class boost::shared_ptr<void>' lvalue <NoOp>
	//  |   |               `-DeclRefExpr 0xb20cca0 <col:113> 'utility::pointer::ReferenceCountOP':'class boost::shared_ptr<void>' lvalue Var 0xb20c720 'data_pointer' 'utility::pointer::ReferenceCountOP':'class boost::shared_ptr<void>'

	finder.addMatcher(
		callExpr(
			allOf(
				has(
					declRefExpr( isDynamicPointerCast() ).bind( "dpc_node" )
				),
				has(
					declRefExpr( isSharedPtrToVoid() ).bind( "ptr_to_void" )
				)
			)
		),
		new DynamicPointerCastFromVoidFinder( replacements, "declRefExpr" ));
}
