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

#include "typedef.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;

/*
	Replace owning_ptr/access_ptr in declarations

	Type definitions:
		typedef utility::pointer::access_ptr< Some > SomeAP
		typedef utility::pointer::owning_ptr< Some > SomeOP
		typedef utility::pointer::access_ptr< Some const > SomeCAP
		typedef utility::pointer::owning_ptr< Some const > SomeCOP

	Parameter Declarations:
		void foo( utility::vector1< utility::pointer::owning_ptr<class ClassA> > & ops_ );
		void foo( utility::vector1< utility::pointer::access_ptr<class ClassA> > & aps_ );

	Field Declarations (in class header):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > ops_;
		utility::vector1< utility::pointer::access_ptr<class ClassA> > aps_;

	Function Declarations (return type):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > foo();
		utility::vector1< utility::pointer::access_ptr<class ClassA> > foo();

	Method Declarations (return type):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > foo();
		utility::vector1< utility::pointer::access_ptr<class ClassA> > foo();

	Variable Declarations (in method):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > ops_;
		utility::vector1< utility::pointer::access_ptr<class ClassA> > aps_;
*/

RewriteDecl::RewriteDecl(
	tooling::Replacements *Replace,
	const char *tag,
	const char *delim
) :
	ReplaceMatchCallback(Replace, tag),
	delim(delim)
{}

void RewriteDecl::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const Decl *node = Result.Nodes.getStmtAs<Decl>("decl");

	if(!rewriteThisFile(node, sm))
		return;

	const std::string origCode( getText(sm, node) );
	std::string newCode( origCode );
	std::string suffix;

	if(delim && *delim) {
		size_t p = newCode.find(delim);
		if(p != std::string::npos) {
			suffix = newCode.substr(p);
			newCode = newCode.substr(0, p);
		}
	}

	replace(newCode, "owning_ptr", "shared_ptr");
	replace(newCode, "access_ptr", "weak_ptr");

	doRewrite(sm, node, origCode, newCode + suffix);
}

void
add_typedef_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;
	// Typedefs for access_ptr and owning_ptr
	finder.addMatcher(
		decl( isTypedefDecl() ).bind("decl"),
		new RewriteDecl(replacements, "Decl:isTypedefDecl"));

	// Parameter declaration in methods
	finder.addMatcher(
		parmVarDecl().bind("decl"),
		new RewriteDecl(replacements, "Decl:paramVarDecl"));

	// Field (variable) declaration in classes
	finder.addMatcher(
		fieldDecl().bind("decl"),
		new RewriteDecl(replacements, "Decl:fieldDecl"));

	// Method return type declaration
	finder.addMatcher(
		methodDecl().bind("decl"),
		new RewriteDecl(replacements, "Decl:methodDecl", "{"));

	// Function return type declaration
	finder.addMatcher(
		functionDecl().bind("decl"),
		new RewriteDecl(replacements, "Decl:functionDecl", "{"));

	// Local method varialble declaration
	finder.addMatcher(
		varDecl().bind("decl"),
		new RewriteDecl(replacements, "Decl:varDecl"));
}
