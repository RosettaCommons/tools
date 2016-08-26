// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/*
	Find constructors
*/

#include "constructor_decl.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include <string>

using namespace clang;

ConstructorDeclFinder::ConstructorDeclFinder(tooling::Replacements *Replace) :
	ReplaceMatchCallback(Replace, "ConstructorDeclFinder")
	{}

// Main callback for all matches
void ConstructorDeclFinder::run(const ast_matchers::MatchFinder::MatchResult &Result) {

	SourceManager *sm = Result.SourceManager;
	const CXXConstructorDecl *decl = Result.Nodes.getStmtAs<CXXConstructorDecl>("constructordecl");
	if(!rewriteThisFile(decl, *sm))
		return;

	const std::string name = decl->getQualifiedNameAsString();
	const std::string cls = decl->getParent()->getQualifiedNameAsString();
	const std::string loc = decl->getSourceRange().getBegin().printToString(*sm);

	const CharSourceRange range = CharSourceRange::getTokenRange(decl->getSourceRange());
	SourceLocation SpellingBegin = sm->getSpellingLoc(range.getBegin());
	SourceLocation SpellingEnd = sm->getSpellingLoc(range.getEnd());
	std::pair<FileID, unsigned> Start = sm->getDecomposedLoc(SpellingBegin);
	std::pair<FileID, unsigned> End = sm->getDecomposedLoc(SpellingEnd);

	llvm::outs()
		<< "constructor" << "\t"
		<< name << "\t"
		<< cls << "\t"
		<< loc << "\t"
		<< Start.second << "-" << End.second << "\t"
		<< decl->getAccess() << "\t"
		<< decl->isDefaultConstructor() << "\t"
		<< decl->isCopyConstructor() << "\t"
		<< decl->getNumCtorInitializers() << "\t"
		<< decl->getMinRequiredArguments() << "\t"
		;
	llvm::outs() << "\n";
}


void
add_constructor_decl_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang;
	using namespace clang::ast_matchers;

	finder.addMatcher(
		constructorDecl().bind("constructordecl"),
		new ConstructorDeclFinder(replacements));

}
