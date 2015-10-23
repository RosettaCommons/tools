// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

FieldDeclFinder::FieldDeclFinder(tooling::Replacements *Replace) :
	ReplaceMatchCallback(Replace, "FieldDeclFinder")
	{}

// Main callback for all matches
void FieldDeclFinder::run(const ast_matchers::MatchFinder::MatchResult &Result) {

	SourceManager *sm = Result.SourceManager;
	const FieldDecl *decl = Result.Nodes.getStmtAs<FieldDecl>("fielddecl");
	if(!rewriteThisFile(decl, *sm))
		return;

	const std::string type(
		QualType::getAsString( decl->getType().split() )
	);
	const std::string typeD(
		QualType::getAsString( decl->getType().getSplitDesugaredType() )
	);

	const std::string name = decl->getQualifiedNameAsString();
	const std::string cls = decl->getParent()->getQualifiedNameAsString();
	const std::string loc = decl->getSourceRange().getBegin().printToString(*sm);

	const CharSourceRange range = CharSourceRange::getTokenRange(decl->getSourceRange());
	SourceLocation SpellingBegin = sm->getSpellingLoc(range.getBegin());
	SourceLocation SpellingEnd = sm->getSpellingLoc(range.getEnd());
	std::pair<FileID, unsigned> Start = sm->getDecomposedLoc(SpellingBegin);
	std::pair<FileID, unsigned> End = sm->getDecomposedLoc(SpellingEnd);

	llvm::outs()
		<< "field" << "\t"
		<< name << "\t"
		<< cls << "\t"
		<< loc << "\t"
		<< Start.second << "-" << End.second << "\t"
		<< type << "\t"
		<< typeD << "\t"
		;
	llvm::outs() << "\n";
}

void
add_field_decl_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{

	using namespace clang::ast_matchers;

	FieldDeclFinder *cb = new FieldDeclFinder(replacements);

	finder.addMatcher(
		fieldDecl().bind("fielddecl"),
		cb);
}
