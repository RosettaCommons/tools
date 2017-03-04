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

#include "record_decl.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include <string>

using namespace clang;

RecordDeclFinder::RecordDeclFinder(tooling::Replacements *Replace) :
	ReplaceMatchCallback(Replace, "RecordDeclFinder")
	{}

// Main callback for all matches
void RecordDeclFinder::run(const ast_matchers::MatchFinder::MatchResult &Result) {

	SourceManager *sm = Result.SourceManager;
	const CXXRecordDecl *decl = Result.Nodes.getStmtAs<CXXRecordDecl>("recorddecl");
	const Decl *parent = Result.Nodes.getStmtAs<Decl>("parent");
	if(!rewriteThisFile(decl, *sm))
		return;
	if(!decl->isCompleteDefinition())
		return;

	const std::string name = decl->getQualifiedNameAsString();
	const std::string loc = decl->getSourceRange().getBegin().printToString(*sm);

	const CharSourceRange range = CharSourceRange::getTokenRange(decl->getSourceRange());
	SourceLocation SpellingBegin = sm->getSpellingLoc(range.getBegin());
	SourceLocation SpellingEnd = sm->getSpellingLoc(range.getEnd());
	std::pair<FileID, unsigned> Start = sm->getDecomposedLoc(SpellingBegin);
	std::pair<FileID, unsigned> End = sm->getDecomposedLoc(SpellingEnd);

	llvm::outs()
		<< decl->getKindName() << "\t"
		<< name << "\t"
		<< "" << "\t"
		<< loc << "\t"
		<< Start.second << "-" << End.second << "\t"
		<< (parent && parent->getKind() == Decl::ClassTemplate) << "\t"
		<< decl->getAccess() << "\t"
		<< decl->isPolymorphic() << "\t"
		;
	llvm::outs() << "\n";
	if ( decl->getNumBases() > 0 ) {
	  for ( CXXRecordDecl::base_class_const_iterator iter = decl->bases_begin(), iter_end = decl->bases_end(); iter != iter_end; ++iter ) {
	    llvm::outs()
	      << "parent" << "\t"
	      << name << "\t"
	      << "" << "\t"
	      << loc << "\t"
	      << sm->getDecomposedLoc( iter->getLocStart() ).second << "-" << sm->getDecomposedLoc( iter->getLocEnd() ).second << "\t"
	      //<< iter->getType().getAsString() << "\t"
	      //<< iter->getType().getUnqualifiedType().getAsString() << "\t"
	      << iter->getType().getCanonicalType().getAsString() << "\t"
	      << iter->getAccessSpecifier() << "\t"
	      << iter->isVirtual();
	    llvm::outs() << "\n";
	  }
	}
}


void
add_record_decl_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;
	finder.addMatcher(
		recordDecl(
			hasParent(
				decl().bind("parent")
			)
		).bind("recorddecl"),
		new RecordDeclFinder(replacements));
}
