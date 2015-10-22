// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_matchers_base_HH
#define INCLUDED_matchers_base_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"
#include "clang/Tooling/Refactoring.h"

#include "utils.hh"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Base classes

class ReplaceMatchCallback : public clang::ast_matchers::MatchFinder::MatchCallback {
public:
	ReplaceMatchCallback(clang::tooling::Replacements *Replace, const char *tag)
		: Replace(Replace), tag(tag), debug_( false ) {}
	ReplaceMatchCallback(clang::tooling::Replacements *Replace, const char *tag, bool debug)
		: Replace(Replace), tag(tag), debug_( debug ) {}

private:
	clang::tooling::Replacements *Replace;

protected:
	std::string tag;
	bool debug_;

	template <typename T>
	void doRewrite(
		clang::SourceManager &sm,
		T * node,
		const std::string & origCode,
		const std::string & newCode
	) {
		if(origCode == newCode)
			return;

		if(!checkAndDumpRewrite(tag, sm, node, newCode))
			return;
		Replace->insert(Replacement(sm, node, newCode));
	}

	template <typename T>
	bool rewriteThisFile(T * node, clang::SourceManager & sm) {
		using namespace clang;
		if(!node)
			return false;
		const FullSourceLoc FullLocation = FullSourceLoc(node->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID()) {
			// llvm::errs() << tag << " skipping file: " << node->getSourceRange().getBegin().printToString(sm) << "\n";
			return false;
		}

		if (debug_)
			llvm::errs() << color("blue") << tag << " matched: " << node->getSourceRange().getBegin().printToString(sm) << color("") << "\n";

		return true;
	}

};

#endif
