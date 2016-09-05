// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
	ReplaceMatchCallback(clang::tooling::Replacements *Replace, const char *tag);
	ReplaceMatchCallback(clang::tooling::Replacements *Replace, const char *tag, bool debug);
	virtual ~ReplaceMatchCallback();

	bool debug() const { return debug_; }
	clang::tooling::Replacements * Replace() { return replace_; }
	std::string const & tag() const { return tag_; }

private:
	clang::tooling::Replacements * replace_;
	std::string tag_;
	bool debug_;

protected:

	template <typename T>
	void doRewrite(
		clang::SourceManager &sm,
		T * node,
		const std::string & origCode,
		const std::string & newCode
	) {
		if(origCode == newCode)
			return;

		if(!checkAndDumpRewrite(tag_, sm, node, newCode))
			return;
		replace_->insert(clang::tooling::Replacement(sm, node, newCode));
	}

	template <typename T>
	bool rewriteThisFile(T * node, clang::SourceManager & sm) {
 		using namespace clang;
		if(!node)
			return false;
		const FullSourceLoc FullLocation = FullSourceLoc(node->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID()) {
			// llvm::errs() << tag_ << " skipping file: " << node->getSourceRange().getBegin().printToString(sm) << "\n";
			return false;
		}

		if (debug_)
			llvm::errs() << color("blue") << tag_ << " matched: " << node->getSourceRange().getBegin().printToString(sm) << color("") << "\n";

		return true;
	}

};

#endif
