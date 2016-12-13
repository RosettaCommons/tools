// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_matchers_rewrite_add_serialization_code_HH
#define INCLUDED_matchers_rewrite_add_serialization_code_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"
#include "clang/AST/DeclCXX.h"


#include "../../original_tool/RosettaRefactorTool.hh"
#include "../../matchers_base.hh"

#include <string>

/*
	Auto add Cereal serialization code to a classes and structs.
	Note: this will fail if with overlapping rewrites, i.e.
	a class definition containing another (private) class definition.
*/

class AddSerializationCode : public ReplaceMatchCallback {

	typedef std::vector< std::string > Strings;
	typedef std::shared_ptr< Strings > StringsOP;
	typedef std::map< std::string, StringsOP > ClassFields;
	typedef std::pair< const clang::CXXRecordDecl *, const clang::Decl * > Decl_Parent_Pair;

	ClassFields class_fields;
	std::vector< Decl_Parent_Pair > class_recorddecls;
	clang::SourceManager * sm;

	RosettaRefactorTool * refactor_tool;

public:
	AddSerializationCode(
		clang::tooling::Replacements *Replace,
		RosettaRefactorTool *t);
	virtual ~AddSerializationCode();

	// Main callback for all matches
	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

	// Collect class field declarations
	void handle_field_decl(const clang::ast_matchers::MatchFinder::MatchResult &Result);

	// Collect class record declaration
	void handle_record_decl(const clang::ast_matchers::MatchFinder::MatchResult &Result);

	// Do the actual rewrite at the end of the translation unit
	virtual void onEndOfTranslationUnit();

};

void
add_serialization_code_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements, RosettaRefactorTool * rrt );

#endif
