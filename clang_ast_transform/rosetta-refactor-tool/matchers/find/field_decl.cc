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
#include <iostream>

using namespace clang;

// Default behavior for the class is to report to the screen the fields declared in a single
// target file.  With the other constructor, it's possible to indicate that all files that are
// #included from the target file should be examined, and it's possible to say that the fields
// should not be written out to the screen.
FieldDeclFinder::FieldDeclFinder(tooling::Replacements *Replace) :
	ReplaceMatchCallback(Replace, "FieldDeclFinder"),
	verbose_( true ),
	match_target_file_only_( true )
{}

FieldDeclFinder::FieldDeclFinder(tooling::Replacements *Replace, bool verbose, bool match_target_file_only ) :
	ReplaceMatchCallback(Replace, "FieldDeclFinder"),
	verbose_( verbose ),
	match_target_file_only_( match_target_file_only )
{}

// Main callback for all matches
void FieldDeclFinder::run(const ast_matchers::MatchFinder::MatchResult & Result) {

	SourceManager *sm = Result.SourceManager;
	const FieldDecl *decl = Result.Nodes.getStmtAs<FieldDecl>("fielddecl");

	if( match_target_file_only_ && !rewriteThisFile(decl, *sm)) {
		return;
	}

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

	if ( verbose_ ) {
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
	FieldDeclaration fd;
	fd.full_name_ = name;
	fd.var_name_ = name.substr( cls.size()+2 );
	fd.type_ = type;
	fd.type_desugared_ = typeD;
	fd.cls_ = cls;
	fd.loc_ = loc;
	fd.range_ = range;
	fd.start_ = Start.second;
	fd.end_ = End.second;

	class_fields_[ cls ].push_back( fd );
}

std::list< std::string >
FieldDeclFinder::all_classes_with_fields() const
{
	std::list< std::string > classnames;
	for ( std::map< std::string, std::list< FieldDeclaration > >::const_iterator
					iter = class_fields_.begin(), iter_end = class_fields_.end();
				iter != iter_end; ++iter ) {
		classnames.push_back( iter->first );
	}
	return classnames;
}


std::list< FieldDeclaration >
FieldDeclFinder::fields_for_class( std::string const & classname ) const
{
	if ( class_fields_.find( classname ) == class_fields_.end() ) {
		//std::cout << "Didn't find any fields for the requested class! " << classname << std::endl;
		std::list< FieldDeclaration > empty_list;
		return empty_list;
	} else {
		return class_fields_.find( classname )->second;
	}
}

clang::ast_matchers::DeclarationMatcher
match_to_field_decl()
{
	using namespace clang::ast_matchers;
	return fieldDecl().bind("fielddecl");

}

void
add_field_decl_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	FieldDeclFinder *cb = new FieldDeclFinder(replacements);
	finder.addMatcher( match_to_field_decl(),	cb);
}
