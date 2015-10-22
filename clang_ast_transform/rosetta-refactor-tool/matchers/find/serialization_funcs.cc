// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include "serialization_funcs.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"
#include "clang/Lex/Lexer.h"
#include "clang/Tooling/CompilationDatabase.h"
#include "clang/Tooling/Refactoring.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Frontend/TextDiagnosticPrinter.h"

#include "llvm/Support/CommandLine.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/Path.h"
#include "llvm/Support/Signals.h"
#include "llvm/Support/raw_ostream.h"
//#include "llvm/Support/system_error.h"

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>


/*
	Find records with save and load functions and also the member variables
	that are mentioned inside of these functions.
*/


SerializationFuncFinder::SerializationFuncFinder(clang::tooling::Replacements *Replace) :
	ReplaceMatchCallback(Replace, "RecordDeclFinder")
{}

// Main callback for all matches
void SerializationFuncFinder::run(const clang::ast_matchers::MatchFinder::MatchResult &Result) {
	using namespace clang;

	SourceManager *sm = Result.SourceManager;
	const CXXOperatorCallExpr * opcall = Result.Nodes.getStmtAs<CXXOperatorCallExpr>("op_call");
	const MemberExpr * member_var = Result.Nodes.getStmtAs<MemberExpr>("member");
	//const CXXMethodDecl * save_method = Result.Nodes.getStmtAs<CXXMethodDecl>("savemethod");
	const FunctionTemplateDecl * method = Result.Nodes.getStmtAs<FunctionTemplateDecl>("savemethod");

	if(!rewriteThisFile(opcall, *sm))
		return;
	//std::cout << "found one" << std::endl; opcall->dump(); std::cout << "\n"; member_var->dump(); std::cout << "\n" << std::endl;
	std::cout << "method:" << method << " " << member_var << " " << opcall << " ";
	method->dump();

	std::cout << std::endl;
	//std::cout << save_method->getThisType( *Result.Context )->getPointeeType().getCanonicalType().getUnqualifiedType().getAsString();
	//
	std::cout << " " << member_var->getMemberNameInfo().getAsString();
	std::cout << std::endl;

}

//   | |-CXXMethodDecl 0x4aab890 <line:11:3, line:14:3> line:11:8 save 'void (class cereal::BinaryOutputArchive &) const'
//   | | |-TemplateArgument type 'class cereal::BinaryOutputArchive'
//   | | |-ParmVarDecl 0x4aab7d0 <col:14, col:24> col:24 arc 'class cereal::BinaryOutputArchive &'
//   | | `-CompoundStmt 0x4aac1b8 <col:36, line:14:3>
//   | |   |-CXXOperatorCallExpr 0x4aabd40 <line:12:5, col:17> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue
//   | |   | |-ImplicitCastExpr 0x4aabd28 <col:5, col:17> 'class cereal::BinaryOutputArchive &(*)(const int &&)' <FunctionToPointerDecay>
//   | |   | | `-DeclRefExpr 0x4aabca8 <col:5, col:17> 'class cereal::BinaryOutputArchive &(const int &&)' lvalue CXXMethod 0x4aabba0 'operator()' 'class cereal::BinaryOutputArchive &(const int &&)'
//   | |   | |-ImplicitCastExpr 0x4aabd88 <col:5> 'class cereal::OutputArchive<class cereal::BinaryOutputArchive, 1>' lvalue <UncheckedDerivedToBase (OutputArchive)>
//   | |   | | `-DeclRefExpr 0x4aab998 <col:5> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue ParmVar 0x4aab7d0 'arc' 'class cereal::BinaryOutputArchive &'
//   | |   | `-MemberExpr 0x4aab220 <col:10> 'const int' lvalue ->myint_ 0x4aaae80
//   | |   |   `-CXXThisExpr 0x4aab208 <col:10> 'const class MyClass *' this
//   | |   `-CXXOperatorCallExpr 0x4aac150 <line:13:5, col:19> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue
//   | |     |-ImplicitCastExpr 0x4aac138 <col:5, col:19> 'class cereal::BinaryOutputArchive &(*)(const float &&)' <FunctionToPointerDecay>
//   | |     | `-DeclRefExpr 0x4aac0b8 <col:5, col:19> 'class cereal::BinaryOutputArchive &(const float &&)' lvalue CXXMethod 0x4aabfb0 'operator()' 'class cereal::BinaryOutputArchive &(const float &&)'
//   | |     |-ImplicitCastExpr 0x4aac198 <col:5> 'class cereal::OutputArchive<class cereal::BinaryOutputArchive, 1>' lvalue <UncheckedDerivedToBase (OutputArchive)>
//   | |     | `-DeclRefExpr 0x4aabda8 <col:5> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue ParmVar 0x4aab7d0 'arc' 'class cereal::BinaryOutputArchive &'
//   | |     `-MemberExpr 0x4aab2c0 <col:10> 'const float' lvalue ->myfloat_ 0x4aaaee0
//   | |       `-CXXThisExpr 0x4aab2a8 <col:10> 'const class MyClass *' this

// This works fine, but cannot report on multiple data members serialized within a single call operator, e.g. arc( myfloat_, mychar_ );
// Finder.addMatcher(
// 	operatorCallExpr(
// 		hasDescendant( memberExpr().bind("member")),
// 		hasAncestor( methodDecl(
// 			hasName( "save" )
// 			//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::BinaryOutputArchive")))))
// 			//hasTemplateArgument(0,refersToType(asString("class cereal::BinaryOutputArchive")))
// 		).bind("savemethod") )
// 	).bind( "op_call" ), new SerializationFuncFinder(Replacements) );

// Finder.addMatcher(
// 	cxxOperatorCallExpr(
// 		hasDescendant( memberExpr().bind("member")),
// 		hasAncestor( functionTemplateDecl(
// 			hasName( "save" )
// 			//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::BinaryOutputArchive")))))
// 			,has( templateTypeParmDecl())
// 		).bind("savemethod") )
// 	).bind( "op_call" ), new SerializationFuncFinder(Replacements) );



// This for some reason reports each match multiple times! Grrr
//

void
add_serialization_func_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements ) {
	using namespace clang;
	using namespace clang::ast_matchers;

	finder.addMatcher(
		memberExpr( hasParent( cxxOperatorCallExpr(
			hasAncestor( functionTemplateDecl(
				hasName("save"),
				//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::JSONOutputArchive")))))
				has( parmVarDecl( hasName( "Archive" ) ) )
			).bind("savemethod") )).bind("op_call"))
		).bind("member"),
		new SerializationFuncFinder( replacements ) );
}
//Finder.addMatcher(
//	recordDecl(
//		hasParent(
//			decl().bind("parent")
//		)
//	).bind("recorddecl"),
//	new RecordDeclFinder(Replacements));
