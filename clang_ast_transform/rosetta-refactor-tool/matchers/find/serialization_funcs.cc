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

SerializationFunctionFinder::SerializationFunctionFinder( clang::tooling::Replacements * replace, bool verbose ) :
	ReplaceMatchCallback( replace, "SerializationFuncFinder" ),
	verbose_( verbose )
{}

SerializationFunctionFinder::~SerializationFunctionFinder() {}

// Main callback for all matches
void
SerializationFunctionFinder::run( clang::ast_matchers::MatchFinder::MatchResult const & result)
{
	using namespace clang;

	SourceManager *sm = result.SourceManager;
	FunctionTemplateDecl const * temp_decl = result.Nodes.getStmtAs<FunctionTemplateDecl>("temp_decl");
	CXXMethodDecl const * meth_decl = result.Nodes.getStmtAs<CXXMethodDecl>("method_decl");
	ReferenceType const * ref_to_serialized_class = result.Nodes.getStmtAs<ReferenceType>("serialized_class_type");

	if(!rewriteThisFile(temp_decl, *sm))
		return;

	std::string classname;
	if ( meth_decl ) {
		classname = meth_decl->getThisType( *result.Context )->getPointeeType().getCanonicalType().getUnqualifiedType().getAsString();
	} else {
		classname = ref_to_serialized_class->getPointeeType().getCanonicalType().getUnqualifiedType().getAsString();
	}

	if ( verbose_ ) {
		if ( classes_w_serialization_funcs_.find( classname ) == classes_w_serialization_funcs_.end() ) {
			std::cout << "Found serialization function for record: " << classname << std::endl;
		}
		//temp_decl->dump();
	} else { std::cout << "SFF being very quiet" << std::endl; }

	classes_w_serialization_funcs_.insert( classname );
}

SerializationFunctionFinder::class_names const &
SerializationFunctionFinder::classes_w_serialization_funcs() const
{
	return classes_w_serialization_funcs_;
}

clang::ast_matchers::DeclarationMatcher
match_to_serialization_method_definition()
{
	using namespace clang::ast_matchers;
	return functionTemplateDecl(
		anyOf( hasName( "save"), hasName("load"),hasName( "load_and_construct") )
		,hasDescendant( parmVarDecl( hasType( referenceType(pointee(asString("Archive" ))))))
		,anyOf(
			hasDescendant( methodDecl(isDefinition()).bind( "method_decl" )),
			hasDescendant( functionDecl(
				parameterCountIs(2)
				,hasParameter(1,parmVarDecl(hasType( referenceType().bind("serialized_class_type"))))
			)))
	).bind( "temp_decl" );

}


/*
	Find records with save and load functions and also the member variables
	that are mentioned inside of these functions.
*/


SerializedMemberFinder::SerializedMemberFinder(clang::tooling::Replacements *replace, bool verbose ) :
	ReplaceMatchCallback(replace, "RecordDeclFinder"),
	verbose_( verbose )
{}

SerializedMemberFinder::~SerializedMemberFinder() {}

// Main callback for all matches
void SerializedMemberFinder::run(const clang::ast_matchers::MatchFinder::MatchResult &result) {
	using namespace clang;

	SourceManager *sm = result.SourceManager;
	// const CXXOperatorCallExpr * opcall = result.Nodes.getStmtAs<CXXOperatorCallExpr>("op_call");
	const MemberExpr * member_var = result.Nodes.getStmtAs<MemberExpr>("member");
	const CXXMethodDecl * methdecl = result.Nodes.getStmtAs<CXXMethodDecl>("methdecl");
	const FunctionTemplateDecl * tempfunc = result.Nodes.getStmtAs<FunctionTemplateDecl>("tempfunc");
	const ReferenceType * ref_to_serialized_class = result.Nodes.getStmtAs<ReferenceType>("serialized_class_type");


	if(!rewriteThisFile(member_var, *sm))
		return;
	std::string classname;
	if ( methdecl ) {
	 classname = methdecl->getThisType( *result.Context )->getPointeeType().getCanonicalType().getUnqualifiedType().getAsString();
	} else {
		classname = ref_to_serialized_class->getPointeeType().getCanonicalType().getUnqualifiedType().getAsString();
		//std::cout << "ref_to_serialized_class:" << classname << std::endl;
		//vardecl->dump();
		//tempfunc->dump();
		//std::cout << "\n\n";
	}
	std::string varname = member_var->getMemberNameInfo().getAsString();
	std::string funcname = tempfunc->getNameAsString();


  if ( funcname == "save" && save_variables_.find( std::make_pair( classname, varname )) != save_variables_.end() ) {
		return;
	}
  if ( funcname == "load" && load_variables_.find( std::make_pair( classname, varname )) != load_variables_.end() ) {
		return;
	}

	//std::cout << "found one" << std::endl; opcall->dump(); std::cout << "\n"; member_var->dump(); std::cout << "\n" << std::endl;
	if ( verbose_ ) {
		//std::cout << "HIT!" << std::endl;
		//std::cout << "templated function:" << tempfunc << " " << member_var << " " << opcall << " " << methdecl << std::endl;
		std::cout << tempfunc->getNameAsString() << " " << classname << " " << varname << std::endl;
		// tempfunc->dump();
		// std::cout << "\n";
		// methdecl->dump();
		// std::cout << std::endl;
		// std::cout << "parm: " << parmvardecl->getOriginalType().getAsString() << " " << parmvardecl->isTemplateDecl() << std::endl;
		//
		//
		// std::cout << " " << varname << " " << classname;
		//std::cout << "----\n" << std::endl;
	}// else { std::cout << "SMF being very quiet" << std::endl; }

	if ( funcname == "save" )	save_variables_.insert( std::make_pair(classname, varname ));
	if ( funcname == "load" ) load_variables_.insert( std::make_pair(classname, varname ));
}

SerializedMemberFinder::class_names  const & SerializedMemberFinder::classes_w_serialization_funcs() const
{
	return classes_w_serialization_funcs_;
}

SerializedMemberFinder::data_members const & SerializedMemberFinder::members_saved() const
{
	return save_variables_;
}

SerializedMemberFinder::data_members const & SerializedMemberFinder::members_loaded() const
{
	return load_variables_;
}


//   | |-CXXMethodDecl 0x4aab890 <line:11:3, line:14:3> line:11:8 save 'void (class cereal::BinaryOutputArchive &) const'
//   | | |-TemplateArgument type 'class cereal::BinaryOutputArchive'
//   | | |-ParmVarDecl 0x4aab7d0 <col:14, col:24> col:24 arc 'class cereal::BinaryOutputArchive &'
//   | | `-CompoundStmt 0x4aac1b8 <col:36, line:14:3>
//   | |   |-CXXOperatorCallExpr 0x4aabd40 <line:12:5, col:17> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue
//   | |   | |-ImplicitCastExpr 0x4aabd28 <col:5, col:17> 'class cereal::BinaryOutputArchive &(*)(const int &&)' <FunctionToPointerDecay>
//   | |   | | `-DeclRefExpr 0x4aabca8 <col:5, col:17> 'class cereal::BinaryOutputArchive &(const int &&)' lvalue CXXMethod 0x4aabba0 'operator()' 'class cereal::BinaryOutputArchive &(const int &&)'
//   | |   | |-ImplicitCastExpr 0x4aabd88 <col:5> 'class cereal::OutputArchive<class cereal::BinaryOutputArchive, 1>' lvalue <UncheckedDerivedToBase (OutputArchive)>
//   | |   | | `-DeclRefExpr 0x4aab998 <col:5> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArc<hive' lvalue ParmVar 0x4aab7d0 'arc' 'class cereal::BinaryOutputArchive &'
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
// 	).bind( "op_call" ), new SerializedMemberFinder(Replacements) );

// Finder.addMatcher(
// 	cxxOperatorCallExpr(
// 		hasDescendant( memberExpr().bind("member")),
// 		hasAncestor( functionTemplateDecl(
// 			hasName( "save" )
// 			//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::BinaryOutputArchive")))))
// 			,has( templateTypeParmDecl())
// 		).bind("savemethod") )
// 	).bind( "op_call" ), new SerializedMemberFinder(Replacements) );



// This for some reason reports each match multiple times! Grrr
//

clang::ast_matchers::StatementMatcher
match_to_serialized_data_members() {
	using namespace clang::ast_matchers;

	return
		memberExpr(
			hasAncestor( functionTemplateDecl(
				anyOf(hasName("save"),hasName("load")),
				//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::JSONOutputArchive")))))
				//,has( methodDecl( hasDescendant( parmVarDecl( hasType(referenceType()) ).bind("vardecl") )).bind("methdecl") )
				has( methodDecl( hasDescendant( parmVarDecl( hasType(referenceType(pointee(asString("Archive")) )) ).bind("vardecl") )).bind("methdecl") )
			).bind("tempfunc") )
		).bind("member");
}

// look for member variables that are referred to in "external" (i.e. non-member function) save and load methods
clang::ast_matchers::StatementMatcher
match_to_externally_serialized_data_members() {
	using namespace clang::ast_matchers;

	return
		memberExpr(
			hasAncestor( functionTemplateDecl(
				anyOf(hasName("save"),hasName("load")),
				//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::JSONOutputArchive")))))
				//,has( methodDecl( hasDescendant( parmVarDecl( hasType(referenceType()) ).bind("vardecl") )).bind("methdecl") )
				hasDescendant( functionDecl(
					hasDescendant( parmVarDecl( hasType(referenceType(pointee(asString("Archive")) )) ) )
					,parameterCountIs( 2 )
					,hasParameter(1,parmVarDecl(hasType( referenceType().bind("serialized_class_type"))))
				))
			).bind("tempfunc") )
		).bind("member");
}


void
add_serialization_func_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements ) {
	using namespace clang;
	std::cout << "adding serialized member finder" << std::endl;
	SerializedMemberFinder * smf = new SerializedMemberFinder( replacements, true );
	finder.addMatcher( match_to_serialized_data_members(), smf );
	finder.addMatcher( match_to_externally_serialized_data_members(), smf );
}


//Finder.addMatcher(
//	recordDecl(
//		hasParent(
//			decl().bind("parent")
//		)
//	).bind("recorddecl"),
//	new RecordDeclFinder(Replacements));
