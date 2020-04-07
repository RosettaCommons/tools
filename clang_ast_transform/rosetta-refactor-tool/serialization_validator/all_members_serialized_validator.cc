// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//===- RosettaRefactorTool.cpp - Clean up for boost ptr conversion -===//
//
// This tool generates a log of source code changes to be applies to
// the original source once all changes have been collected.
// The log is writted to stdout, one change per line.
//
//===-----------------------------------------------------------===//

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

// Code includes
#include "../utils.hh"
#include "../ast_matchers.hh"
#include "../matchers_base.hh"
// #include "RosettaRefactorTool.hh"

// #include "matchers/find/calls.hh"
// #include "matchers/find/constructor_decl.hh"
#include "../matchers/find/field_decl.hh"
// #include "matchers/find/record_decl.hh"
// #include "matchers/find/self_ptr_in_ctor.hh"
#include "../matchers/find/serialization_funcs.hh"

// #include "matchers/code_quality/bad_pointer_casts.hh"
// #include "matchers/code_quality/naked_ptr_op_casts.hh"
// #include "matchers/code_quality/obj_on_stack.hh"

// #include "matchers/rewrite/add_serialization_code.hh"
// #include "matchers/rewrite/call_operator.hh"
// #include "matchers/rewrite/cast_from_new.hh"
// #include "matchers/rewrite/cast_from_new_expr.hh"
// #include "matchers/rewrite/cast_from_new_vardecl.hh"
// #include "matchers/rewrite/cast_in_assignment.hh"
// #include "matchers/rewrite/ctor_initializer.hh"
// #include "matchers/rewrite/datamap_get.hh"
// #include "matchers/rewrite/dynamic_cast.hh"
// #include "matchers/rewrite/member_calls.hh"
// #include "matchers/rewrite/pointer_name.hh"
// #include "matchers/rewrite/real_comparison.hh"
// #include "matchers/rewrite/pose_dynamic_cast.hh"
// #include "matchers/rewrite/typedef.hh"

using namespace clang;
using namespace clang::ast_matchers;
using namespace llvm;
using clang::tooling::Replacement;

////////////////////////////////////////////////////////////////////////////////////////////////////
// Command line options

cl::opt<bool> Debug(
	"debug",
	cl::desc("Enable debugging output"),
	cl::init(false));

cl::opt<bool> Verbose(
	"verbose",
	cl::desc("Enable verbose output"),
	cl::init(false));

cl::opt<bool> Colors(
	"colors",
	cl::desc("Enable color output"),
	cl::init(true));

cl::opt<bool> DangerousRewrites(
	"dangenous-rewrites",
	cl::desc("Enable dangarous matchers in the rewriter"),
	cl::init(false));

// cl::list<std::string> Matchers(
// 	"matchers",
// 	cl::CommaSeparated,
// 	cl::OneOrMore,
// 	cl::desc("Comma-separated list of matchers to apply"));

cl::opt<std::string> BuildPath(
	cl::Positional,
	cl::desc("<build-path>"));

cl::list<std::string> SourcePaths(
	cl::Positional,
	cl::desc("<source0> [... <sourceN>]"),
	cl::OneOrMore);

std::unique_ptr< clang::tooling::CompilationDatabase >
compilation_database_from_commandline( int argc, const char **argv )
{
	using namespace clang::tooling;
	std::unique_ptr< CompilationDatabase > compilations;
	llvm::sys::PrintStackTraceOnErrorSignal();
	{
		std::unique_ptr< CompilationDatabase > compdb_from_cl( FixedCompilationDatabase::loadFromCommandLine(argc, argv));
		compilations.swap( compdb_from_cl );
	}

	cl::ParseCommandLineOptions(argc, argv);
	if(!compilations) {
		std::string ErrorMessage;
		std::unique_ptr< CompilationDatabase > compdb_from_dir = CompilationDatabase::loadFromDirectory(BuildPath, ErrorMessage);
		compilations.swap( compdb_from_dir );
		if(!compilations)
			llvm::report_fatal_error(ErrorMessage);
	}
	return compilations;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

int
identify_datamembers_not_serialized(
	SerializationFunctionFinder const & sff,
	SerializedMemberFinder const & smf,
	FieldDeclFinder const & fdf
)
{
	SerializationFunctionFinder::class_names const & classes_w_serialization_funcs =
		sff.classes_w_serialization_funcs();

	SerializedMemberFinder::data_members const & saved_variables( smf.members_saved() );
	SerializedMemberFinder::data_members const & saved_w_opts_variables( smf.members_saved_w_opts() );
	SerializedMemberFinder::data_members const & loaded_variables( smf.members_loaded() );
	SerializedMemberFinder::data_members const & loaded_w_opts_variables( smf.members_loaded_w_opts() );

	SerializationFunctionFinder::data_members const & exempted_save_variables( sff.exempted_members_from_save() );
	SerializationFunctionFinder::data_members const & exempted_save_w_opts_variables( sff.exempted_members_from_save_w_opts() );
	SerializationFunctionFinder::data_members const & exempted_load_variables( sff.exempted_members_from_load() );
	SerializationFunctionFinder::data_members const & exempted_load_w_opts_variables( sff.exempted_members_from_load_w_opts() );


	std::list< std::string > classes_w_fields = fdf.all_classes_with_fields();
	//for ( std::list< std::string >::const_iterator
	//				iter = classes_w_fields.begin(), iter_end = classes_w_fields.end();
	//			iter != iter_end; ++iter ) {
	//	std::cout << "Class with fields: " << *iter << std::endl;
	//}

	// for ( std::set< std::pair< std::string, std::string > >::const_iterator
	// 				iter = saved_variables.begin(), iter_end = saved_variables.end();
	// 			iter != iter_end; ++iter ) {
	// 	std::cout << "Saved variable " << iter->first << "::" << iter->second << std::endl;
	// }
	//
	// for ( std::set< std::pair< std::string, std::string > >::const_iterator
	// 				iter = loaded_variables.begin(), iter_end = loaded_variables.end();
	// 			iter != iter_end; ++iter ) {
	// 	std::cout << "Loaded variable " << iter->first << "::" << iter->second << std::endl;
	// }
	// for ( std::set< std::pair< std::string, std::string > >::const_iterator
	// 				iter = saved_w_opts_variables.begin(), iter_end = saved_w_opts_variables.end();
	// 			iter != iter_end; ++iter ) {
	// 	std::cout << "Saved_w_opts variable " << iter->first << "::" << iter->second << std::endl;
	// }
	//
	// for ( std::set< std::pair< std::string, std::string > >::const_iterator
	// 				iter = loaded_w_opts_variables.begin(), iter_end = loaded_w_opts_variables.end();
	// 			iter != iter_end; ++iter ) {
	// 	std::cout << "Loaded_w_opts variable " << iter->first << "::" << iter->second << std::endl;
	// }

	bool any_missed = false;
	for ( SerializationFunctionFinder::class_names::const_iterator
			class_iter = classes_w_serialization_funcs.begin(),
			class_iter_end = classes_w_serialization_funcs.end();
			class_iter != class_iter_end; ++class_iter ) {
		//std::cout << "Examining class " << *class_iter << std::endl;
		std::list< FieldDeclaration > fields = fdf.fields_for_class( *class_iter );
		for ( std::list< FieldDeclaration >::const_iterator
				field_iter = fields.begin(), field_iter_end = fields.end();
				field_iter != field_iter_end; ++field_iter ) {
			//std::cout << "  Examining field " << field_iter->var_name_ << std::endl;
			std::pair< std::string, std::string > class_field_pair( std::make_pair( *class_iter, field_iter->var_name_ ) );
			if ( saved_variables.find( class_field_pair ) == saved_variables.end() &&
					exempted_save_variables.find( class_field_pair ) == exempted_save_variables.end() ) {
				std::cout << "Data member of class " << *class_iter << " named " << field_iter->var_name_ << " was not saved" << std::endl;
				any_missed = true;
			}
			if ( loaded_variables.find( class_field_pair ) == loaded_variables.end() &&
					exempted_load_variables.find( class_field_pair ) == exempted_load_variables.end() ) {
				std::cout << "Data member of class " << *class_iter << " named " << field_iter->var_name_ << " was not loaded" << std::endl;
				any_missed = true;
			}
		}
	}

	if ( ! saved_w_opts_variables.empty() || ! loaded_w_opts_variables.empty() ||
		! exempted_save_w_opts_variables.empty() || ! exempted_load_w_opts_variables.empty() ) {

		std::set< std::string > classes_w_save_load_w_opts;
		for ( auto class_var_pair: saved_w_opts_variables ) { classes_w_save_load_w_opts.insert( class_var_pair.first ); }
		for ( auto class_var_pair: loaded_w_opts_variables ) { classes_w_save_load_w_opts.insert( class_var_pair.first ); }
		for ( auto class_var_pair: exempted_load_w_opts_variables ) { classes_w_save_load_w_opts.insert( class_var_pair.first ); }
		for ( auto class_var_pair: exempted_save_w_opts_variables ) { classes_w_save_load_w_opts.insert( class_var_pair.first ); }

  	for ( SerializationFunctionFinder::class_names::const_iterator
  			class_iter = classes_w_serialization_funcs.begin(),
  			class_iter_end = classes_w_serialization_funcs.end();
  			class_iter != class_iter_end; ++class_iter ) {
  		//std::cout << "Examining class " << *class_iter << std::endl;

			if ( ! classes_w_save_load_w_opts.count( *class_iter ) ) continue;

  		std::list< FieldDeclaration > fields = fdf.fields_for_class( *class_iter );
  		for ( std::list< FieldDeclaration >::const_iterator
  				field_iter = fields.begin(), field_iter_end = fields.end();
  				field_iter != field_iter_end; ++field_iter ) {
  			//std::cout << "  Examining field " << field_iter->var_name_ << std::endl;
  			std::pair< std::string, std::string > class_field_pair( std::make_pair( *class_iter, field_iter->var_name_ ) );
  			if ( saved_w_opts_variables.find( class_field_pair ) == saved_w_opts_variables.end() &&
  					exempted_save_w_opts_variables.find( class_field_pair ) == exempted_save_w_opts_variables.end() ) {
  				std::cout << "Data member of class " << *class_iter << " named " << field_iter->var_name_ << " was not saved in save_with_options" << std::endl;
  				any_missed = true;
  			}
  			if ( loaded_w_opts_variables.find( class_field_pair ) == loaded_w_opts_variables.end() &&
  					exempted_load_w_opts_variables.find( class_field_pair ) == exempted_load_w_opts_variables.end() ) {
  				std::cout << "Data member of class " << *class_iter << " named " << field_iter->var_name_ << " was not loaded in load_with_options" << std::endl;
  				any_missed = true;
  			}
  		}
		}
	}


	return any_missed;
}

int main(int argc, const char **argv) {

	// std::cout << "RUNNING!" << std::endl;

	using namespace clang::tooling;
	std::unique_ptr< CompilationDatabase > compilations = compilation_database_from_commandline( argc, argv );
	clang::tooling::RefactoringTool * tool = new RefactoringTool(*compilations, SourcePaths);

	ast_matchers::MatchFinder finder;
	tooling::Replacements * replacements = &(tool->getReplacements());

	SerializationFunctionFinder sff( replacements, false );
	SerializedMemberFinder smf( replacements, false );
	FieldDeclFinder fdf(
		replacements,
		false /* not verbose */,
		false /*look for fields in all files, not just the target file */ );

	finder.addMatcher( match_to_serialization_method_definition(), &sff );
	finder.addMatcher( match_to_serialized_data_members(), &smf );
	finder.addMatcher( match_to_externally_serialized_data_members(), &smf );
	finder.addMatcher( match_to_field_decl(), &fdf );

	int run_result = tool->run( newFrontendActionFactory(&finder).get() );
	if ( run_result != 0 ) {
		return run_result;
	}

	int some_members_forgotten = identify_datamembers_not_serialized( sff, smf, fdf );
	return some_members_forgotten;
}
