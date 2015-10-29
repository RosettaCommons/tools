// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
// #include "matchers/find/field_decl.hh"
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

int main(int argc, const char **argv) {

	using namespace clang::tooling;
	std::unique_ptr< CompilationDatabase > compilations = compilation_database_from_commandline( argc, argv );
	clang::tooling::RefactoringTool * tool = new RefactoringTool(*compilations, SourcePaths);

	ast_matchers::MatchFinder finder;
	tooling::Replacements * replacements = &(tool->getReplacements());

	SerializedMemberFinder smf( replacements, true );
	SerializationFunctionFinder sff( replacements, true );

	finder.addMatcher( match_to_serialization_method_definition(), &sff );
	finder.addMatcher( match_to_serialized_data_members(), &smf );
	finder.addMatcher( match_to_externally_serialized_data_members(), &smf );

	return tool->run( newFrontendActionFactory(&finder).get() );
}
