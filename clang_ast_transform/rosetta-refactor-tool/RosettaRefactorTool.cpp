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
#include "llvm/Support/system_error.h"

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// #define AST_TEST
// #define DEBUG 

bool verbose = true;

using namespace clang;
using namespace clang::ast_matchers;
using namespace llvm;
using clang::tooling::Replacement;

// Code includes
#include "utils.hh"
#include "ast_matchers.hh"
#include "matchers_base.hh"

cl::opt<std::string> BuildPath(
	cl::Positional,
	cl::desc("<build-path>"));

cl::list<std::string> SourcePaths(
	cl::Positional,
	cl::desc("<source0> [... <sourceN>]"),
	cl::OneOrMore);


////////////////////////////////////////////////////////////////////////////////////////////////////

int runMatchers(clang::tooling::RefactoringTool & Tool) {

	using namespace clang::tooling;

	ast_matchers::MatchFinder Finder;
	tooling::Replacements *Replacements = &Tool.getReplacements();

	/////////////////////////////////
	// Include matchers to apply here
	/////////////////////////////////
	
#ifndef AST_TEST
	// Good  matchers
	#include "matchers_typedef.hh"
	#include "matchers_pointer_name.hh"
	#include "matchers_cast_in_assignment.hh"
	#include "matchers_cast_from_new.hh"
	#include "matchers_member_calls.hh"
	#include "matchers_call_operator.hh"
	#include "matchers_dynamic_cast.hh"
#else
	// #include "matchers_op_to_void_ptr.hh"
	#include "matchers_match_test.hh"
#endif

	// Run tool and generate change log
	return Tool.run(newFrontendActionFactory(&Finder));
}
	
////////////////////////////////////////////////////////////////////////////////////////////////////

int saveOutput(clang::tooling::RefactoringTool & Tool) {
	
	using namespace clang::tooling;
	
	LangOptions DefaultLangOptions;
	IntrusiveRefCntPtr<DiagnosticOptions> DiagOpts = new DiagnosticOptions();
	TextDiagnosticPrinter DiagnosticPrinter(llvm::errs(), &*DiagOpts);
	DiagnosticsEngine Diagnostics(
		IntrusiveRefCntPtr<DiagnosticIDs>(new DiagnosticIDs()),
		&*DiagOpts, &DiagnosticPrinter, false);
	SourceManager Sources(Diagnostics, Tool.getFiles());
	Rewriter Rewrite(Sources, DefaultLangOptions);
	const Replacements & Replaces = Tool.getReplacements();
	int result = 0;

	if(BuildPath.empty())
	{
		// Output change log to stdout to be applied later
		for(Replacements::const_iterator I = Replaces.begin(), E = Replaces.end(); I != E; ++I) {
			const Replacement &r = *I;
    	if (r.isApplicable()) {
      	std::string replacementText = r.getReplacementText();
      	replace(replacementText, "\n", "\\n");
      	llvm::outs()
      		<< r.getFilePath() << "\t"
      		<< r.getOffset() << "\t"
      		<< r.getLength() << "\t"
      		<< replacementText << "\n";
      }
		}
	}
	else
	{
		// Don't use runAndSave here to not to overwrite original files
		applyAllReplacements(Replaces, Rewrite);

		size_t slashes = 0;
		std::string outputBaseDir;
		for(size_t i = 0; i < BuildPath.size(); ++i) {
			if(BuildPath[i] == '/') {
				slashes++;
				outputBaseDir = std::string(BuildPath, 0, i);
			}
		}

		// llvm::errs() << "Output directory: " << outputBaseDir << "\n";
				
		for (Rewriter::buffer_iterator I = Rewrite.buffer_begin(),
				E = Rewrite.buffer_end(); I != E; ++I) {

			const FileEntry *Entry = Rewrite.getSourceMgr().getFileEntryForID(I->first);
			const std::string origFileName = Entry->getName();

			// Busywork to determine output path and create directories
			std::string origFileNameRelPath(origFileName);
			size_t this_slashes = 0;
			for(size_t i = 0; i < origFileName.size() && this_slashes < slashes; ++i) {
				if(origFileName[i] == '/') {
					this_slashes++;
					if(this_slashes == slashes)
						origFileNameRelPath = std::string(origFileNameRelPath, i);
				}
			}
			
			std::string outputFileName = outputBaseDir + origFileNameRelPath;

			// Create dir
			for(size_t i = 0; i < outputFileName.size(); ++i) {
				if(i > 0 && outputFileName[i] == '/') {
						std::string path(outputFileName, 0, i);
						mkdir(path.c_str(), 0755);
				}
			}

			// Adapted from clang's Rewriter.cc
			{
				SmallString<256> TempFilename(outputFileName.c_str());
				TempFilename += "-%%%%%%%%";
				std::unique_ptr<llvm::raw_fd_ostream> FileStream;
				int FD;
				bool ok = false;

				if(!llvm::sys::fs::createUniqueFile(TempFilename.str(), FD, TempFilename))
					FileStream.reset(new llvm::raw_fd_ostream(FD, /*shouldClose=*/true));
				if(FileStream) {
					I->second.write(*FileStream); // no error checking on raw_ostream
					ok = !llvm::sys::fs::rename(TempFilename.str(), outputFileName);
					llvm::sys::fs::remove(TempFilename.str()); // if rename fails
				}

				llvm::errs()
					<< origFileName << " -> " << outputFileName << ": "
					<< (ok ? "OK" : "Failed!") << "\n";
					
				if(!ok)
					result++;
			}
		}
	}
	
	return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {

	using namespace clang::tooling;

	llvm::sys::PrintStackTraceOnErrorSignal();
	std::unique_ptr<CompilationDatabase> Compilations(
			FixedCompilationDatabase::loadFromCommandLine(argc, argv));

	cl::ParseCommandLineOptions(argc, argv);
	if(!Compilations) {
		std::string ErrorMessage;
		Compilations.reset(
			CompilationDatabase::loadFromDirectory(BuildPath, ErrorMessage));
		if(!Compilations)
			llvm::report_fatal_error(ErrorMessage);
	}

	
	RefactoringTool Tool(*Compilations, SourcePaths);
	if(int r = runMatchers(Tool))
		return r;
	return saveOutput(Tool);
}
