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

cl::opt<bool> DangarousRewrites(
	"danganous-rewrites",
	cl::desc("Enable dangarous matchers in the rewriter"),
	cl::init(false));

cl::list<std::string> Matchers(
	"matchers",
	cl::CommaSeparated,
	cl::OneOrMore,
	cl::desc("Comma-separated list of matchers to apply"));

cl::opt<std::string> BuildPath(
	cl::Positional,
	cl::desc("<build-path>"));

cl::list<std::string> SourcePaths(
	cl::Positional,
	cl::desc("<source0> [... <sourceN>]"),
	cl::OneOrMore);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Code includes

#include "utils.hh"
#include "ast_matchers.hh"
#include "matchers_base.hh"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Main tool class

class RosettaRefactorTool {

public:
	RosettaRefactorTool(int argc, const char **argv);
	int Run();
	int runMatchers();
	int saveOutput();

private:
	clang::tooling::RefactoringTool * Tool;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {

	RosettaRefactorTool rrt(argc, argv);
	return rrt.Run();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tool c'tor
RosettaRefactorTool::RosettaRefactorTool(int argc, const char **argv)
{
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

	Tool = new RefactoringTool(*Compilations, SourcePaths);
}

/// Main tool runner
int RosettaRefactorTool::Run() {
	int r = runMatchers();
	if(r)
		return r;
	return saveOutput();
}
	
/// Run selected matchers
int RosettaRefactorTool::runMatchers() {

	ast_matchers::MatchFinder Finder;
	tooling::Replacements *Replacements = &(Tool->getReplacements());

	for(cl::list<std::string>::const_iterator it(Matchers.begin()), end(Matchers.end());
			it != end; ++it) {
		std::string matcher( *it );
		if(Verbose)
			llvm::errs() << color("gray") << "Applying matcher: " << matcher << color("") << "\n";

		/////////////////////////////////
		// Include matchers to apply here
		/////////////////////////////////

		// Finders
		if(matcher == "find" || matcher == "find_calls") {
			#include "matchers/find_calls.hh"
		}
		if(matcher == "find" || matcher == "find_naked_ptr_op_casts") {
			#include "matchers/find_naked_ptr_op_casts.hh"
		}

		// Code quality checkers
		if(matcher == "quality_check" || matcher == "find_bad_pointer_casts") {
			#include "matchers/find_bad_pointer_casts.hh"
		}

		// Rewriters
		if(matcher == "rewrite" || matcher == "rewrite_typedef") {
			#include "matchers/typedef.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_pointer_name") {
			#include "matchers/pointer_name.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_cast_from_new") {
			#include "matchers/cast_from_new.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_cast_in_assignment") {
			#include "matchers/cast_in_assignment.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_ctor_initializer") {
			#include "matchers/ctor_initializer.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_dynamic_cast") {
			#include "matchers/dynamic_cast.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_datamap_get") {
			#include "matchers/datamap_get.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_pose_dynamic_cast") {
			#include "matchers/pose_dynamic_cast.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_call_operator") {
			#include "matchers/call_operator.hh"
		}
		if(matcher == "rewrite" || matcher == "rewrite_member_calls") {
			#include "matchers/member_calls.hh"
		}
	
		if(matcher == "rewrite_not_operator") {
			// Not needed; see comment in file
			#include "matchers/not_operator.hh"
		}
	}

	// Run tool and generate change log
	return Tool->run(clang::tooling::newFrontendActionFactory(&Finder));
}

//////////////////////////////////////////////////////////////////////////////////////////////////

/// Save rewritten output to files or output to stdout
int RosettaRefactorTool::saveOutput() {
	
	using namespace clang::tooling;
	
	LangOptions DefaultLangOptions;
	IntrusiveRefCntPtr<DiagnosticOptions> DiagOpts = new DiagnosticOptions();
	TextDiagnosticPrinter DiagnosticPrinter(llvm::errs(), &*DiagOpts);
	DiagnosticsEngine Diagnostics(
		IntrusiveRefCntPtr<DiagnosticIDs>(new DiagnosticIDs()),
		&*DiagOpts, &DiagnosticPrinter, false);
	SourceManager Sources(Diagnostics, Tool->getFiles());
	Rewriter Rewrite(Sources, DefaultLangOptions);
	const Replacements & Replaces = Tool->getReplacements();
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
