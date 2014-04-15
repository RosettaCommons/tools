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
#include <iostream>
#include <fstream>
#include <sstream>

#include <sys/stat.h>

using namespace clang;
using namespace clang::ast_matchers;
using namespace llvm;
using clang::tooling::Replacement;

bool verbose = true;

cl::opt<std::string> BuildPath(
	cl::Positional,
	cl::desc("<build-path>"));

cl::list<std::string> SourcePaths(
	cl::Positional,
	cl::desc("<source0> [... <sourceN>]"),
	cl::OneOrMore);

////////////////////////////////////////////////////////////////////////////////////////////////////
// String utils

namespace {

void replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
}

std::string trim(const std::string & s, const char *whitespace =0) {
	size_t start, end;
	if(!whitespace)
		whitespace = " \r\n\t";
	for(start = 0; start < s.length() && strchr(whitespace, s[start]); start++);
	for(end = s.length() -1; end > 0  && strchr(whitespace, s[end]);   end--);
	return std::string(s, start, end);
}

bool endsWith(const std::string& a, const std::string& b) {
	if (b.size() > a.size()) return false;
	return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}

bool beginsWith(const std::string& a, const std::string& b) {
	return (a.compare(0, b.length(), b) == 0);
}

} // namespace

////////////////////////////////////////////////////////////////////////////////////////////////////
// Own matchers

namespace clang {
namespace ast_matchers {

AST_MATCHER(CastExpr, isNullToPointer) {
  return Node.getCastKind() == CK_NullToPointer ||
    Node.getCastKind() == CK_NullToMemberPointer;
}

AST_MATCHER(CastExpr, isNonNoopCast) {
  return Node.getCastKind() != CK_NoOp;
}

AST_MATCHER(CastExpr, isLValueToRValueCast) {
  return Node.getCastKind() == CK_LValueToRValue;
}

AST_MATCHER(CastExpr, isConstructorConversionCast) {
  return Node.getCastKind() == CK_ConstructorConversion;
}

AST_MATCHER(Decl, isTypedefDecl) {
  //return (Node.getKind() == NK_Typedef);
  return !strcmp(Node.getDeclKindName(), "Typedef");
}

AST_MATCHER(Expr, isUtilityPointer) {
	std::string type = QualType::getAsString(Node.getType().split());
	if(beginsWith(type, "const ")) // trim const
		type = std::string(type, 6);
	return
		beginsWith(type, "class utility::pointer::owning_ptr") ||
		beginsWith(type, "class utility::pointer::access_ptr");
}

AST_MATCHER(Type, sugaredNullptrType) {
  const Type *DesugaredType = Node.getUnqualifiedDesugaredType();
  if (const BuiltinType *BT = dyn_cast<BuiltinType>(DesugaredType))
    return BT->getKind() == BuiltinType::NullPtr;
  return false;
}

} // end namespace ast_matchers
} // end namespace clang

////////////////////////////////////////////////////////////////////////////////////////////////////
// Clang Utils

namespace {
	
// Returns the text that makes up 'node' in the source.
// Returns an empty string if the text cannot be found.
static std::string getText(
	const SourceManager &SourceManager,
	const SourceLocation &StartSpellingLocation, 
	const SourceLocation &EndSpellingLocation
) {
  if (!StartSpellingLocation.isValid() || !EndSpellingLocation.isValid()) {
    return std::string();
  }
  bool Invalid = true;
  const char *Text =
      SourceManager.getCharacterData(StartSpellingLocation, &Invalid);
  if (Invalid) {
    return std::string();
  }
  std::pair<FileID, unsigned> Start =
      SourceManager.getDecomposedLoc(StartSpellingLocation);
  std::pair<FileID, unsigned> End =
      SourceManager.getDecomposedLoc(Lexer::getLocForEndOfToken(
          EndSpellingLocation, 0, SourceManager, LangOptions()));
  if (Start.first != End.first) {
    // Start and end are in different files.
    return std::string();
  }
  if (End.second < Start.second) {
    // Shuffling text with macros may cause this.
    return std::string();
  }
  return std::string(Text, End.second - Start.second);
}

template <typename T>
static std::string getText(const SourceManager &SourceManager, const T *Node) {
	if(!Node)
		return std::string();
	return getText(SourceManager,
		SourceManager.getSpellingLoc(Node->getLocStart()), 
		SourceManager.getSpellingLoc(Node->getLocEnd()));
}

template <typename T, typename U>
static std::string getTextToDelim(const SourceManager &SourceManager, const T *Node, const U *DelimNode) {
	if(!Node || !DelimNode)
		return std::string();
	return getText(SourceManager,
		SourceManager.getSpellingLoc(Node->getLocStart()), 
		SourceManager.getSpellingLoc(DelimNode->getLocStart()));
}

template <typename T>
void dumpRewrite(SourceManager & sm, T * node, const std::string & newCodeStr) {
	
	if(!verbose)
		return;
		
	const std::string origCodeStr = getText(sm, node);
	const std::string locStr( node->getSourceRange().getBegin().printToString(sm) );
		
	llvm::errs() 
		<< "LOCATION:\t" << locStr << "\n" 
		<< "OLD CODE:\t" << origCodeStr << "\n"
		<< "NEW CODE:\t" << newCodeStr << "\n"
		<< "\n";
}

} // anon namespace

////////////////////////////////////////////////////////////////////////////////////////////////////

class ReplaceMatchCallback : public ast_matchers::MatchFinder::MatchCallback {
public:
	ReplaceMatchCallback(tooling::Replacements *Replace)
			: Replace(Replace) {}

private:
	tooling::Replacements *Replace;

protected:
	template <typename T>
	void doRewrite(SourceManager &sm, T * node, const std::string & origCode, const std::string & newCode) {
		if(origCode == newCode)
			return;
		dumpRewrite(sm, node, newCode);
		Replace->insert(Replacement(sm, node, newCode));
	}
};
	
////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace owning_ptr/access_ptr in typedefs
//   typedef utility::pointer::access_ptr< Some > SomeAP
//   typedef utility::pointer::owning_ptr< Some > SomeOP
//   typedef utility::pointer::access_ptr< Some const > SomeCAP
//   typedef utility::pointer::owning_ptr< Some const > SomeCOP

class RewriteTypedefDecl : public ReplaceMatchCallback {
public:
	RewriteTypedefDecl(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Decl *decl = Result.Nodes.getStmtAs<Decl>("typedefdecl");

		const FullSourceLoc FullLocation = FullSourceLoc(decl->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;
		
		std::string origCode( getText(sm, decl) );
		std::string newCode( origCode );

		replace(newCode, "owning_ptr", "shared_ptr");
		replace(newCode, "access_ptr", "weak_ptr");

		doRewrite(sm, decl, origCode, newCode);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in assignments
//   SomeAP a = new Some;
//   SomeAP a = 0;
//   SomeAP a = NULL;

class RewriteImplicitAssignmentOPCast : public ReplaceMatchCallback {
public:
	RewriteImplicitAssignmentOPCast(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const CXXOperatorCallExpr *opercallexpr = Result.Nodes.getStmtAs<CXXOperatorCallExpr>("opercallexpr");
		const CastExpr *castexpr = Result.Nodes.getStmtAs<CastExpr>("castexpr");
		const MemberExpr *memberexpr = Result.Nodes.getStmtAs<MemberExpr>("memberexpr");
		
		const FullSourceLoc FullLocation = FullSourceLoc(opercallexpr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;

		// Get original code			
		const std::string type = QualType::getAsString(memberexpr->getType().split());
		const std::string origCode = getText(sm, opercallexpr);
		std::string leftSideCode = getTextToDelim(sm, opercallexpr, castexpr);
		std::string rightSideCode(origCode, origCode.length() - leftSideCode.length());
		
		leftSideCode = trim(leftSideCode, " \t\n=");
		rightSideCode = trim(rightSideCode, " \t\n=");

		std::string newCode;
		if(rightSideCode == "0" || rightSideCode == "NULL") {
			//newCode = leftCode + " " + typedef_name + "()";
			newCode = leftSideCode + ".reset()";
		} else {
			newCode = leftSideCode + " = " + type + "( " + rightSideCode + " )";
		}

		doRewrite(sm, opercallexpr, origCode, newCode);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in returns
//   SomeAP function() { return new Some; }

class RewriteImplicitReturnOPCast : public ReplaceMatchCallback {
public:
	RewriteImplicitReturnOPCast(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const CXXBindTemporaryExpr *temporaryexpr = Result.Nodes.getStmtAs<CXXBindTemporaryExpr>("temporaryexpr");

		const FullSourceLoc FullLocation = FullSourceLoc(temporaryexpr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;
		
		std::string type = QualType::getAsString( temporaryexpr->getType().split() );
		std::string origCode( getText(sm, temporaryexpr) );
		
		// If original code was 0 then just leave it empty
		std::string newCode(origCode);
		if(origCode == "0" || origCode == "NULL")
			newCode = type + "()";
		else
			newCode = type + "( " + origCode + " )";

		doRewrite(sm, temporaryexpr, origCode, newCode);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////

int runMatchers(clang::tooling::RefactoringTool & Tool) {

	using namespace clang::tooling;

	ast_matchers::MatchFinder Finder;

	// Typedefs for access_ptr and owning_ptr
	RewriteTypedefDecl RewriteTypedefDeclCallback(&Tool.getReplacements());
	Finder.addMatcher(
		decl(isTypedefDecl()).bind("typedefdecl"),
		&RewriteTypedefDeclCallback);

	// operator= calls with implicit conversions to OPs/APs
	RewriteImplicitAssignmentOPCast RewriteImplicitAssignmentOPCastCallback(&Tool.getReplacements());
	Finder.addMatcher(
		operatorCallExpr(
			allOf(
				hasDescendant(
					implicitCastExpr().bind("castexpr")
				),
				hasDescendant(
					memberExpr().bind("memberexpr")
				),
				hasDescendant(
					implicitCastExpr( isLValueToRValueCast() )
				),
				isUtilityPointer()
			)
		).bind("opercallexpr"),
		&RewriteImplicitAssignmentOPCastCallback);

	// Return statements with implicit conversions to OPs/APs
	RewriteImplicitReturnOPCast RewriteImplicitReturnOPCastCallback(&Tool.getReplacements());
	Finder.addMatcher(
		returnStmt(
			hasDescendant(
				implicitCastExpr(
					allOf(
						hasDescendant(
							bindTemporaryExpr(
								hasDescendant(
									implicitCastExpr( isNonNoopCast() )
								)
							).bind("temporaryexpr")
						),
						isUtilityPointer()
					)
				)
			)
		),
		&RewriteImplicitReturnOPCastCallback);

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

		llvm::errs() << "Output directory: " << outputBaseDir << "\n";
				
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
						llvm::errs() << path << "\n";
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
