// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/*
	Find all function calls in class methods
	To get a parsable list, run:
	./run.sh test-cases.cc -matchers=find_calls -colors=false|grep ' => '|sed 's/ => /\t/'
*/

#include "calls.hh"

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

CallsFinder::CallsFinder(
	clang::tooling::Replacements *Replace,
	const char *tag ) :
	ReplaceMatchCallback(Replace, tag) {}

void CallsFinder::run(const clang::ast_matchers::MatchFinder::MatchResult &Result) {
	using namespace clang;
	SourceManager &sm = *Result.SourceManager;
	const CXXMethodDecl *caller = Result.Nodes.getStmtAs<CXXMethodDecl>("caller");
	const CallExpr *call = Result.Nodes.getStmtAs<CallExpr>("call");
	const CXXMemberCallExpr *memberCall = Result.Nodes.getStmtAs<CXXMemberCallExpr>("call");
	const DeclRefExpr *declref = Result.Nodes.getStmtAs<DeclRefExpr>("declref");

	if(!rewriteThisFile(caller, sm))
		return;

	if(!caller->hasBody())
		return;

	// TODO: retrieve enclosing namespace of caller and call node -- how?

	if ( ! call->getDirectCallee() ) return;
	const std::string callingMethodName = caller->getNameAsString();
	const std::string calledMethodName = call->getDirectCallee()->getNameAsString();

	const std::string callingMethodClassName = caller->getParent()->getNameAsString();
	std::string calledMethodClassName;
	if(memberCall) {
		calledMethodClassName = memberCall->getRecordDecl()->getNameAsString();
	} else if(declref) {
		calledMethodClassName = declref->getFoundDecl()->getQualifiedNameAsString();
		size_t p = calledMethodClassName.find_last_of(':');
		if(p != std::string::npos && p > 0)
			calledMethodClassName = calledMethodClassName.substr(0, p-1);
	}

	const std::string locStr( call->getSourceRange().getBegin().printToString(sm) );
	const std::string callCode = getText(sm, call);

	llvm::outs()
		<< "@ " << locStr << /* color("cyan") << " (" << tag << ")" << color("") << */ "\n"
		<< color("red") << callCode << color("") << "\n"
		<< color("purple") << callingMethodClassName << color("") << "::"
		<< color("purple") << callingMethodName << color("") << " => "
		<< color("green") << calledMethodClassName << color("") << "::"
		<< color("green") << calledMethodName << color("") << "\n"
		;
}

// Member calls: object->method()
/*
CXXMethodDecl 0x43a0370 </data/rosetta/tools/clang_ast_transform/test-cases.cc:209:2, line:212:2> line:209:7 caller 'void (ClassBOP)'
|-ParmVarDecl 0x43a02f0 <col:14, col:23> col:23 b 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
|-CompoundStmt 0x440b688 <col:26, line:212:2>
| |-CXXMemberCallExpr 0x440b548 <line:210:3, col:12> 'void'
| | `-MemberExpr 0x440b518 <col:3, col:6> '<bound member function type>' ->new_a 0x439f1e0
| |   `-CXXOperatorCallExpr 0x440b4d8 <col:3> 'pointer':'class ClassB *'
| |     |-ImplicitCastExpr 0x440b4c0 <col:4> 'pointer (*)(void) const' <FunctionToPointerDecay>
| |     | `-DeclRefExpr 0x440b498 <col:4> 'pointer (void) const' lvalue CXXMethod 0x42df6e0 'operator->' 'pointer (void) const'
| |     `-ImplicitCastExpr 0x440b480 <col:3> 'const class utility::pointer::owning_ptr<class ClassB>' lvalue <NoOp>
| |       `-DeclRefExpr 0x440b458 <col:3> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>' lvalue ParmVar 0x43a02f0 'b' 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
| `-CXXMemberCallExpr 0x440b660 <line:211:3, col:17> 'void'
|   `-MemberExpr 0x440b630 <col:3, col:6> '<bound member function type>' ->set_a_null 0x439eba0
|     `-CXXOperatorCallExpr 0x440b5f0 <col:3> 'pointer':'class ClassB *'
|       |-ImplicitCastExpr 0x440b5d8 <col:4> 'pointer (*)(void) const' <FunctionToPointerDecay>
|       | `-DeclRefExpr 0x440b5b0 <col:4> 'pointer (void) const' lvalue CXXMethod 0x42df6e0 'operator->' 'pointer (void) const'
|       `-ImplicitCastExpr 0x440b598 <col:3> 'const class utility::pointer::owning_ptr<class ClassB>' lvalue <NoOp>
|         `-DeclRefExpr 0x440b570 <col:3> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>' lvalue ParmVar 0x43a02f0 'b' 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
*/

void
add_calls_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang;
	using namespace clang::ast_matchers;

	finder.addMatcher(
		methodDecl(
			forEachDescendant(
				memberCallExpr().bind("call")
			)
		).bind("caller"),
		new CallsFinder( replacements, "memberCallExpr"));

	finder.addMatcher(
		methodDecl(
			forEachDescendant(
				callExpr(
					has(
						declRefExpr().bind("declref")
					)
				).bind("call")
			)
		).bind("caller"),
		new CallsFinder(replacements, "callExpr"));

}


// (static) Calls: object::method()

/*
CXXMethodDecl 0x380a590 </local/luki/tools/clang_ast_transform/test-cases.cc:213:2, line:220:2> line:213:7 caller 'void (ClassBOP)'
| |-ExprWithCleanups 0x29a0518 <line:219:3, col:19> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
| | `-CXXBindTemporaryExpr 0x29a04f8 <col:3, col:19> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>' (CXXTemporary 0x29a04f0)
| |   `-CallExpr 0x29a04c0 <col:3, col:19> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
| |     `-ImplicitCastExpr 0x29a04a8 <col:3, col:11> 'ClassBOP (*)(void)' <FunctionToPointerDecay>
| |       `-DeclRefExpr 0x29a0418 <col:3, col:11> 'ClassBOP (void)' lvalue CXXMethod 0x292cdd0 'factory' 'ClassBOP (void)'
| `-DeclStmt 0x29a33a0 <line:220:3, col:36>
|   `-VarDecl 0x29a0540 <col:3, col:35> col:12 my_b 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
|     `-ExprWithCleanups 0x29a3388 <col:19, col:35> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
|       `-CXXConstructExpr 0x29a3350 <col:19, col:35> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>' 'void (const class utility::pointer::owning_ptr<class ClassB> &)' elidable
|         `-MaterializeTemporaryExpr 0x29a3330 <col:19, col:35> 'const class utility::pointer::owning_ptr<class ClassB>' lvalue
|           `-ImplicitCastExpr 0x29a3318 <col:19, col:35> 'const class utility::pointer::owning_ptr<class ClassB>' <NoOp>
|             `-CXXBindTemporaryExpr 0x29a0658 <col:19, col:35> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>' (CXXTemporary 0x29a0650)
|               `-CallExpr 0x29a0620 <col:19, col:35> 'ClassBOP':'class utility::pointer::owning_ptr<class ClassB>'
|                 `-ImplicitCastExpr 0x29a0608 <col:19, col:27> 'ClassBOP (*)(void)' <FunctionToPointerDecay>
|                   `-DeclRefExpr 0x29a05d0 <col:19, col:27> 'ClassBOP (void)' lvalue CXXMethod 0x292cdd0 'factory' 'ClassBOP (void)'
*/
