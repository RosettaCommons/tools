/*
	Find all function calls in class methods
	To get a parsable list, run:
	./run.sh test-cases.cc -matchers=find_calls -colors=false|grep ' => '|sed 's/ => /\t/'
*/

class CallsFinder : public ReplaceMatchCallback {
public:
	CallsFinder(
		tooling::Replacements *Replace,
		const char *tag = "CallsFinder") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const CXXMethodDecl *caller = Result.Nodes.getStmtAs<CXXMethodDecl>("caller");
		const CXXMemberCallExpr *call = Result.Nodes.getStmtAs<CXXMemberCallExpr>("call");

		if(!rewriteThisFile(caller, sm))
			return;

		if(!caller->hasBody())
			return;

		const std::string callingMethodClassName = caller->getParent()->getNameAsString();
		const std::string callingMethodName = caller->getNameAsString();
		const std::string calledMethodClassName = call->getRecordDecl()->getNameAsString();
		const std::string calledMethodName = call->getDirectCallee()->getNameAsString();

		const std::string locStr( call->getSourceRange().getBegin().printToString(sm) );
		const std::string callCode = getText(sm, call);

		llvm::outs()
			<< "@ " << locStr << /* color("cyan") << " (" << tag << ")" << color("") << */ "\n"
			<< color("red") << callCode << color("") << "\n"
			<< color("purple") << callingMethodClassName << color("") << ":"
			<< color("purple") << callingMethodName << color("") << " => "
			<< color("green") << calledMethodClassName << color("") << ":"
			<< color("green") << calledMethodName << color("") << "\n"
			;
	}
};

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

Finder.addMatcher(
	methodDecl(
		forEachDescendant(
			memberCallExpr().bind("call")
		)
	).bind("caller"),
	new CallsFinder(Replacements));
