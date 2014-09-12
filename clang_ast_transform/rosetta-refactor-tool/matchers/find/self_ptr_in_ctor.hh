/*
	Find instances where get_self_ptr() or get_self_weak_ptr() is used in a c'tor.
	This is illegal because the weak self-pointer isn't set yet, and will
	result in bad_weak_ptr exception at runtime.
*/

class SelfPtrInCtorFinder : public ReplaceMatchCallback {
public:
	SelfPtrInCtorFinder(
		tooling::Replacements *Replace,
		const char *tag = "SelfPtrInCtorFinder") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const CXXConstructorDecl *ctor = Result.Nodes.getStmtAs<CXXConstructorDecl>("ctor");

		if(!rewriteThisFile(ctor, sm))
			return;

		if(!ctor->hasBody())
			return;

		// TODO: retrieve enclosing namespace of caller and call node -- how?
		// DeclRefExpr::getFoundDecl()->getQualifiedNameAsString() 
		
		const std::string locStr( ctor->getSourceRange().getBegin().printToString(sm) );
		const std::string code = getText(sm, ctor);
		
		if(
			code.find("get_self_ptr()") == std::string::npos &&
			code.find("get_self_weak_ptr()") == std::string::npos
		)
			return;
		
		llvm::outs()
			<< "@ " << locStr << /* color("cyan") << " (" << tag << ")" << color("") << */ "\n"
			;
	}
};


Finder.addMatcher(
	constructorDecl(
		// Doesn't work... so ugly "solution" above
		//forEachDescendant(
		//	memberCallExpr(
		//		hasName("get_self_ptr")
		//	)
		//)
	).bind("ctor"),
	new SelfPtrInCtorFinder(Replacements));
