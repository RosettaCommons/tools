/*
	Replace owning_ptr/access_ptr name in templates, etc.
*/

class RewritePointerName : public ReplaceMatchCallback {
	
public:
	RewritePointerName(
		tooling::Replacements *Replace,
		const char *tag ="RewritePointerName") :
		ReplaceMatchCallback(Replace, tag) {}
	
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Stmt *node = Result.Nodes.getStmtAs<Stmt>("stmt");

		if(!rewriteThisFile(node, sm))
			return;
		
		std::string origCode( getText(sm, node) );
		std::string newCode( origCode );

		replace(newCode, "owning_ptr", "shared_ptr");
		replace(newCode, "access_ptr", "weak_ptr");

		doRewrite(sm, node, origCode, newCode);
	}
};


// CXXUnresolvedConstructExpr in templates
RewritePointerName RewritePointerNameCallback(Replacements);
Finder.addMatcher(
	unresolvedConstructExpr( isUtilityPointer() ).bind("stmt"),
	&RewritePointerNameCallback);

