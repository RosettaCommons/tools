////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace owning_ptr/access_ptr name in templates, etc.

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

////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace owning_ptr/access_ptr name in declarations

class RewritePointerNameDecl : public ReplaceMatchCallback {

public:
	RewritePointerNameDecl(
			tooling::Replacements *Replace,
			const char *tag ="RewritePointerNameDecl") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Decl *node = Result.Nodes.getStmtAs<Decl>("decl");

		if(!rewriteThisFile(node, sm))
			return;

		std::string origCode( getText(sm, node) );
		std::string part1;
		std::string part2;
		size_t p = origCode.find('(');
		if(p != std::string::npos) {
			part1 = std::string(origCode, 0, p);
			part2 = std::string(origCode, p);
		} else {
			part1 = origCode;
			part2 = "";
		}

		replace(part1, "owning_ptr", "shared_ptr");
		replace(part1, "access_ptr", "weak_ptr");

		std::string newCode( part1 + part2 );
		doRewrite(sm, node, origCode, newCode);
	}
};

// CXXMethodDecl in method return type
RewritePointerNameDecl RewritePointerNameDeclCallback(Replacements);
Finder.addMatcher(
	methodDecl().bind("decl"),
	&RewritePointerNameDeclCallback);
