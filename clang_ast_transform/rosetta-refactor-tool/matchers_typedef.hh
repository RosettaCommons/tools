////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace owning_ptr/access_ptr in typedefs
//   typedef utility::pointer::access_ptr< Some > SomeAP
//   typedef utility::pointer::owning_ptr< Some > SomeOP
//   typedef utility::pointer::access_ptr< Some const > SomeCAP
//   typedef utility::pointer::owning_ptr< Some const > SomeCOP
// OK

class RewriteTypedefDecl : public ReplaceMatchCallback {
	
public:
	RewriteTypedefDecl(
			tooling::Replacements *Replace,
			const char *tag ="RewriteTypedefDecl") :
		ReplaceMatchCallback(Replace, tag) {}
	
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Decl *node = Result.Nodes.getStmtAs<Decl>("decl");
		const FullSourceLoc FullLocation = FullSourceLoc(node->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;
		
		std::string origCode( getText(sm, node) );
		std::string newCode( origCode );

		replace(newCode, "owning_ptr", "shared_ptr");
		replace(newCode, "access_ptr", "weak_ptr");

		doRewrite(sm, node, origCode, newCode);
	}
};


// Typedefs for access_ptr and owning_ptr
RewriteTypedefDecl RewriteTypedefDeclCallback1(Replacements, "RewriteTypedefDecl:isTypedefDecl");
Finder.addMatcher(
	decl( isTypedefDecl() ).bind("decl"),
	&RewriteTypedefDeclCallback1);

// ParmVarDecl in templates and function defs
RewriteTypedefDecl RewriteTypedefDeclCallback2(Replacements, "RewriteTypedefDecl:paramVarDecl");
Finder.addMatcher(
	parmVarDecl().bind("decl"),
	&RewriteTypedefDeclCallback2);
