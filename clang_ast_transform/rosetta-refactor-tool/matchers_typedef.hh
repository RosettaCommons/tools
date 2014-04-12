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

		if(!rewriteThisFile(node, sm))
			return;
		
		const std::string origCode( getText(sm, node) );
		if(!checkIsUtilityPointer(origCode))
			return;

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

// Parameter declaration in methods
RewriteTypedefDecl RewriteTypedefDeclCallback2(Replacements, "RewriteTypedefDecl:paramVarDecl");
Finder.addMatcher(
	parmVarDecl().bind("decl"),
	&RewriteTypedefDeclCallback2);

// Field (variable) declaration in classes
RewriteTypedefDecl RewriteTypedefDeclCallback3(Replacements, "RewriteTypedefDecl:fieldDecl");
Finder.addMatcher(
	fieldDecl().bind("decl"),
	&RewriteTypedefDeclCallback3);

// Method return type declatation
RewriteTypedefDecl RewriteTypedefDeclCallback4(Replacements, "RewriteTypedefDecl:methodDecl");
Finder.addMatcher(
	methodDecl().bind("decl"),
	&RewriteTypedefDeclCallback4);

// Local method varialble declaration
RewriteTypedefDecl RewriteTypedefDeclCallback5(Replacements, "RewriteTypedefDecl:varDecl");
Finder.addMatcher(
	varDecl().bind("decl"),
	&RewriteTypedefDeclCallback5);
