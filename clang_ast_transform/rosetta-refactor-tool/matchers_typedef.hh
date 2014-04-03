////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace owning_ptr/access_ptr in typedefs
//   typedef utility::pointer::access_ptr< Some > SomeAP
//   typedef utility::pointer::owning_ptr< Some > SomeOP
//   typedef utility::pointer::access_ptr< Some const > SomeCAP
//   typedef utility::pointer::owning_ptr< Some const > SomeCOP
// OK

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

		doRewrite("RewriteTypedefDecl", sm, decl, origCode, newCode);
	}
};


// Typedefs for access_ptr and owning_ptr
RewriteTypedefDecl RewriteTypedefDeclCallback(Replacements);
Finder.addMatcher(
	decl( isTypedefDecl() ).bind("typedefdecl"),
	&RewriteTypedefDeclCallback);
