/*
	Find instances where get_self_ptr() or get_self_weak_ptr() is used in a c'tor.
	This is illegal because the weak self-pointer isn't set yet, and will
	result in bad_weak_ptr exception at runtime.
*/

class FieldDeclFinder : public ReplaceMatchCallback {
public:
	FieldDeclFinder(
		tooling::Replacements *Replace,
		const char *tag = "FieldDeclFinder") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const FieldDecl *fielddecl = Result.Nodes.getStmtAs<FieldDecl>("fielddecl");

		if(!rewriteThisFile(fielddecl, sm))
			return;

		const std::string locStr( fielddecl->getSourceRange().getBegin().printToString(sm) );
		std::string code = getText(sm, fielddecl);
		const std::string type(
			stripQualifiers(
				QualType::getAsString( fielddecl->getType().getSplitDesugaredType() )
			)
		);		
		const std::string ns = fielddecl->getQualifiedNameAsString();
		
		replace(code, "\n", " ");
		llvm::outs() << locStr << "\t" << ns << "\t" << type << "\t" << code << "\n";
	}
};


Finder.addMatcher(
	fieldDecl().bind("fielddecl"),
	new FieldDeclFinder(Replacements));
