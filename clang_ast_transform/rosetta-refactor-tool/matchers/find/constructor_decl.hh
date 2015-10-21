/*
	Find constructors
*/

class ConstructorDeclFinder : public ReplaceMatchCallback {

public:
	ConstructorDeclFinder(tooling::Replacements *Replace) :
		ReplaceMatchCallback(Replace, "ConstructorDeclFinder")
		{}

	// Main callback for all matches
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {

		SourceManager *sm = Result.SourceManager;
		const CXXConstructorDecl *decl = Result.Nodes.getStmtAs<CXXConstructorDecl>("constructordecl");
		if(!rewriteThisFile(decl, *sm))
			return;

		const std::string name = decl->getQualifiedNameAsString();
		const std::string cls = decl->getParent()->getQualifiedNameAsString();
		const std::string loc = decl->getSourceRange().getBegin().printToString(*sm);

		const CharSourceRange range = CharSourceRange::getTokenRange(decl->getSourceRange());
		SourceLocation SpellingBegin = sm->getSpellingLoc(range.getBegin());
		SourceLocation SpellingEnd = sm->getSpellingLoc(range.getEnd());
		std::pair<FileID, unsigned> Start = sm->getDecomposedLoc(SpellingBegin);
		std::pair<FileID, unsigned> End = sm->getDecomposedLoc(SpellingEnd);

		llvm::outs()
			<< "constructor" << "\t"
			<< name << "\t"
			<< cls << "\t"
			<< loc << "\t"
			<< Start.second << "-" << End.second << "\t"
			<< decl->getAccess() << "\t"
			<< decl->isDefaultConstructor() << "\t"
			<< decl->isCopyConstructor() << "\t"
			<< decl->getNumCtorInitializers() << "\t"
			<< decl->getMinRequiredArguments() << "\t"
			;
		llvm::outs() << "\n";
	}
};

Finder.addMatcher(
	cxxConstructorDecl().bind("constructordecl"),
	new ConstructorDeclFinder(Replacements));
