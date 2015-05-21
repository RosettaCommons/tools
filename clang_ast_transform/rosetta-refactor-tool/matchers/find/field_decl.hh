/*
	Find instances where get_self_ptr() or get_self_weak_ptr() is used in a c'tor.
	This is illegal because the weak self-pointer isn't set yet, and will
	result in bad_weak_ptr exception at runtime.
*/

class FieldDeclFinder : public ReplaceMatchCallback {

public:
	FieldDeclFinder(tooling::Replacements *Replace) :
		ReplaceMatchCallback(Replace, "FieldDeclFinder")
		{}

	// Main callback for all matches
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {

		SourceManager *sm = Result.SourceManager;
		const FieldDecl *decl = Result.Nodes.getStmtAs<FieldDecl>("fielddecl");
		if(!rewriteThisFile(decl, *sm))
			return;

		const std::string type(
			QualType::getAsString( decl->getType().split() )
		);
		const std::string typeD(
			QualType::getAsString( decl->getType().getSplitDesugaredType() )
		);

		const std::string name = decl->getQualifiedNameAsString();
		const std::string cls = decl->getParent()->getQualifiedNameAsString();
		const std::string loc = decl->getSourceRange().getBegin().printToString(*sm);

		const CharSourceRange range = CharSourceRange::getTokenRange(decl->getSourceRange());
		SourceLocation SpellingBegin = sm->getSpellingLoc(range.getBegin());
		SourceLocation SpellingEnd = sm->getSpellingLoc(range.getEnd());
		std::pair<FileID, unsigned> Start = sm->getDecomposedLoc(SpellingBegin);
		std::pair<FileID, unsigned> End = sm->getDecomposedLoc(SpellingEnd);

		llvm::outs()
			<< "field" << "\t"
			<< name << "\t"
			<< cls << "\t"
			<< loc << "\t"
			<< Start.second << "-" << End.second << "\t"
			<< type << "\t"
			<< typeD << "\t"
			;
		llvm::outs() << "\n";
	}

};


FieldDeclFinder *cb = new FieldDeclFinder(Replacements);

Finder.addMatcher(
	fieldDecl().bind("fielddecl"),
	cb);
