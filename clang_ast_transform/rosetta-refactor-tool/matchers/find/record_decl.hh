/*
	Find CXX records (class, struct, ...)
*/

class RecordDeclFinder : public ReplaceMatchCallback {

public:
	RecordDeclFinder(tooling::Replacements *Replace) :
		ReplaceMatchCallback(Replace, "RecordDeclFinder")
		{}

	// Main callback for all matches
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {

		SourceManager *sm = Result.SourceManager;
		const CXXRecordDecl *decl = Result.Nodes.getStmtAs<CXXRecordDecl>("recorddecl");
		const Decl *parent = Result.Nodes.getStmtAs<Decl>("parent");
		if(!rewriteThisFile(decl, *sm))
			return;
		if(!decl->isCompleteDefinition())
			return;

		const std::string name = decl->getQualifiedNameAsString();
		const std::string loc = decl->getSourceRange().getBegin().printToString(*sm);

		const CharSourceRange range = CharSourceRange::getTokenRange(decl->getSourceRange());
		SourceLocation SpellingBegin = sm->getSpellingLoc(range.getBegin());
		SourceLocation SpellingEnd = sm->getSpellingLoc(range.getEnd());
		std::pair<FileID, unsigned> Start = sm->getDecomposedLoc(SpellingBegin);
		std::pair<FileID, unsigned> End = sm->getDecomposedLoc(SpellingEnd);

		llvm::outs()
			<< decl->getKindName() << "\t"
			<< name << "\t"
			<< "" << "\t"
			<< loc << "\t"
			<< Start.second << "-" << End.second << "\t"
			<< (parent && parent->getKind() == Decl::ClassTemplate) << "\t"
			<< decl->getAccess() << "\t"
			<< decl->isPolymorphic() << "\t"
			;
		llvm::outs() << "\n";
		if ( decl->getNumBases() > 0 ) {
		  for ( CXXRecordDecl::base_class_const_iterator iter = decl->bases_begin(), iter_end = decl->bases_end(); iter != iter_end; ++iter ) {
		    llvm::outs()
		      << "parent" << "\t"
		      << name << "\t"
		      << "" << "\t"
		      << loc << "\t"
		      << sm->getDecomposedLoc( iter->getLocStart() ).second << "-" << sm->getDecomposedLoc( iter->getLocEnd() ).second << "\t"
		      //<< iter->getType().getAsString() << "\t"
		      //<< iter->getType().getUnqualifiedType().getAsString() << "\t"
		      << iter->getType().getCanonicalType().getAsString() << "\t"
		      << iter->getAccessSpecifier() << "\t"
		      << iter->isVirtual();
		    llvm::outs() << "\n";
		  }
		}
	}
};

Finder.addMatcher(
	recordDecl(
		hasParent(
			decl().bind("parent")
		)
	).bind("recorddecl"),
	new RecordDeclFinder(Replacements));
	
