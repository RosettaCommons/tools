////////////////////////////////////////////////////////////////////////////////////////////////////
// Base classes

class ReplaceMatchCallback : public ast_matchers::MatchFinder::MatchCallback {
public:
	ReplaceMatchCallback(tooling::Replacements *Replace, const char *tag)
			: Replace(Replace), tag(tag) {}
	
private:
	tooling::Replacements *Replace;

protected:
	std::string tag;

	template <typename T>
	void doRewrite(
		SourceManager &sm, T * node,
		const std::string & origCode,
		const std::string & newCode
	) {
		if(origCode == newCode)
			return;

		if(!checkAndDumpRewrite(tag, sm, node, newCode))
			return;
		Replace->insert(Replacement(sm, node, newCode));
	}

	template <typename T>
	bool rewriteThisFile(T * node, SourceManager & sm) {
		if(!node)
			return false;
		const FullSourceLoc FullLocation = FullSourceLoc(node->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID()) {
			// llvm::errs() << tag << " skipping file: " << node->getSourceRange().getBegin().printToString(sm) << "\n";
			return false;
		}

		if(Debug)
			llvm::errs() << color("blue") << tag << " matched: " << node->getSourceRange().getBegin().printToString(sm) << color("") << "\n";

		return true;
	}

};
