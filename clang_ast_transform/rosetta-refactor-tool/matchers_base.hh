////////////////////////////////////////////////////////////////////////////////////////////////////
// Base classes

class ReplaceMatchCallback : public ast_matchers::MatchFinder::MatchCallback {
public:
	ReplaceMatchCallback(tooling::Replacements *Replace)
			: Replace(Replace) {}

private:
	tooling::Replacements *Replace;

protected:
	template <typename T>
	void doRewrite(
		const std::string & tag,
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
};
