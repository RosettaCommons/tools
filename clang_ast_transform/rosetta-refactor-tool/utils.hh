////////////////////////////////////////////////////////////////////////////////////////////////////
// String utils

namespace {

void replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
}

std::string trim(const std::string & s, const char *whitespace =0) {
	size_t start, end;
	if(!whitespace)
		whitespace = " \r\n\t";
	for(start = 0; start < s.length() && strchr(whitespace, s[start]); start++);
	for(end = s.length() -1; end > 0  && strchr(whitespace, s[end]);   end--);
	return std::string(s, start, end - start + 1);
}

bool endsWith(const std::string& a, const std::string& b) {
	if (b.size() > a.size()) return false;
	return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}

bool beginsWith(const std::string& a, const std::string& b) {
	return (a.compare(0, b.length(), b) == 0);
}

bool checkIsUtilityPointer(std::string const & type) {
	if(beginsWith(type, "const ")) // trim const
		return checkIsUtilityPointer(std::string(type, 6));

	return
		beginsWith(type, "class utility::pointer::owning_ptr") ||
		beginsWith(type, "class utility::pointer::access_ptr");
}

bool checkContainsUtilityPointer(std::string const & type) {
	if(
		beginsWith(type, "utility::vector0<") ||
		beginsWith(type, "utility::vector1<")
	) {
		size_t p = type.find('<');
		if(p == std::string::npos)
			return false;
		std::string sub_type = trim(std::string(type, p), "<> ");
		return checkIsUtilityPointer(sub_type);
	}
	return false;
}

} // namespace


////////////////////////////////////////////////////////////////////////////////////////////////////
// Clang Utils

namespace {
	
std::set<std::string> RewrittenLocations_;
	
// Returns the text that makes up 'node' in the source.
// Returns an empty string if the text cannot be found.
static std::string getText(
	const SourceManager &SourceManager,
	const SourceLocation &StartSpellingLocation, 
	const SourceLocation &EndSpellingLocation
) {
  if (!StartSpellingLocation.isValid() || !EndSpellingLocation.isValid()) {
		llvm::errs() << "getText: invalid locations\n";
    return std::string();
  }
  bool Invalid = true;
  const char *Text =
      SourceManager.getCharacterData(StartSpellingLocation, &Invalid);
  if (Invalid) {
		llvm::errs() << "getText: can't get character data\n";
    return std::string();
  }
  std::pair<FileID, unsigned> Start =
      SourceManager.getDecomposedLoc(StartSpellingLocation);
  std::pair<FileID, unsigned> End =
      SourceManager.getDecomposedLoc(Lexer::getLocForEndOfToken(
          EndSpellingLocation, 0, SourceManager, LangOptions()));
  if (Start.first != End.first) {
    // Start and end are in different files.
		llvm::errs() << "getText: Start/end in different files\n";
    return std::string();
  }
  if (End.second < Start.second) {
    // Shuffling text with macros may cause this.
		llvm::errs() << "getText: end before start\n";
    return std::string();
  }
  return std::string(Text, End.second - Start.second);
}

template <typename T>
static std::string getText(const SourceManager &SourceManager, const T *Node) {
	if(!Node)
		return std::string();
	return getText(SourceManager,
		SourceManager.getSpellingLoc(Node->getLocStart()), 
		SourceManager.getSpellingLoc(Node->getLocEnd()));
}

template <typename T, typename U>
static std::string getTextToDelim(const SourceManager &SourceManager, const T *Node, const U *DelimNode) {
	if(!Node || !DelimNode)
		return std::string();
	return getText(SourceManager,
		SourceManager.getSpellingLoc(Node->getLocStart()), 
		SourceManager.getSpellingLoc(DelimNode->getLocStart()));
}

template <typename T>
bool checkAndDumpRewrite(
	const std::string & tag,
	SourceManager & sm, T * node,
	const std::string & newCodeStr
) {

	const std::string locStr( node->getSourceRange().getBegin().printToString(sm) );
	bool notYetSeen = (RewrittenLocations_.find(locStr) == RewrittenLocations_.end());
	if(notYetSeen)
		RewrittenLocations_.insert(locStr);
	
	if(!verbose)
		return notYetSeen;

	const std::string origCodeStr = getText(sm, node);
		
	llvm::errs() 
		<< "@ " << locStr << " \033[36m(" << tag << ")\033[0m" 
		<< (!notYetSeen ? " \033[31m[skipped]\033[0m" : "") <<	"\n" 
		<< "- \033[31m" << origCodeStr << "\033[0m\n"
		<< "+ \033[32m" << newCodeStr << "\033[0m\n"
		<< "\n";
		
	return notYetSeen;
}

} // anon namespace
