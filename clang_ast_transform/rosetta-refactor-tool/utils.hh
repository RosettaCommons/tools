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

////////////////////////////////////////////////////////////////////////////////////////////////////
// C++ code parsing string utils

// Container types from which we extract the contained type
bool checkContainerType(const std::string& container_type) {
	return
		container_type == "utility::vector0" ||
		container_type == "utility::vector1" ||
		container_type == "std::vector" ||
		container_type == "std::list" ||
		container_type == "std::set" ||
		container_type == "std::map" ||
		container_type == "std::deque";
}

// Qualifiers we don't care about and drop
std::string stripQualifiers(const std::string& type) {
	if(beginsWith(type, "const ")) // trim const (typedefs)
		return stripQualifiers(std::string(type, strlen("const ")));
	if(beginsWith(type, "class ")) // trim class (templates)
		return stripQualifiers(std::string(type, strlen("class ")));
	if(beginsWith(type, "static ")) // trim static in var decls
		return stripQualifiers(std::string(type, strlen("static ")));
	return trim(type);
}

// Extract container type from code, i.e. vector, map, etc.
std::string extractContainerType(const std::string& container) {
	// Strip utility::vector1<utility::pointer::ownining_ptr<class ClassX>, allocator ... >
	size_t p = container.find('<');
	if(p == std::string::npos)
		return "";
		
	std::string container_type = stripQualifiers(std::string(container, 0, p));
	//llvm::errs() << "container_type: " << container_type << "\n";
	return container_type;
}

// Extract contained type in a vector, map, etc.
std::string extractContainedType(const std::string& container) {
	
	// Strip utility::vector1<utility::pointer::ownining_ptr<class ClassX>, allocator ... >
	size_t p = container.find('<');
	size_t q = container.find_last_of('>');
	if(p == std::string::npos)
		return "";
		
	std::string container_type(container, 0, p);	
	std::string contained_type = std::string(container, p + 1, q - p - 1); 
	
	if(container_type == "std::map") {
		// Handle maps differently: map key, mapped type, allocator -- we want mapped type
		q = contained_type.find(',');
		if(q != std::string::npos)
			contained_type = trim(std::string(contained_type, q+1));
	}
	
	// Use first def (drop allocator, etc)
	q = contained_type.find(',');
	if(q != std::string::npos)
		contained_type = trim(std::string(contained_type, 0, q));

	contained_type = stripQualifiers(contained_type);
	
	//llvm::errs() << "contained_type: " << contained_type << "\n";
	return contained_type;
}

// Is it our utility pointer that we care about?
bool checkIsUtilityPointer(const std::string& fulltype) {

	std::string type = extractContainerType(fulltype);
	return
		type == "utility::pointer::owning_ptr" ||
		type == "utility::pointer::access_ptr";
}

// Is it our utility access pointer that we care about?
bool checkIsUtilityAccessPointer(const std::string& fulltype) {
	std::string type = extractContainerType(fulltype);
	return type == "utility::pointer::access_ptr";
}

// Is it our utility owning pointer that we care about?
bool checkIsUtilityOwningPointer(const std::string& fulltype) {
	std::string type = extractContainerType(fulltype);
	return type == "utility::pointer::owning_ptr";
}

// Does it contain a pointer that we care about?
bool checkContainsUtilityPointer(const std::string& fulltype) {

	std::string type = stripQualifiers(fulltype);
	std::string container_type = extractContainerType(fulltype);

	if(container_type.empty())
		return false;
		
	if(checkContainerType(container_type))
		return checkIsUtilityPointer(extractContainedType(type));

	return false;
}

bool checkIsClassOperator(const std::string& type) {
	return
		type == "operator()" ||
		type == "operator=";
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
	std::string locStrPartial( locStr );
	size_t p = locStrPartial.find_last_of(':');
	if(p != std::string::npos)
		locStrPartial = std::string(locStrPartial, 0, p);
		
	bool notYetSeen = (RewrittenLocations_.find(locStrPartial) == RewrittenLocations_.end());
	if(notYetSeen)
		RewrittenLocations_.insert(locStrPartial);
	
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
