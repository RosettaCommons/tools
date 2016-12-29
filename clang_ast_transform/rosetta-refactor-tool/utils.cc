// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include "utils.hh"

#include <string>
#include <cstring>
#include <set>

#include "clang/Basic/SourceManager.h"
#include "clang/Lex/Lexer.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// String utils

// luki created an anonymous namespace and I don't know why
// namespace {

std::string color(const char *_color, bool show_colors) {

	if(!show_colors)
		return "";

	std::string color(_color);
	if(color == "")       return "\033[0m";
	if(color == "black")  return "\033[0m";
	if(color == "red")    return "\033[31m";
	if(color == "green")  return "\033[32m";
	if(color == "brown")  return "\033[33m";
	if(color == "blue")   return "\033[34m";
	if(color == "purple") return "\033[35m";
	if(color == "cyan")   return "\033[36m";
	if(color == "gray")   return "\033[37m";
	return "";
}

void replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
}

std::string trim(const std::string & s, const char *whitespace ) {
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
		type == "utility::pointer::access_ptr" ||
		type == "utility::pointer::shared_ptr" ||
		type == "utility::pointer::weak_ptr" ||
		type == "std::shared_ptr" ||
		type == "std::weak_ptr";
}

// Is it our utility access pointer that we care about?
bool checkIsUtilityAccessPointer(const std::string& fulltype) {
	std::string type = extractContainerType(fulltype);
	return
		type == "utility::pointer::access_ptr" ||
		type == "utility::pointer::weak_ptr" ||
		type == "std::weak_ptr";
}

// Is it our utility owning pointer that we care about?
bool checkIsUtilityOwningPointer(const std::string& fulltype) {
	std::string type = extractContainerType(fulltype);
	return
		type == "utility::pointer::owning_ptr" ||
		type == "utility::pointer::shared_ptr" ||
		type == "std::shared_ptr";
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

// end luki's anonymous namespace } // namespace


////////////////////////////////////////////////////////////////////////////////////////////////////
// Clang Utils

// namespace {

// GRRR
// static std::set<std::string> RewrittenLocations_;

// Returns the text that makes up 'node' in the source.
// Returns an empty string if the text cannot be found.
std::string getText(
	const clang::SourceManager &SourceManager,
	const clang::SourceLocation &StartSpellingLocation,
	const clang::SourceLocation &EndSpellingLocation
) {
  if (!StartSpellingLocation.isValid() || !EndSpellingLocation.isValid()) {
    //llvm::errs() << "getText: invalid locations\n";
    return std::string();
  }
  bool Invalid = true;
  const char *Text =
      SourceManager.getCharacterData(StartSpellingLocation, &Invalid);
  if (Invalid) {
    //lvm::errs() << "getText: can't get character data\n";
    return std::string();
  }
  std::pair<clang::FileID, unsigned> Start =
		SourceManager.getDecomposedLoc(StartSpellingLocation);
  std::pair<clang::FileID, unsigned> End =
		SourceManager.getDecomposedLoc(clang::Lexer::getLocForEndOfToken(
		EndSpellingLocation, 0, SourceManager, clang::LangOptions()));
  if (Start.first != End.first) {
    // Start and end are in different files.
    //llvm::errs() << "getText: Start/end in different files\n";
    return std::string();
  }
  if (End.second < Start.second) {
    // Shuffling text with macros may cause this.
    //llvm::errs() << "getText: end before start\n";
    return std::string();
  }
  return std::string(Text, End.second - Start.second);
}


// } // anon namespace


