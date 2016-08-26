// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_utils_HH
#define INCLUDED_utils_HH

#include <string>
#include <cstring>
#include <set>

#include "clang/Basic/SourceManager.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// String utils

//namespace {

std::string color(const char *_color, bool show_colors=false);

void replace(std::string& str, const std::string& from, const std::string& to);

std::string trim(const std::string & s, const char *whitespace =0);

bool endsWith(const std::string& a, const std::string& b);

bool beginsWith(const std::string& a, const std::string& b);

////////////////////////////////////////////////////////////////////////////////////////////////////
// C++ code parsing string utils

// Container types from which we extract the contained type
bool checkContainerType(const std::string& container_type);

// Qualifiers we don't care about and drop
std::string stripQualifiers(const std::string& type);

// Extract container type from code, i.e. vector, map, etc.
std::string extractContainerType(const std::string& container);

// Extract contained type in a vector, map, etc.
std::string extractContainedType(const std::string& container);

// Is it our utility pointer that we care about?
bool checkIsUtilityPointer(const std::string& fulltype);

// Is it our utility access pointer that we care about?
bool checkIsUtilityAccessPointer(const std::string& fulltype);

// Is it our utility owning pointer that we care about?
bool checkIsUtilityOwningPointer(const std::string& fulltype);

// Does it contain a pointer that we care about?
bool checkContainsUtilityPointer(const std::string& fulltype);

bool checkIsClassOperator(const std::string& type);

//} // namespace


////////////////////////////////////////////////////////////////////////////////////////////////////
// Clang Utils

//namespace {

static std::set<std::string> RewrittenLocations_;

// Returns the text that makes up 'node' in the source.
// Returns an empty string if the text cannot be found.
std::string getText(
	const clang::SourceManager &SourceManager,
	const clang::SourceLocation &StartSpellingLocation,
	const clang::SourceLocation &EndSpellingLocation
);

template <typename T>
std::string getText(const clang::SourceManager &SourceManager, const T *Node);

template <typename T, typename U>
std::string getTextToDelim(const clang::SourceManager &SourceManager, const T *Node, const U *DelimNode);

template <typename T>
bool checkAndMarkSourceLocation(
	T * node,
	clang::SourceManager & sm
);


template <typename T>
bool checkAndDumpRewrite(
	const std::string & tag,
	clang::SourceManager & sm, T * node,
	const std::string & newCodeStr,
	bool verbose = false
);

template <typename T>
std::string getText(const clang::SourceManager &SourceManager, const T *Node) {
	if(!Node)
		return std::string();
	return getText(SourceManager,
		SourceManager.getSpellingLoc(Node->getLocStart()),
		SourceManager.getSpellingLoc(Node->getLocEnd()));
}

template <typename T, typename U>
std::string getTextToDelim(const clang::SourceManager &SourceManager, const T *Node, const U *DelimNode) {
	if(!Node || !DelimNode)
		return std::string();
	return getText(SourceManager,
		SourceManager.getSpellingLoc(Node->getLocStart()),
		SourceManager.getSpellingLoc(DelimNode->getLocStart()));
}

template <typename T>
bool checkAndMarkSourceLocation(
	T * node,
	clang::SourceManager & sm
) {
	const std::string locStr( node->getSourceRange().getBegin().printToString(sm) );
	std::string locStrPartial( locStr );
	size_t p = locStrPartial.find_last_of(':');
	if(p != std::string::npos)
		locStrPartial = std::string(locStrPartial, 0, p);

	bool notYetSeen = (RewrittenLocations_.find(locStrPartial) == RewrittenLocations_.end());
	if(notYetSeen)
		RewrittenLocations_.insert(locStrPartial);
	return notYetSeen;
}

template <typename T>
bool checkAndDumpRewrite(
	const std::string & tag,
	clang::SourceManager & sm, T * node,
	const std::string & newCodeStr,
	bool verbose
) {

	bool notYetSeen = checkAndMarkSourceLocation(node, sm);
	if(!verbose)
		return notYetSeen;

	const std::string locStr( node->getSourceRange().getBegin().printToString(sm) );
	const std::string origCodeStr = getText(sm, node);

	llvm::errs() << "@ " << locStr << color("cyan") << " (" << tag << ")" << color("");
	if(!notYetSeen)
		llvm::errs()  << " " << color("red") << "[skipped]" << color("");
	llvm::errs() <<	"\n"
		<< "- " << color("red") << origCodeStr << color("") << "\n"
		<< "+ " << color("green") << newCodeStr << color("") << "\n"
		<< "\n";

	return notYetSeen;
}

//} // anon namespace

#endif
