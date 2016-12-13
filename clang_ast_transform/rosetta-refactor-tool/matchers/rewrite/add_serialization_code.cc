// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/*
	Find instances where get_self_ptr() or get_self_weak_ptr() is used in a c'tor.
	This is illegal because the weak self-pointer isn't set yet, and will
	result in bad_weak_ptr exception at runtime.
*/

#include "add_serialization_code.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"


#include <string>

using namespace clang;

/*
	Auto add Cereal serialization code to a classes and structs.
	Note: this will fail if with overlapping rewrites, i.e.
	a class definition containing another (private) class definition.
*/


AddSerializationCode::AddSerializationCode(
	tooling::Replacements *Replace,
	RosettaRefactorTool *t) :
	ReplaceMatchCallback(Replace, "AddSerializationCode"),
	refactor_tool(t)
{}

AddSerializationCode::~AddSerializationCode() {}


// Main callback for all matches
void AddSerializationCode::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	sm = Result.SourceManager;
	handle_field_decl(Result);
	handle_record_decl(Result);
}

// Collect class field declarations
void AddSerializationCode::handle_field_decl(const ast_matchers::MatchFinder::MatchResult &Result) {

	const FieldDecl *fielddecl = Result.Nodes.getStmtAs<FieldDecl>("fielddecl");
	if(!rewriteThisFile(fielddecl, *sm))
		return;

	const std::string cls = fielddecl->getParent()->getQualifiedNameAsString();

	if( class_fields.find(cls) == class_fields.end() )
		class_fields[cls] = StringsOP( new Strings() );

	// TODO: don't add raw pointers?

	class_fields[cls]->push_back( fielddecl->getName() );
}

// Collect class record declaration
void AddSerializationCode::handle_record_decl(const ast_matchers::MatchFinder::MatchResult &Result) {

	const CXXRecordDecl *recorddecl = Result.Nodes.getStmtAs<CXXRecordDecl>("recorddecl");
	const Decl *parent = Result.Nodes.getStmtAs<Decl>("parent");
	if(!rewriteThisFile(recorddecl, *sm))
		return;
	if(!recorddecl->isCompleteDefinition())
		return;
	if(!recorddecl->isClass() && !recorddecl->isStruct())
		return;

	Decl_Parent_Pair key = std::make_pair(recorddecl, parent);
	if( std::find(class_recorddecls.begin(), class_recorddecls.end(), key) == class_recorddecls.end() )
		class_recorddecls.push_back( key );
}

// Do the actual rewrite at the end of the translation unit
void AddSerializationCode::onEndOfTranslationUnit() {

	bool has_polymorphic = false;
	std::string suffix;

	// For each class in the translation unit
	for( auto it = class_recorddecls.begin(), end = class_recorddecls.end(); it != end; ++it) {
		const CXXRecordDecl *recorddecl = it->first;
		const Decl *parent = it->second;
		const std::string cls = recorddecl->getQualifiedNameAsString();
		bool is_template = parent && parent->getKind() == Decl::ClassTemplate;

		llvm::outs() << "Processing"
			<< ( is_template ? " template" : "" )
			<< ": " << cls << " (" << recorddecl->getName() << ")\n";

		auto fields = class_fields.find(cls);
		if(fields == class_fields.end())
			continue;

		// Make list of variables to serialize
		std::string class_var_stubs;
		for( auto it2 = fields->second->begin(), end2 = fields->second->end(); it2 != end2; ++it2) {
			llvm::outs() << "\tField: " << *it2 << "\n";
			class_var_stubs += "\t\tCEREAL_NVP( " + *it2 + " ),\n";
		}
		if(class_var_stubs.empty())
			continue;

		class_var_stubs = class_var_stubs.substr(0, class_var_stubs.size()-2) + "\n";

		std::string inline_stub, serialize_code;

		// Add code inline if this is a template or struct
		if(is_template || recorddecl->isStruct()) {
			// Inline
			inline_stub =
				"\ttemplate< class Archive > void serialize( Archive & ar ) {\n"
				"\t\tar(\n" + class_var_stubs + "\t\t);\n"
				"\t}\n";
		} else {
			// Not inline
			inline_stub = "\ttemplate< class Archive > void serialize( Archive & ar );\n";
			serialize_code =
				"template<class Archive>\n"
				"void " + cls + "::serialize( Archive & ar ) {\n"
				"\tar(\n" + class_var_stubs + "\t);\n}\n";
		}

		std::string origCode = getText(*sm, recorddecl);

		// Write new code -- inline stub
		std::string newCode = origCode.substr(0, origCode.size() -1) +
			"\n"
				"#ifdef SERIALIZATION\n"
				"public:\n"
				+ inline_stub +
				"#endif\n"
			"\n}"
			;

		// Write new code -- suffix at end of file, if any
		if(!serialize_code.empty()) {
			std::string guard_define = "DEFINED_" + cls;
			replace(guard_define, "::", "_");

			suffix +=
				"#ifndef " + guard_define + "\n"
				"#define " + guard_define + "\n";

			if(recorddecl->isPolymorphic()) {
				// TODO: To avoid multiple definition linker errors
				// polymorphic class registrations must go into .cc
				// files and not .hh files
				suffix +=
					"// CEREAL_REGISTER_TYPE( " + cls + " );\n\n";
				has_polymorphic = true;
			}

			suffix +=
				serialize_code +
				"#endif // " + guard_define + "\n"
				;
		}

		doRewrite(*sm, recorddecl, origCode, newCode);
	}

	// Write #includes and entire suffix with serialization code
	if(!suffix.empty()) {
		refactor_tool->suffix_ += "\n#ifdef SERIALIZATION\n";
		refactor_tool->suffix_ += has_polymorphic ?
			"#include <cereal/types/polymorphic.hpp>\n" :
			"#include <cereal/cereal.hpp>\n";
		refactor_tool->suffix_ += suffix +
			"#endif // SERIALIZATION\n";
	}

	class_recorddecls.clear();
	class_fields.clear();
}

void
add_serialization_code_rewriter( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements, RosettaRefactorTool * rrt )
{
	using namespace clang::ast_matchers;

	AddSerializationCode *cb = new AddSerializationCode(replacements, rrt );

	finder.addMatcher(
		fieldDecl().bind("fielddecl"),
		cb);

	finder.addMatcher(
		recordDecl(
			hasParent(
				decl().bind("parent")
			)
		).bind("recorddecl"),
		cb);
}
