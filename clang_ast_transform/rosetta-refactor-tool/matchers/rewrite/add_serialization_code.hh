/*
	Auto add Cereal serialization code to a classes and structs.
*/


class AddSerializationCode : public ReplaceMatchCallback {

typedef std::vector< std::string > Strings;
typedef std::shared_ptr< Strings > StringsOP;
typedef std::map< std::string, StringsOP > ClassFields;

ClassFields class_fields;
std::vector<const CXXRecordDecl *> class_recorddecls;
SourceManager *sm;

RosettaRefactorTool *refactor_tool;

public:
	AddSerializationCode(
		tooling::Replacements *Replace,
		RosettaRefactorTool *t) :
		ReplaceMatchCallback(Replace, "AddSerializationCode"),
		refactor_tool(t)
		{}

		class_fields[cls]->push_back( fielddecl->getName() );
	}

	// Collect class record declaration
	void handle_record_decl(const ast_matchers::MatchFinder::MatchResult &Result) {

		const CXXRecordDecl *recorddecl = Result.Nodes.getStmtAs<CXXRecordDecl>("recorddecl");
		if(!rewriteThisFile(recorddecl, *sm))
			return;

		if(!recorddecl->isClass())
			return;

		if( std::find(class_recorddecls.begin(), class_recorddecls.end(), recorddecl) == class_recorddecls.end() )
			class_recorddecls.push_back( recorddecl );
	}

	// Do the actual rewrite at the end of the translation unit
	virtual void onEndOfTranslationUnit() {

		bool cereal_includes = false;

		for( auto it = class_recorddecls.begin(), end = class_recorddecls.end(); it != end; ++it) {
			const CXXRecordDecl *recorddecl = *it;
			const std::string cls = recorddecl->getQualifiedNameAsString();

			llvm::outs() << "Processing class: " << cls << " (" << recorddecl->getName() << ")\n";

			auto fields = class_fields.find(cls);
			if(fields == class_fields.end())
				continue;

			std::string class_var_stubs;
			for( auto it2 = fields->second->begin(), end2 = fields->second->end(); it2 != end2; ++it2) {
				llvm::outs() << "\tField: " << *it2 << "\n";
				class_var_stubs += "\t\tCEREAL_NVP( " + *it2 + " ),\n";
			}
			if(class_var_stubs.empty())
				continue;

			class_var_stubs = class_var_stubs.substr(0, class_var_stubs.size()-2) + "\n";
			std::string serialize_code =
				"template<class Archive>\n"
				"void " + cls + "::serialize( Archive & ar ) {\n"
				"\tar(\n" + class_var_stubs + "\t);\n}\n";

			// Write new code
			std::string origCode = getText(*sm, recorddecl);
			std::string newCode = origCode.substr(0, origCode.size() -1) +
				"\n\n"
					"#ifdef SERIALIZATION\n"
					"public:\n"
					"\ttemplate< class Archive > void serialize( Archive & ar );\n"
					"#endif\n"
				"\n}"
				;

			std::string guard_define = "DEFINED_" + cls;
			replace(guard_define, "::", "_");

			refactor_tool->suffix_ +=
				"\n\n"
				"#ifdef SERIALIZATION\n"
				"#ifndef " + guard_define + "\n"
				"#define " + guard_define + "\n";

			if(!cereal_includes) {
				refactor_tool->suffix_ += recorddecl->isPolymorphic() ?
					"#include <cereal/types/polymorphic.hpp>\n" :
					"#include <cereal/cereal.hpp>\n";
				cereal_includes = true;
			}

			if(recorddecl->isPolymorphic())
			refactor_tool->suffix_ +=
				"CEREAL_REGISTER_TYPE( " + cls + " );\n\n";

			refactor_tool->suffix_ +=
				serialize_code +
				"#endif // " + guard_define + "\n"
				"#endif // SERIALIZATION\n"
				;

			doRewrite(*sm, recorddecl, origCode, newCode);
		}

		class_recorddecls.clear();
		class_fields.clear();
	}

};

AddSerializationCode *cb = new AddSerializationCode(Replacements, this);

Finder.addMatcher(
	fieldDecl().bind("fielddecl"),
	cb);

Finder.addMatcher(
	recordDecl().bind("recorddecl"),
	cb);
