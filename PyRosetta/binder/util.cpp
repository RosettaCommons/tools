// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/util.cpp
/// @brief  Various helper functions
/// @author Sergey Lyskov

#include <util.hpp>

#include <clang/AST/ASTContext.h>
#include <clang/AST/ExprCXX.h>


using namespace llvm;
using namespace clang;
using std::string;
using std::vector;

namespace binder {


/// Split string using given separator
vector<string> split(string const &buffer, string const & separator)
{
	string line;
	vector<string> lines;

	for(uint i=0; i<buffer.size(); ++i) {
		if( buffer.compare(i, separator.size(), separator) ) line.push_back( buffer[i] );
		else {
			lines.push_back(line);
			line.resize(0);
		}
	}

	if( line.size() ) lines.push_back(line);

	return lines;
}


/// Replace all occurrences of string
string replace(string const &s, string const & from, string const &to)
{
	string r{s};
	size_t i = r.size();
	while( ( i = r.rfind(from, i) ) != string::npos) {
		r.replace(i, from.size(), to);
		--i;
	}
	return r;
}


string indent(string const &code, string const &indentation)
{
	auto lines = split(code);
	string r;

	for(auto & l : lines) r += indentation + l + '\n';

	return r;
}


string namespace_from_named_decl(NamedDecl *decl)
{
	string qn = decl->getQualifiedNameAsString();
	string n  = decl->getNameAsString();
	//name = decl->getQualifiedNameAsString();

	int namespace_len = qn.size() - n.size();

	string path = decl->getQualifiedNameAsString().substr(0, namespace_len > 1 ? namespace_len-2 : namespace_len );  // removing trailing '::'

	return path;
}



/// Calculate base (upper) namespace for given one: core::pose::motif --> core::pose
string base_namespace(string const & ns)
{
	size_t f = ns.rfind("::");
	if( f == string::npos ) return "";
	else return ns.substr(0, f);
}


/// Calculate last namespace for given one: core::pose::motif --> motif
string last_namespace(string const & ns)
{
	size_t f = ns.rfind("::");
	if( f == string::npos ) return ns;
	else return ns.substr(f+2, ns.size()-f-2);
}


// replace all _Bool types with bool
void fix_boolean_types(string &type)
{
	string B("_Bool");
	size_t i = 0;
	while( ( i = type.find(B, i) ) != string::npos ) {
		if( ( i==0  or ( !std::isalpha(type[i-1]) and  !std::isdigit(type[i-1]) ) ) and
			( i+B.size() == type.size()  or ( !std::isalpha(type[i+B.size()]) and  !std::isdigit(type[i+B.size()]) ) ) ) type.replace(i, B.size(), "bool");
		++i;
	}
}


// Generate string representation of given expression
string expresion_to_string(clang::Expr *e)
{
	clang::LangOptions lang_opts;
	lang_opts.CPlusPlus = true;
	clang::PrintingPolicy Policy(lang_opts);

	std::string _;
	llvm::raw_string_ostream s(_);
	e->printPretty(s, 0, Policy);
	return s.str();
}


// Generate string representation of given TemplateArgument
string template_argument_to_string(clang::TemplateArgument const &t)
{
	clang::LangOptions lang_opts;
	lang_opts.CPlusPlus = true;
	clang::PrintingPolicy Policy(lang_opts);

	std::string _;
	llvm::raw_string_ostream s(_);
	t.print(Policy, s);
	return s.str();
}




// extract include needed for declaration and add it to includes
bool add_relevant_include(NamedDecl *decl, vector<string> &includes)
{
	ASTContext & ast_context( decl->getASTContext() );
	SourceManager & sm( ast_context.getSourceManager() );

	// outs() << "  Source location: " << F->getLocation().printToString(sm) << "\n";
	//outs() << "  Source location: " << sm.getFilename( decl->getLocation() )  << "\n";

	FileID fid = sm.getFileID( decl->getLocation() );
	SourceLocation include = sm.getIncludeLoc(fid);

	//outs() << "  Include location: " << include.printToString( ast_context.getSourceManager() ) << "\n";
	//outs() << "  printToString: " << include.printToString( ast_context.getSourceManager() ) << "\n";
	//include.dump(ast_context.getSourceManager());
	//outs() << "  Line: " << sm.getCharacterData(include) << "\n";
	// outs() << "  Source location: " << sm.getFileEntryForID( FullSourceLoc(F->getLocation(), sm ).getFileID() )->getName() << "\n";

	string include_string;
	if( include.isValid() ) {
		char const *data = sm.getCharacterData(include);

		char terminator = *data == '"' ? '"' : '>';

		for(; *data and  *data != terminator; ++data ) include_string.push_back(*data);
		if(*data == terminator) include_string.push_back(*data);

		includes.push_back(include_string);
		//outs() << "  Include for " << decl->getNameAsString() << ": " << include_string << "\n";
	}

	return include.isValid();
}




} // namespace binder
