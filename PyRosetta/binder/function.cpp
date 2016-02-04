// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/function.cpp
/// @brief  Binding generation for static and member functions
/// @author Sergey Lyskov

#include <function.hpp>


#include <util.hpp>
#include <fmt/format.h>

#include <clang/AST/DeclCXX.h>
#include <clang/AST/ExprCXX.h>
#include <clang/AST/ASTContext.h>

#include <vector>


using namespace llvm;
using namespace clang;

using std::string;
using std::vector;
using std::unordered_map;

using namespace fmt::literals;

namespace binder {


// Generate function argument list separate by comma: int, bool, std::sting
string function_arguments(clang::FunctionDecl *record)
{
	string r;

	for(uint i=0; i<record->getNumParams(); ++i) {
		r += record->getParamDecl(i)->getOriginalType().getCanonicalType().getAsString();
		if( i != record->getNumParams()-1 ) r += ", ";
	}

	fix_boolean_types(r);

	return r;
}


// Generate function pointer type string for given function: void (*)(int, doule)_ or  void (ClassName::*)(int, doule)_ for memeber function
string function_pointer_type(FunctionDecl *F)
{
	string r;
	string prefix { F->isCXXClassMember() ? cast<CXXRecordDecl>( F->getParent() )->getQualifiedNameAsString() + "::" : "" };

	r += F->getReturnType().getAsString();  r+= " ({}*)("_format(prefix);

	r += function_arguments(F);

	r += ") ";

	fix_boolean_types(r);

	return r;
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


// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
string bind_function(FunctionDecl *F)
{
	string function_name { F->getNameAsString() };
	string function_qualified_name { F->getQualifiedNameAsString() };

	string r = R"(.def("{}", ({}) &{}, "doc")"_format(function_name, function_pointer_type(F), function_qualified_name);

	//F->dump();

	for(auto p = F->param_begin(); p != F->param_end(); ++p) {
		string defalt_argument = (*p)->hasDefaultArg() ? " = " + expresion_to_string( (*p)->getDefaultArg() ) : "";
		r += ", pybind11::arg(\"{}\"){}"_format( string( (*p)->getName() ), defalt_argument );

		//add_relevant_include(*p, includes);

		//outs() << (*p)->getDefaultArg()->getAsString();
		//(*p)->dump();
		// if( (*p)->hasDefaultArg() ) {
		// 	outs() << "  CXXDefaultArgExpr: " << (*p)->getName() << " = " << expresion_to_string( (*p)->getDefaultArg() ) << "\n";

		// 	//SourceRange s = (*p)->getDefaultArgRange();
		// 	//s.getBegin().dump( F->getASTContext().getSourceManager() );


		// 	//if( auto e = dyn_cast<clang::CXXDefaultArgExpr>( (*p)->getDefaultArg() ) ) {
		// 		//e->dump();
		// 	// }

		// 	// Expr::EvalResult result;
		// 	// if( (*p)->getDefaultArg()->EvaluateAsRValue(result, F->getASTContext() ) ) {
		// 	// 	outs() << "  Default for: " << (*p)->getName() << " = " << result.Val.getAsString(F->getASTContext(), (*p)->getOriginalType() ) << "\n";
		// 	// }
		// }

	}

	r += ')';

	return r;
}


/// check if generator can create binding
bool FunctionBinder::is_bindable() const
{
	return true;
}


/// check if generator can create binding
bool is_bindable(FunctionDecl *F)
{
	QualType ret( F->getReturnType() );

	//if( )

	return true;
}


/// generate binding code
string FunctionBinder::operator()(string const &module_variable_name, string const &indentation) const
{
	//return indentation+"// Function: " + F->getQualifiedNameAsString() + "\n";
	string c = indentation + module_variable_name + bind_function(F) + ";\n\n";

	return c;
}


} // namespace binder
