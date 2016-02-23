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

#include <class.hpp>
#include <util.hpp>
#include <fmt/format.h>

#include <clang/AST/DeclCXX.h>
#include <clang/AST/ASTContext.h>

#include <clang/AST/ExprCXX.h>


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
	string prefix, maybe_const;
	if( auto m = dyn_cast<CXXMethodDecl>(F) ) {
		prefix = m->isStatic() ? "" : class_qualified_name( cast<CXXRecordDecl>( F->getParent() ) ) + "::";
	    maybe_const = m->isConst() ? " const" : "";
	}

	r += F->getReturnType().getCanonicalType().getAsString();  r+= " ({}*)("_format(prefix);

	r += function_arguments(F);

	r += ")" + maybe_const;

	fix_boolean_types(r);

	return r;
}


// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
string bind_function(FunctionDecl *F)
{
	string function_name { F->getNameAsString() };
	string function_qualified_name { F->getQualifiedNameAsString() };

	CXXMethodDecl * m = dyn_cast<CXXMethodDecl>(F);
	string maybe_static = m and m->isStatic() ? "_static" : "";

	string r = R"(.def{}("{}", ({}) &{}, "doc")"_format(maybe_static, function_name, function_pointer_type(F), function_qualified_name);

	for(auto p = F->param_begin(); p != F->param_end(); ++p) {
		// .def("myFunction", py::arg("arg") = SomeType(123));
		//string defalt_argument = (*p)->hasDefaultArg() ? " = ({})({})"_format( (*p)->getOriginalType().getCanonicalType().getAsString(), expresion_to_string( (*p)->getDefaultArg() ) ) : "";
		//string defalt_argument = (*p)->hasDefaultArg() ? expresion_to_string( (*p)->getDefaultArg() ) : "";
		//r += ", pybind11::arg(\"{}\"){}"_format( string( (*p)->getName() ), defalt_argument );

		// .def("myFunction", py::arg_t<int>("arg", 123, "123"));
		if( (*p)->hasDefaultArg()  and  !(*p)->hasUninstantiatedDefaultArg() ) {
			string arg_type = (*p)->getOriginalType().getCanonicalType().getAsString();  fix_boolean_types(arg_type);
			r += ", pybind11::arg_t<{}>(\"{}\", {}, \"{}\")"_format(arg_type,
																	string( (*p)->getName() ),
																	expresion_to_string( (*p)->getDefaultArg() ),
																	replace(expresion_to_string( (*p)->getDefaultArg() ), "\"", "\\\"") );

			// clang::Expr *e = (*p)->getDefaultArg();
			// clang::Expr::EvalResult result;
			// outs() << "ex: " << expresion_to_string(e) << " EvaluateAsRValue: " << e->EvaluateAsRValue(result, F->getASTContext());
			// outs() << " res: " << result.Val.getAsString(F->getASTContext(), (*p)->getOriginalType()) << "\n";
		}
		else r += ", pybind11::arg(\"{}\")"_format( string( (*p)->getName() ) );

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
// bool FunctionBinder::is_bindable() const
// {
// 	return true;
// }


bool is_bindable(QualType const &qt)
{
	bool r = true;

	r &= !qt->isFunctionPointerType()  and  !qt->isInstantiationDependentType()  and  !qt->isArrayType();  //and  !qt->isConstantArrayType()  and  !qt->isIncompleteArrayType()  and  !qt->isVariableArrayType()  and  !qt->isDependentSizedArrayType()

	if( PointerType const *pt = dyn_cast<PointerType>( qt.getTypePtr() ) ) {
		if( pt->getPointeeType()->isPointerType() ) return false;  // refuse to bind 'value**...' types
		if( pt->getPointeeType()->isArrayType() or pt->getPointeeType()->isConstantArrayType() ) return false;  // refuse to bind 'T* v[]...' types
		//qt->dump();
		r &= is_bindable( pt->getPointeeType() );
	}

	if( ReferenceType const *rt = dyn_cast<ReferenceType>( qt.getTypePtr() ) ) {
		//rt->dump();
		//outs() << "Ref " << qt.getAsString() << " -> " << is_bindable( rt->getPointeeType() ) << "\n";
		r &= is_bindable( rt->getPointeeType() );
	}

	if( Type const *tp = qt.getTypePtrOrNull() ) {
		if( CXXRecordDecl *rd = tp->getAsCXXRecordDecl() ) r &= is_bindable(rd);
	}


	// if( auto rt = dyn_cast<RecordType>( qt.getTypePtr() ) ) {
	// 	if( auto cxxr = dyn_cast<CXXRecordDecl>( rt->getDecl() ) ) r &= is_bindable(cxxr);
	// }

	//outs() << " isIncompleteArrayType(): " << qt->isIncompleteArrayType() << " r:" << r << "\n";
	//qt->dump();

	//outs() << "Qt " << qt.getAsString() << " -> " << r << " isInstantiationDependentType: " << qt->isInstantiationDependentType() << "\n";
	return r;
}


/// check if generator can create binding
bool is_bindable(FunctionDecl const *F)
{
	bool r = true;

	// todo: bindging for operators and type conversion
	r &= !F->isOverloadedOperator()  and  !isa<CXXConversionDecl>(F)  and  !F->isDeleted();

	QualType rt( F->getReturnType() );

	r &= is_bindable(rt);

	for(auto p = F->param_begin(); p != F->param_end(); ++p) r &= is_bindable( (*p)->getOriginalType() );

	return r;
}


/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
string FunctionBinder::id() const
{
	string maybe_const;
	if( auto m = dyn_cast<CXXMethodDecl>(F) ) maybe_const = m->isConst() ? " const" : "";

	string r = F->getReturnType().getCanonicalType().getAsString() + " "+ F->getQualifiedNameAsString() + "(" + function_arguments(F) + ")" + maybe_const;
	fix_boolean_types(r);
	return r;
}


bool FunctionBinder::bindable() const
{
	return binder::is_bindable(F);
}


/// generate binding code for this object and all its dependencies
void FunctionBinder::bind(Context &context)
{
	if( is_binded() ) return;

	string const indentation="\t";
	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(F) );
	string const include = relevant_include(F);

	code()  = indentation+"// " + F->getQualifiedNameAsString() + " " + function_arguments(F) + " file:" + include.substr(1, include.size()-2) + " line:" + line_number(F) + "\n";
	code() += indentation + module_variable_name + bind_function(F) + ";\n\n";
}


} // namespace binder
