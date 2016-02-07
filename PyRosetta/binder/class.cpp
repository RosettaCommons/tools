// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/class.cpp
/// @brief  Binding generation for C++ struct and class objects
/// @author Sergey Lyskov

#include <class.hpp>
#include <function.hpp>
#include <util.hpp>

#include <fmt/format.h>

//#include <vector>

#include <clang/AST/DeclTemplate.h>

using namespace llvm;
using namespace clang;

using std::string;
//using std::vector;
//using std::unordered_map;

using namespace fmt::literals;

namespace binder {

// generate classtemplate specialization for ClassTemplateSpecializationDecl or empty string otherwise
string template_specialization(clang::CXXRecordDecl *C)
{
	string templ;

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		templ += "<";
		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
			//outs() << " template argument: " << template_argument_to_string(t->getTemplateArgs()[i]) << "\n";
			templ += template_argument_to_string(t->getTemplateArgs()[i]) + ",";

			//if( t->getTemplateArgs()[i].ArgKind() == TemplateArgument::ArgKind::Integral ) outs() << " template arg:" << t->getTemplateArgs()[i].<< "\n";
			//outs() << expresion_to_string( t->getTemplateArgs()[i].getAsExpr() ) << "\n";
		}
		templ.back() = '>';
	}

	return templ;
}


// generate class name that could be used in bindings code indcluding template specialization if any
string class_name(clang::CXXRecordDecl *C)
{
	return C->getNameAsString() + template_specialization(C);
}


// generate qualified class name that could be used in bindings code indcluding template specialization if any
string class_qualified_name(clang::CXXRecordDecl *C)
{
	return C->getQualifiedNameAsString() + template_specialization(C);
}


// Generate bindings for class data member
string bind_data_member(FieldDecl *d, string const &class_qualified_name)
{
	return ".def_readwrite(\"{}\", &{}::{})"_format(d->getNameAsString(), class_qualified_name, d->getNameAsString());
}


/// check if generator can create binding
bool is_bindable(clang::CXXRecordDecl *C)
{
	bool r = true;

	r &= C->isCompleteDefinition()  and  !C->isDependentType(); //C->hasDefinition();

	// todo: bindging for abstract classes
	//if(r) r &= !C->isAbstract();  // need an 'if' here or clang assert got triggered on classed with incomplete definitions

	// if(r) {
	// 	outs() << C->getQualifiedNameAsString() << " isCXXClassMember:" << C->isCXXClassMember() << " isCompleteDefinition:" << C->isCompleteDefinition() //<< " isBeingDefined:" << C->isBeingDefined()
	// 		   << " isDependentType:" << C->isDependentType() << " isCXXInstanceMember:" << C->isCXXInstanceMember() << " isExternallyVisible:" << C->isExternallyVisible()
	// 		   << " getNumTemplateParameterLists:" << C->getNumTemplateParameterLists()  << " isa<ClassTemplateSpecializationDecl>:" << isa<ClassTemplateSpecializationDecl>(C) << "\n";
	// 	C->dump();
	// 	//auto l = C->getTemplateParameterList();
	// 	//for(auto p = l.begin(); p!=l.end(); ++p) outs() << (*p)->getNameAsString() << "\n";
	// }

	return r;
}


/// generate binding code
string ClassBinder::operator()(string const &module_variable_name, string const &indentation) const
{
	string c;

	//class_<A>(module_a, "A")
	c += R"(pybind11::class_<{}>({}, "{}"))"_format(class_qualified_name(C), module_variable_name, class_name(C)) + '\n';

	if( !C->isAbstract() ) {
		if( C->ctor_begin() == C->ctor_end() ) {  // No constructors defined, adding default constructor
			c+= "\t.def(pybind11::init<>())\n\n";
		}
		else {
			bool added=false;
			for(auto t = C->ctor_begin(); t != C->ctor_end(); ++t) {
				if( t->getAccess() == AS_public  and  !t->isMoveConstructor() ) { added=true;  c+= "\t.def(pybind11::init<{}>())\n"_format( function_arguments(*t) ); }
			}
			if(added) c += '\n';
		}
	}

	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(FieldDecl *f = dyn_cast<FieldDecl>(*d) ) {
			if( f->getAccess() == AS_public ) c+= '\t' + bind_data_member(f, class_qualified_name(C)) + '\n';
		}
	}

	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( is_bindable(*m)  and  m->getAccess() == AS_public  and   !isa<CXXConstructorDecl>(*m)  and   !isa<CXXDestructorDecl>(*m)  /*and  !m->isStatic()*/) {
			//(*m)->dump();
			c += '\t' + bind_function(*m) + '\n';
		}
	}

	c += ";\n\n";

	return indent(c, indentation);
}


} // namespace binder
