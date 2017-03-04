// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/*
	Find all function calls in class methods
	To get a parsable list, run:
	./run.sh test-cases.cc -matchers=find_calls -colors=false|grep ' => '|sed 's/ => /\t/'
*/

#include "bad_pointer_casts.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include "../../ast_matchers.hh"

#include <string>

using namespace clang;



/*
	Code quality checker finder:
	- Find all cases where "this" is being put into an OP or AP
	- Find all cases where reference is being put into an OP or AP

	Example:

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

class X;
typedef utility::pointer::access_ptr< class X > XAP;
typedef utility::pointer::owning_ptr< class X > XOP;

class X : public utility::pointer::ReferenceCount {
private:
	XAP partner_;
public:
	X() {}
	void setPartner(XAP p) { partner_ = p; }

	void foo1() {
		XOP x1 = new X;         // OK
		X & x1_r = *x1;         // OK
		X * xp1 = x1.get();     // OK
		X * xp2 = &x1_r;        // OK
		XOP x2 = &x1_r;         // BAD
		XAP x3 = &x1_r;         // BAD
		x1->setPartner( this ); // BAD
		setPartner( x1 );       // OK
		setPartner( &x1_r );    // BAD
		setPartner( xp1 );      // OK
	}
};

*/

BadPointerCastFinder::BadPointerCastFinder(
	tooling::Replacements *Replace,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag),
	verbose_( false )
{}

BadPointerCastFinder::BadPointerCastFinder(
	clang::tooling::Replacements * Replace,
	bool verbose,
	const char *tag
) :
	ReplaceMatchCallback(Replace, tag),
	verbose_( verbose )
{}


void
BadPointerCastFinder::run( ast_matchers::MatchFinder::MatchResult const & Result ) {
	SourceManager &sm = *Result.SourceManager;
	const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
	const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");

	if(!rewriteThisFile(castFrom, sm))
		return;

	const std::string locStr = castFrom->getSourceRange().getBegin().printToString(sm);
	const std::string castFromType = QualType::getAsString( castFrom->getType().split() );
	const std::string castToType = QualType::getAsString( castTo->getType().split() );
	const std::string castFromTypeD = QualType::getAsString( castFrom->getType().getSplitDesugaredType() );
	const std::string castToTypeD = QualType::getAsString( castTo->getType().getSplitDesugaredType() );

	if (verbose_) {
		llvm::errs() << tag() << ": " << locStr << "\n";

		llvm::errs() << "castFrom: " << color("green") << castFromType << color("") << "\n";
		if(castFromType != castFromTypeD && !castFromTypeD.empty())
			llvm::errs() << "          " << color("green") << castFromTypeD << color("") << "\n";

		llvm::errs() << "castTo:   " << color("red") << castToType << color("") << "\n";
		if(castToType != castToTypeD && !castToTypeD.empty())
			llvm::errs() << "          " << color("red") << castToTypeD << color("") << "\n";
	}

	llvm::outs() << tag() << "\t" << locStr << "\t" << castFromType << "\t" << castToType << "\n";
}


void
add_bad_pointer_cast_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;

	finder.addMatcher(
		constructExpr(
			has(
				thisExpr().bind("castFrom")
			),
			isUtilityPointer()
		).bind("castTo"),
		new BadPointerCastFinder(replacements, "BadPointerCastFinder:thisExpr"));


	// XOP x2 = &x1_r;
	// XAP x3 = &x1_r;
	// setPartner( &x1_r );

	finder.addMatcher(
		unaryOperator(
			hasParent(
				constructExpr( isUtilityPointer() ).bind("castTo")
			),
			has(
				expr().bind("castFrom")
			)
		),
		new BadPointerCastFinder(replacements, "BadPointerCastFinder:unaryOperator"));
}
