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
#ifndef INCLUDED_matchers_code_quality_bad_pointer_casts_HH
#define INCLUDED_matchers_code_quality_bad_pointer_casts_HH

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Tooling/Refactoring.h"

#include "../../matchers_base.hh"

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

class BadPointerCastFinder : public ReplaceMatchCallback {
public:
	BadPointerCastFinder(
		clang::tooling::Replacements *Replace,
		const char *tag = "BadPointerCastFinder" );

	BadPointerCastFinder(
		clang::tooling::Replacements *Replace,
		bool verbose,
		const char *tag = "BadPointerCastFinder" );


	virtual void run(const clang::ast_matchers::MatchFinder::MatchResult &Result);

private:
	bool verbose_;
};

void
add_bad_pointer_cast_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements );

#endif
