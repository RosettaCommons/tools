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

#include "obj_on_stack.hh"

#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/Basic/SourceManager.h"

#include <string>

using namespace clang;

/*
	Code quality checker finder:
	- Find locations where an object that uses enable_shared_from_this<>
	  is declared on the stack; see matchers at end of this file
	  for a list of classes as of Sept 2014.

	Example:

class X : public std::enable_shared_from_this< X > {
public:
	X() {}
	~X() {}
};

void foo() {
	X x;          // BAD
	X & x_r = x;  // OK
}

*/

ObjOnStackFinder::ObjOnStackFinder(
	tooling::Replacements *Replace,
	const char *tag )
:
	ReplaceMatchCallback(Replace, tag),
	verbose_( false )
{}

ObjOnStackFinder::ObjOnStackFinder(
	tooling::Replacements *Replace,
	bool verbose,
	const char *tag )
:
	ReplaceMatchCallback(Replace, tag),
	verbose_( verbose )
{}

void ObjOnStackFinder::run(const ast_matchers::MatchFinder::MatchResult &Result) {
	SourceManager &sm = *Result.SourceManager;
	const VarDecl *vardecl = Result.Nodes.getStmtAs<VarDecl>("vardecl");

	if(!rewriteThisFile(vardecl, sm))
		return;

	const std::string locStr = vardecl->getSourceRange().getBegin().printToString(sm);
	const std::string vardeclType = QualType::getAsString( vardecl->getType().split() );
	const std::string vardeclTypeD = QualType::getAsString( vardecl->getType().getSplitDesugaredType() );

	if(checkIsUtilityPointer(vardeclTypeD)) // OPs/APs OK
		return;
	if(endsWith(vardeclTypeD, "&") || endsWith(vardeclTypeD, "*")) // Refs and Ptr* OK
		return;

	const std::string origCode = getText(sm, vardecl);

	if(verbose_) {
		llvm::errs() << tag() << ": " << locStr << "\n";
		llvm::errs() << "vardecl: " << color("red") << vardeclType << color("");
		if(vardeclType != vardeclTypeD && !vardeclTypeD.empty())
			llvm::errs() << " :: " << color("red") << vardeclTypeD << color("");
		llvm::errs() << ": " << origCode << "\n";
	}

	llvm::outs() << tag() << "\t" << locStr << "\t" << vardeclType << "\t" << origCode << "\n";
}

void
add_obj_on_stack_finder( clang::ast_matchers::MatchFinder & finder, clang::tooling::Replacements * replacements )
{
	using namespace clang::ast_matchers;

	/*
	|-DeclStmt 0xaac2d48 <line:210:7, col:28>
	| `-VarDecl 0xaac2cc0 <col:7, col:24> col:24 pose 'core::pose::Pose':'class core::pose::Pose'
	|   `-CXXConstructExpr 0xaac2d18 <col:24> 'core::pose::Pose':'class core::pose::Pose' 'void (void)'
	*/

	// Pose
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::pose::Pose")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Pose"));

	// Movers
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::moves::Mover")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Mover"));

	// Packer
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::graph::Graph")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Graph"));

	// Residues
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::chemical::ResidueType")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::ResidueType"));

	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::chemical::ResidueTypeSet")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::ResidueTypeSet"));

	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::conformation::Residue")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Residue"));

	// Conformation
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::conformation::Conformation")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Conformation"));

	// Datacache
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("basic::datacache::CacheableData")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::CacheableData"));

	// Scoring
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::scoring::ScoreFunction")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::ScoreFunction"));

	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::scoring::hbonds::HBondSet")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::HBondSet"));

	// Graph
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::pack::task::PackerTask")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::PackerTask"));

	// Atom Tree
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::kinematics::tree::Atom")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Atom"));

	// I/O
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::io::silent::SilentStruct")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::SilentStruct"));

	// Tag
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("utility::tag::Tag")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Tag"));

	// Environment
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::environment::Environment")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Environment"));

	// Topology Broker
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::topology_broker::TopologyBroker")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::TopologyBroker"));
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::topology_broker::TopologyClaimer")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::TopologyClaimer"));

	// Fragments
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("core::fragment::FragData")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::FragData"));

	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::frag_picker::FragmentPicker")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::FragmentPicker"));
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::frag_picker::VallChunk")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::VallChunk"));
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::frag_picker::VallProvider")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::VallProvider"));

	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::abinitio::Templates")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::Templates"));

	// EnzDes
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::toolbox::match_enzdes_util::EnzConstraintIO")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::EnzConstraintIO"));
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::toolbox::match_enzdes_util::EnzConstraintParameters")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::EnzConstraintParameters"));
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::toolbox::match_enzdes_util::InvrotTreeNodeBase")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::InvrotTreeNodeBase"));

	// Misc Protocols
	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::flexpack::rotamer_set::FlexbbRotamerSets")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::FlexbbRotamerSets"));

	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::match::MatcherTask")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::MatcherTask"));

	finder.addMatcher(
		varDecl(hasType(recordDecl(isSameOrDerivedFrom("protocols::optimize_weights::OptEMultifunc")))).bind("vardecl"),
		new ObjOnStackFinder(replacements, "ObjOnStackFinder::OptEMultifunc"));

}
