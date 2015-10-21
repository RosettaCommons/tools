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

class ObjOnStackFinder : public ReplaceMatchCallback {
public:
	ObjOnStackFinder(
		tooling::Replacements *Replace,
		const char *tag = "ObjOnStackFinder") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
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
		
		if(Verbose) {
			llvm::errs() << tag << ": " << locStr << "\n";
			llvm::errs() << "vardecl: " << color("red") << vardeclType << color("");
			if(vardeclType != vardeclTypeD && !vardeclTypeD.empty())
				llvm::errs() << " :: " << color("red") << vardeclTypeD << color("");
			llvm::errs() << ": " << origCode << "\n";
		}

		llvm::outs() << tag << "\t" << locStr << "\t" << vardeclType << "\t" << origCode << "\n";
	}
};

/*
|-DeclStmt 0xaac2d48 <line:210:7, col:28>
| `-VarDecl 0xaac2cc0 <col:7, col:24> col:24 pose 'core::pose::Pose':'class core::pose::Pose'
|   `-CXXConstructExpr 0xaac2d18 <col:24> 'core::pose::Pose':'class core::pose::Pose' 'void (void)'
*/

// Pose
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::pose::Pose")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Pose"));

// Movers
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::moves::Mover")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Mover"));

// Packer
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::graph::Graph")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Graph"));
		
// Residues
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::chemical::ResidueType")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::ResidueType"));

Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::chemical::ResidueTypeSet")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::ResidueTypeSet"));

Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::conformation::Residue")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Residue"));

// Conformation
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::conformation::Conformation")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Conformation"));

// Datacache
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("basic::datacache::CacheableData")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::CacheableData"));

// Scoring
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::scoring::ScoreFunction")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::ScoreFunction"));

Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::scoring::hbonds::HBondSet")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::HBondSet"));

// Graph
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::pack::task::PackerTask")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::PackerTask"));

// Atom Tree
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::kinematics::tree::Atom")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Atom"));

// I/O
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::io::silent::SilentStruct")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::SilentStruct"));

// Tag
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("utility::tag::Tag")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Tag"));

// Environment
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::environment::Environment")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Environment"));

// Topology Broker
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::topology_broker::TopologyBroker")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::TopologyBroker"));
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::topology_broker::TopologyClaimer")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::TopologyClaimer"));

// Fragments
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("core::fragment::FragData")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::FragData"));

Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::frag_picker::FragmentPicker")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::FragmentPicker"));
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::frag_picker::VallChunk")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::VallChunk"));
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::frag_picker::VallProvider")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::VallProvider"));

Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::abinitio::Templates")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::Templates"));

// EnzDes
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::toolbox::match_enzdes_util::EnzConstraintIO")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::EnzConstraintIO"));
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::toolbox::match_enzdes_util::EnzConstraintParameters")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::EnzConstraintParameters"));
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::toolbox::match_enzdes_util::InvrotTreeNodeBase")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::InvrotTreeNodeBase"));

// Misc Protocols
Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::flexpack::rotamer_set::FlexbbRotamerSets")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::FlexbbRotamerSets"));

Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::match::MatcherTask")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::MatcherTask"));

Finder.addMatcher(
	varDecl(hasType(cxxRecordDecl(isSameOrDerivedFrom("protocols::optimize_weights::OptEMultifunc")))).bind("vardecl"),
	new ObjOnStackFinder(Replacements, "ObjOnStackFinder::OptEMultifunc"));
