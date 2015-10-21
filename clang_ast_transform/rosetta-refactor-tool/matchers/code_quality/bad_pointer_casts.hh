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
		tooling::Replacements *Replace,
		const char *tag = "BadPointerCastFinder") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
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

		if(Verbose) {
			llvm::errs() << tag << ": " << locStr << "\n";

			llvm::errs() << "castFrom: " << color("green") << castFromType << color("") << "\n";
			if(castFromType != castFromTypeD && !castFromTypeD.empty())
				llvm::errs() << "          " << color("green") << castFromTypeD << color("") << "\n";

			llvm::errs() << "castTo:   " << color("red") << castToType << color("") << "\n";
			if(castToType != castToTypeD && !castToTypeD.empty())
				llvm::errs() << "          " << color("red") << castToTypeD << color("") << "\n";
		}

		llvm::outs() << tag << "\t" << locStr << "\t" << castFromType << "\t" << castToType << "\n";
	}
};

// x1->setPartner( this );

Finder.addMatcher(
	cxxConstructExpr(
		has(
			cxxThisExpr().bind("castFrom")
		),
		isUtilityPointer()
	).bind("castTo"),
	new BadPointerCastFinder(Replacements, "BadPointerCastFinder:thisExpr"));


// XOP x2 = &x1_r;
// XAP x3 = &x1_r;
// setPartner( &x1_r );

Finder.addMatcher(
	unaryOperator(
		hasParent(
			cxxConstructExpr( isUtilityPointer() ).bind("castTo")
		),
		has(
			expr().bind("castFrom")
		)
	),
	new BadPointerCastFinder(Replacements, "BadPointerCastFinder:unaryOperator"));
