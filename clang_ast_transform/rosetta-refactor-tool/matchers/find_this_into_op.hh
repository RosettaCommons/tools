/*
	Find all cases where "this" is being put into an OP or AP

    new claims::TorsionClaim( this, ... );

    `-CXXConstructExpr 0xa3b2490 <col:60> 'ClaimingMoverOP':'class utility::pointer::owning_ptr<class protocols::environment::ClaimingMover>' 'void (pointer)'
      `-ImplicitCastExpr 0xa3b2470 <col:60> 'pointer':'class protocols::environment::ClaimingMover *' <DerivedToBase (ClaimingMover)>
        `-CXXThisExpr 0xa3955a0 <col:60> 'class protocols::abinitio::abscript::AbscriptLoopCloserCM *' this

*/

class ThisIntoOPFinder : public ReplaceMatchCallback {
public:
	ThisIntoOPFinder(
		tooling::Replacements *Replace,
		const char *tag = "ThisIntoOP") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *thisexpr = Result.Nodes.getStmtAs<Expr>("thisexpr");
		const Expr *constructexpr = Result.Nodes.getStmtAs<Expr>("constructexpr");

		if(!rewriteThisFile(thisexpr, sm))
			return;

		const std::string locStr( thisexpr->getSourceRange().getBegin().printToString(sm) );
		const std::string castFrom = QualType::getAsString( thisexpr->getType().split() );
		const std::string castTo = QualType::getAsString( constructexpr->getType().split() );

		llvm::outs()
			<< locStr << "\t"
			<< color("green") << castFrom << color("") << "\t"
			<< color("red") << castTo << color("") << "\n"
			;
	}
};

Finder.addMatcher(
	constructExpr(
		has(
			thisExpr().bind("thisexpr")
		)
	).bind("constructexpr"),
	new ThisIntoOPFinder(Replacements));
