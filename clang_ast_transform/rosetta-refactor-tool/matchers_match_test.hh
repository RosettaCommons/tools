////////////////////////////////////////////////////////////////////////////////////////////////////
// Testing matcher -- just reports matching location and code snippet

class MatchTester : public ReplaceMatchCallback {
public:
	MatchTester(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		const FullSourceLoc FullLocation = FullSourceLoc(expr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;

		const std::string origCode = getText(sm, expr);
		const std::string locStr( expr->getSourceRange().getBegin().printToString(sm) );
			
		llvm::errs() 
			<< locStr << "\n" 
			<< "\t" << origCode << "\n"
			<< "\n";
	}
};

MatchTester MatchTesterCallback(Replacements);

Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDescendant(
				implicitCastExpr(isFunctionToPointerDecayCast()).bind("castexpr")
			),
			has(
				memberExpr()
			),
			isUtilityPointer()
		)
	).bind("expr"),
	&MatchTesterCallback);
