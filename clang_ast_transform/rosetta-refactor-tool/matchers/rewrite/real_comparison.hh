/*
	Replace comparison of Reals (doubles or floats):
	double a, b;
	if(a == b) {}
	if(a != b) {}
*/

class RewriteRealComparisons : public ReplaceMatchCallback {
public:
	RewriteRealComparisons(
		tooling::Replacements *Replace,
		const char *tag,
		const char *fn) :
		ReplaceMatchCallback(Replace, tag), fn_(fn) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *a = Result.Nodes.getStmtAs<Expr>("a");
		const Expr *b = Result.Nodes.getStmtAs<Expr>("b");
		const BinaryOperator *op = Result.Nodes.getStmtAs<BinaryOperator>("op");

		if(!rewriteThisFile(op, sm))
			return;

		const std::string origCode = getText(sm, op);
		const std::string aCode = getText(sm, a);
		const std::string bCode = getText(sm, b);

		if(origCode.empty())
			return;

		std::string newCode = fn_ + "( " + aCode + ", " + bCode + " )";
		doRewrite(sm, op, origCode, newCode);
	}
private:
	std::string fn_;
};

// a == b
Finder.addMatcher(
	binaryOperator(
		hasOperatorName("=="),
		hasLHS(ignoringParenImpCasts(
			anyOf(
				declRefExpr(to(varDecl(hasType(isReal())))).bind("a"),
				//declRefExpr(to(varDecl(hasType(isInteger())))).bind("a"),
				floatLiteral().bind("a"),
				integerLiteral().bind("a")
			)
		)),
		hasRHS(ignoringParenImpCasts(
			anyOf(
				declRefExpr(to(varDecl(hasType(isReal())))).bind("b"),
				declRefExpr(to(varDecl(hasType(isInteger())))).bind("b"),
				floatLiteral().bind("b"),
				integerLiteral().bind("b")
			)
		))
	).bind("op"),
	new RewriteRealComparisons(Replacements, "RewriteRealComparisons:==", "numeric::equal_by_epsilon"));

// a != b
Finder.addMatcher(
	binaryOperator(
		hasOperatorName("!="),
		hasLHS(ignoringParenImpCasts(
			anyOf(
				declRefExpr(to(varDecl(hasType(isReal())))).bind("a"),
				floatLiteral().bind("a"),
				integerLiteral().bind("a")
			)
		)),
		hasRHS(ignoringParenImpCasts(
			anyOf(
				declRefExpr(to(varDecl(hasType(isReal())))).bind("b"),
				floatLiteral().bind("b"),
				integerLiteral().bind("b")
			)
		))
	).bind("op"),
	new RewriteRealComparisons(Replacements, "RewriteRealComparisons:!=", "!numeric::equal_by_epsilon"));
