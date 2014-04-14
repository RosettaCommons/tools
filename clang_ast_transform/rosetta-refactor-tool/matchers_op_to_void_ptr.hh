////////////////////////////////////////////////////////////////////////////////////////////////////
// Find casts from OP or AP to raw pointer, i.e.
//	ClassAAP a(new ClassA);
//	void *p = ap();
//	ClassAOP op(new ClassA);
//	void *p = op();

class MatchTester : public ReplaceMatchCallback {
public:
	MatchTester(
			tooling::Replacements *Replace,
			const char *tag = "OP-to-RawPointer") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
		const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");
		const Stmt *stmt = Result.Nodes.getStmtAs<Stmt>("stmt");

		if(!rewriteThisFile(stmt, sm))
			return;

		const std::string castFromType(
			castFrom ? QualType::getAsString( castFrom->getType().split() ) : ""
		);
		const std::string castToType( 
			castTo ? QualType::getAsString( castTo  ->getType().split() ) : ""
		);

		const std::string castFromTypeD(
			castFrom ? QualType::getAsString( castFrom->getType().getSplitDesugaredType() ) : ""
		);
		const std::string castToTypeD( 
			castTo ? QualType::getAsString( castTo  ->getType().getSplitDesugaredType() ) : ""
		);

		// We only care about cast from utility::pointer to raw (void) pointer
		if(
			!beginsWith(castToTypeD, "pointer (void)") ||
			!checkIsUtilityPointer(castFromTypeD)
		) {
			return;
		}
		
		const std::string origCodeStr = getText(sm, stmt);
		const std::string locStr( stmt->getSourceRange().getBegin().printToString(sm) );
			
		llvm::errs() << locStr << "\n";
		llvm::errs() << "\033[31m" << origCodeStr << "\033[0m\n";

		llvm::errs() << "\033[32mcastFrom: " << castFromType << "\033[0m";
		if(castFromType != castFromTypeD)
			llvm::errs() << " : \033[32m" << castFromTypeD << "\033[0m";
		llvm::errs() << "\n";

		llvm::errs() << "\033[33mcastTo:   " << castToType << "\033[0m";
		if(castToType != castToTypeD)
			llvm::errs() << ":\033[33m" << castToTypeD << "\033[0m";
		llvm::errs() << "\n";

		llvm::outs() << locStr << "\t" << origCodeStr << "\n";
	}
};

MatchTester MatchTesterCallback(Replacements);
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				declRefExpr( isClassOperator() ).bind("castTo")
			),
			has(
				declRefExpr( isNotClassOperator() ).bind("castFrom")
			),
			hasAncestor(
				declStmt().bind("stmt")
			)
		)
	),
	&MatchTesterCallback);
