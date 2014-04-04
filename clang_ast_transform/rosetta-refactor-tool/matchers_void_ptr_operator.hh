////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace calls to operator() on owning_ptrs
/*
CXXOperatorCallExpr 'pointer':'const class X *'
|-ImplicitCastExpr 'pointer (*)(void) const' <FunctionToPointerDecay>
| `-DeclRefExpr 'pointer (void) const' lvalue
`-ImplicitCastExpr 'const class utility::pointer::owning_ptr<X>' lvalue <NoOp>
	`-DeclRefExpr 'TokenCOP':'class utility::pointer::owning_ptr<X>' 
*/

class RewriteVoidPtrOperator : public ReplaceMatchCallback {
public:
	RewriteVoidPtrOperator(
			tooling::Replacements *Replace,
			const char *tag = "RewriteVoidPtrOperator") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		
		const FullSourceLoc FullLocation = FullSourceLoc(expr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;

		// Get original code
		const std::string origCode = getText(sm, expr);
		if(origCode.empty())
			return;
			
		if(origCode.find('*') != std::string::npos) {
			// Hack: try to avoid *e1() etc.
			return;
		}

		// origCode should end with (), so strip that call operator
		std::string newCode = "&(*" + std::string(origCode, 0, origCode.length() -2) + ")";
		doRewrite(sm, expr, origCode, newCode);
	}
};


RewriteVoidPtrOperator RewriteVoidPtrOperatorCallback(Replacements);
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				implicitCastExpr(
					has(
						declRefExpr( isVoidPtrOperator() )
					)
				)
			),
			has(
				implicitCastExpr( isUtilityPointer() )
			)
		)
	).bind("expr"),
	&RewriteVoidPtrOperatorCallback);
	
