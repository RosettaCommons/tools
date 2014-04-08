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

/*
CXXOperatorCallExpr 'pointer':'const class numeric::expression_parser::Expression *'
`-ImplicitCastExpr 'pointer (*)(void) const' <FunctionToPointerDecay>
 `-DeclRefExpr 'pointer (void) const' lvalue CXXMethod 0x3a90cc0 'operator()' 'pointer (void) const'
`-ImplicitCastExpr 'const class utility::pointer::owning_ptr<const class numeric::expression_parser::Expression>' lvalue <NoOp>
  `-DeclRefExpr 'ExpressionCOP':'class utility::pointer::owning_ptr<const class numeric::expression_parser::Expression>' lvalue Var 0x4ab0620 'de1' 'ExpressionCOP':'class utility::pointer::owning_ptr<const class numeric::expression_parser::Expression>'
*/

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

