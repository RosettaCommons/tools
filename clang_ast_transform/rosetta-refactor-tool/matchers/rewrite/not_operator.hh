/*
	Replace calls to operator! on owning_ptrs

		!some_weak_op ==> some_weak_op.expired()
		runtime_assert( !some_weak_op );

	This was designed to catch runtime_assert( ! some_cap ) but this will not work because
	the runtime_assert() macro has a different FileID to it's not rewritten!
*/

class RewriteNotOperator : public ReplaceMatchCallback {
public:
	RewriteNotOperator(
		tooling::Replacements *Replace,
		const char *tag = "NotOperator") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
	
		if(!rewriteThisFile(expr, sm))
			return;

		// Get original code
		const std::string origCode = getText(sm, expr);
		if(origCode.empty())
			return;
		
		if(Debug)
			llvm::errs() << color("red") << "origCode: " << origCode << color("") << "\n";
		
		// origCode should end with (), so strip that call operator
		std::string newCode = trim(std::string(origCode, 1)) + ".expired()";
		doRewrite(sm, expr, origCode, newCode);
	}
};

/*
UnaryOperator 0x2a7ad10 <line:64:7, col:20> '_Bool' prefix '!'
`-ParenExpr 0x2a7acf0 <col:8, col:20> '_Bool'
	`-CXXOperatorCallExpr 0x2a7acb0 </local/luki/main/source/src/core/chemical/ResidueTypeSet.hh:92:19, col:21> '_Bool'
		|-ImplicitCastExpr 0x2a7ac98 <col:19> '_Bool (*)(void) const' <FunctionToPointerDecay>
		| `-DeclRefExpr 0x2a7ac70 <col:19> '_Bool (void) const' lvalue CXXMethod 0x2923680 'operator!' '_Bool (void) const'
		`-ImplicitCastExpr 0x2a7ac58 <col:21> 'const class utility::pointer::access_ptr<const class core::chemical::AtomTypeSet>' lvalue <NoOp>
			`-MemberExpr 0x2a7abf0 <col:21> 'AtomTypeSetCAP':'class utility::pointer::access_ptr<const class core::chemical::AtomTypeSet>' lvalue ->atom_types_ 0x2923b40
				`-CXXThisExpr 0x2a7abd8 <col:21> 'class core::chemical::ResidueTypeSet *' this
*/

Finder.addMatcher(
	operatorCallExpr(
		allOf(
			// CHILD EXPR: operator!
			has(
				declRefExpr( isNotOperator() )
			),
			// CHILD EXPR: castFrom
			anyOf(
				has(
					memberExpr( isUtilityPointer() )
				),
				has(
					bindTemporaryExpr( isUtilityPointer() )
				),
				has(
					declRefExpr( isUtilityPointer() )
				)
			),
			// PARENT Expr/Stmt: castTo
			hasParent(
				unaryOperator().bind("expr")
			)
		)
	),
	new RewriteNotOperator(Replacements, "NotOperator:UnaryOperator"));
