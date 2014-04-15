////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace dynamic casts with smart pointers, i.e.
// From:
//	WrappedRealOP val = dynamic_cast< WrappedReal * > ( data()() );
// To:
//	WrappedRealOP val = std::dynamic_pointer_cast< WrappedReal > ( data() );

class RewriteDynCast : public ReplaceMatchCallback {
public:
	RewriteDynCast(
			tooling::Replacements *Replace,
			const char *tag = "RewriteDynCast") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;

		const Expr *dyncastexpr = Result.Nodes.getStmtAs<Expr>("dyncastexpr");
		const Expr *operatorexpr = Result.Nodes.getStmtAs<Expr>("operatorexpr");
		
		if(!rewriteThisFile(dyncastexpr, sm))
			return;

		// Get original code and cast type
		const std::string origCode = getText(sm, dyncastexpr);
		const std::string origOperatorCode = getText(sm, operatorexpr);
		const std::string origCastType = QualType::getAsString( dyncastexpr->getType().split() );
		
		if(origCode.empty() || origOperatorCode.empty() || origCastType.empty())
			return;

		std::string castType( origCastType );
		if(endsWith(castType, "*"))
			castType = trim(std::string(castType, 0, castType.length() -1));

		// origOperatorCode should end with (), so strip that and hope for the best :)
		std::string newOperatorCode = std::string(origOperatorCode, 0, origOperatorCode.length() -2);
		std::string newCode = 
			"utility::pointer::dynamic_pointer_cast< " + castType + " > "
				+ "( " + newOperatorCode +" )";
			
		doRewrite(sm, dyncastexpr, origCode, newCode);
	}
};


/*
CXXDynamicCastExpr 0x7f5caef5af28 <col:22, col:63> 'WrappedReal *' dynamic_cast<WrappedReal *> <Dynamic>
`-CXXOperatorCallExpr 0x7f5caef59de0 <col:54, col:61> 'pointer':'class utility::pointer::ReferenceCount *'
  |-ImplicitCastExpr 0x7f5caef59dc8 <col:60, col:61> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x7f5caef59d48 <col:60, col:61> 'pointer (void) const' lvalue CXXMethod 0x7f5caf1f97e0 'operator()' 'pointer (void) const'
  `-ImplicitCastExpr 0x7f5caef59e20 <col:54, col:59> 'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' <NoOp>
    `-CXXBindTemporaryExpr 0x7f5caef59d28 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' (CXXTemporary 0x7f5caef59d20)
      `-CXXMemberCallExpr 0x7f5caef59cf0 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>'
        `-MemberExpr 0x7f5caef59cc0 <col:54> '<bound member function type>' ->data 0x7f5caef8a990

CXXDynamicCastExpr 0x54731c0 <col:33, col:82> 'const class numeric::expression_parser::FunctionToken *' dynamic_cast<const class numeric::expression_parser::FunctionToken *> <Dynamic>
`-CXXOperatorCallExpr 0x5473158 <col:72, col:81> 'pointer':'const class numeric::expression_parser::Token *'
  |-ImplicitCastExpr 0x5473140 <col:80, col:81> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x5473118 <col:80, col:81> 'pointer (void) const' lvalue CXXMethod 0x5358c50 'operator()' 'pointer (void) const'
  `-ImplicitCastExpr 0x5473198 <col:72> 'const class utility::pointer::owning_ptr<const class numeric::expression_parser::Token>' lvalue <NoOp>
    `-DeclRefExpr 0x54730f0 <col:72> 'TokenCOP':'class utility::pointer::owning_ptr<const class numeric::expression_parser::Token>' lvalue Var 0x5472ea0 'toptoken' 'TokenCOP':'class utility::pointer::owning_ptr<const class numeric::expression_parser::Token>'

*/

RewriteDynCast RewriteDynCastCallback1(Replacements);
Finder.addMatcher(
	dynamicCastExpr(
		has(
			operatorCallExpr(
				allOf(
					// CHILD EXPR: operator() for owning_ptr::operator()
					has(
						declRefExpr( isCallOperator() )
					),
					// CHILD EXPR: castFrom
					anyOf(
						has(
							memberExpr( isUtilityPointer() )
						),
						has(
							declRefExpr( isUtilityPointer() )
						),
						has(
							bindTemporaryExpr( isUtilityPointer() )
						)
					)
				)
			).bind("operatorexpr")
		)
	).bind("dyncastexpr"),
	&RewriteDynCastCallback1);
