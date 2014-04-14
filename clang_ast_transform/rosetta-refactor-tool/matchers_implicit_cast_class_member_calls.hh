////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in class member calls

class RewriteImplicitCastInMemberCall : public ReplaceMatchCallback {
public:
	RewriteImplicitCastInMemberCall(
			tooling::Replacements *Replace, 
			const char *tag = "RewriteImplicitCastInMemberCall") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *cast = Result.Nodes.getStmtAs<Expr>("cast");
		const Stmt *construct = Result.Nodes.getStmtAs<Stmt>("construct");
		const DeclRefExpr *declrefexpr = Result.Nodes.getStmtAs<DeclRefExpr>("declrefexpr");

		if(!rewriteThisFile(construct, sm))
			return;
			
		// declrefexpr could be vector1<>, map<>, etc.
		
		const std::string origCode = getText(sm, cast);

		std::string newCode;
		std::string type( QualType::getAsString( declrefexpr->getType().split() ) );

		// Clean up type
		if(beginsWith(type, "class "))
			type = std::string(type, 6);
		
		if(beginsWith(type, "utility::vector1<")) {
			// NOTE: will fail on more complex definitions
			type = trim(std::string(type, 16), "<> ");
		}

		if(origCode == "0" || origCode == "NULL")
			newCode = type + "()";
		else
			newCode = type + "( " + origCode + " )";
	
		doRewrite(sm, cast, origCode, newCode);
	}
};


// Class member call with implicit cast
/*
CXXMemberCallExpr 'void'
|-MemberExpr '<bound member function type>' .push_back 0x7f8ef5807e90
| `-ImplicitCastExpr 'class std::vector<xOP>'...
|   `-DeclRefExpr 'utility::vector1<xOP>':...
`-MaterializeTemporaryExpr 'value_type':'class utility::pointer::owning_ptr<X>' xvalue
*/
RewriteImplicitCastInMemberCall RewriteImplicitCastInMemberCallCallback(Replacements,
	"RewriteImplicitCastInMemberCall:implicitCastExpr>declRefExpr>materializeTemporaryExpr");
Finder.addMatcher(
	memberCallExpr(
		allOf(
			hasDirect(
				implicitCastExpr(
					has(
						declRefExpr( containsUtilityPointer() ).bind("declrefexpr")
					)
				)
			),
			has(
				materializeTemporaryExpr( isUtilityPointer() ).bind("cast")
			)
		)
	).bind("construct"),
	&RewriteImplicitCastInMemberCallCallback);


/*
CXXMemberCallExpr 'void'
`-MemberExpr '<bound member function type>' .push_back 0x7f11c3ae6e90
  `-ImplicitCastExpr 'class std::vector<class utility::pointer::owning_ptr<class X>, class std::allocator<class utility::pointer::owning_ptr<class X> > >' lvalue <UncheckedDerivedToBase (vectorL -> vector)>
    `-DeclRefExpr 'utility::vector1<xOP>':'class utility::vector1<class utility::pointer::owning_ptr<class X>, class std::allocator<class utility::pointer::owning_ptr<class X> > >' lvalue Var 'x' 'utility::vector1<xOP>':'class utility::vector1<class utility::pointer::owning_ptr<class xOP>, class std::allocator<class utility::pointer::owning_ptr<class X> > >'
      `-MaterializeTemporaryExpr 'value_type':'class utility::pointer::owning_ptr<class X>' xvalue
*/

RewriteImplicitCastInMemberCall RewriteImplicitCastInMemberCallCallback2(Replacements,
	"RewriteImplicitCastInMemberCall:memberExpr>declRefExpr+materializeTemporaryExpr");
Finder.addMatcher(
	memberCallExpr(
		allOf(
			has(
				memberExpr(
					hasDirect(
						implicitCastExpr(
							has(
								declRefExpr( containsUtilityPointer() ).bind("declrefexpr")
							)
						)
					)
				)
			),
			has(
				materializeTemporaryExpr( isUtilityPointer() ).bind("cast")
			)
		)
	).bind("construct"),
	&RewriteImplicitCastInMemberCallCallback2);

/*
ImplicitCastExpr 0xa56ff60 </data/rosetta/main/source/src/core/conformation/symmetry/SymmetricConformation.cc:744:19> '_Bool' <UserDefinedConversion>
`-CXXMemberCallExpr 0xa56ff38 <col:19> '_Bool'
  `-MemberExpr 0xa56ff08 <col:19> '<bound member function type>' .operator bool 0xa55aab0
    `-ImplicitCastExpr 0xa56fef0 <col:19> 'const class utility::pointer::owning_ptr<const class core::conformation::Residue>' lvalue <NoOp>
      `-DeclRefExpr 0xa56fd58 <col:19> 'ResidueCOP':'class utility::pointer::owning_ptr<const class core::conformation::Residue>' lvalue Var 0xa558bb0 'vrt_res_cop' 'ResidueCOP':'class utility::pointer::owning_ptr<const class core::conformation::Residue>'
*/
// runtime_assert( someop );
// TODO: doesn't match for some reason
RewriteImplicitCastInMemberCall RewriteImplicitCastInMemberCallCallback3(Replacements,
	"RewriteImplicitCastInMemberCall:memberExpr>declRefExpr<implicitCastExpr");
Finder.addMatcher(
	memberCallExpr(
		allOf(
			has(
				memberExpr(
					hasDirect(
						implicitCastExpr(
							has(
								declRefExpr( isUtilityPointer() ).bind("declrefexpr")
							)
						)
					)
				)
			),
			hasParent(
				implicitCastExpr().bind("cast")
			)
		)
	).bind("construct"),
	&RewriteImplicitCastInMemberCallCallback3);
