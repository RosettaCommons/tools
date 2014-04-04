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

		const FullSourceLoc FullLocation = FullSourceLoc(construct->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
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
RewriteImplicitCastInMemberCall RewriteImplicitCastInMemberCallCallback(Replacements);
Finder.addMatcher(
	memberCallExpr(
		allOf(
			hasDescendant(
				declRefExpr( containsUtilityPointer() ).bind("declrefexpr")
			),
			hasDescendant(
				materializeTemporaryExpr( isUtilityPointer() ).bind("cast")
			)
		)
	).bind("construct"),
	&RewriteImplicitCastInMemberCallCallback);
