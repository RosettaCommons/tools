////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts to OP/AP
// (too generic -- don't use)

class RewriteImplicitCastToOP : public ReplaceMatchCallback {
public:
	RewriteImplicitCastToOP(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const CXXBindTemporaryExpr *temporaryexpr = Result.Nodes.getStmtAs<CXXBindTemporaryExpr>("temporaryexpr");

		const FullSourceLoc FullLocation = FullSourceLoc(temporaryexpr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;
		
		std::string type = QualType::getAsString( temporaryexpr->getType().split() );
		std::string origCode( getText(sm, temporaryexpr) );
		
		// If original code was 0 then just leave it empty
		std::string newCode(origCode);
		if(origCode == "0" || origCode == "NULL")
			newCode = type + "()";
		else
			newCode = type + "( " + origCode + " )";

		doRewrite("RewriteImplicitCastToOP", sm, temporaryexpr, origCode, newCode);
	}
};

RewriteImplicitCastToOP RewriteImplicitCastToOPCallback(Replacements);
Finder.addMatcher(
		implicitCastExpr(
			allOf(
				hasDescendant(
					bindTemporaryExpr(
						hasDescendant(
							implicitCastExpr( isNonNoopCast() )
						)
					).bind("temporaryexpr")
				),
				isUtilityPointer()
			)
),
&RewriteImplicitCastToOPCallback);
