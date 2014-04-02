////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in constructors
//   given X(SomeAP), using X(new Some) ==> X(SomeAP(new Some))

class RewriteImplicitCastInConstructor : public ReplaceMatchCallback {
public:
	RewriteImplicitCastInConstructor(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *cast = Result.Nodes.getStmtAs<Expr>("cast");
		const Stmt *construct = Result.Nodes.getStmtAs<Stmt>("construct");

		const FullSourceLoc FullLocation = FullSourceLoc(construct->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;
		
		const std::string origCode = getText(sm, cast);

		std::string newCode;
		std::string type( QualType::getAsString( cast->getType().split() ) );
		if(beginsWith(type, "class "))
			type = std::string(type, 6);

		if(origCode == "0" || origCode == "NULL")
			newCode = type + "()";
		else
			newCode = type + "( " + origCode + " )";
	
		doRewrite("RewriteImplicitCastInConstructor", sm, cast, origCode, newCode);
	}
};


// Implicit casts in constructs
RewriteImplicitCastInConstructor RewriteImplicitCastInConstructorCallback(Replacements);
Finder.addMatcher(
	constructExpr(
		allOf(
			hasDescendant(
				implicitCastExpr( isUtilityPointer() )
			),
			hasDescendant(
				implicitCastExpr( isConstructorConversionCast() ).bind("cast")
//				),
//				hasDescendant(
//					bindTemporaryExpr()
			)
		)
	).bind("construct"),
	&RewriteImplicitCastInConstructorCallback);


