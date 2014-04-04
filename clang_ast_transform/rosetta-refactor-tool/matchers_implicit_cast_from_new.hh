////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in "new" instantiation:
//   SomeOP foo() { return new Some; } -- OK

class RewriteImplicitCastFromNew : public ReplaceMatchCallback {
public:
	RewriteImplicitCastFromNew(
			tooling::Replacements *Replace,
			const char *tag = "RewriteImplicitCastFromNew") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *cast = Result.Nodes.getStmtAs<Expr>("cast");

		const FullSourceLoc FullLocation = FullSourceLoc(cast->getLocStart(), sm);
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
	
		doRewrite(sm, cast, origCode, newCode);
	}
};

// SomeOP foo() { return new Some; }
/*
ImplicitCastExpr 'SomeOP':'class utility::pointer::owning_ptr<class Some>' <ConstructorConversion>
`-CXXConstructExpr 'SomeOP':'class utility::pointer::owning_ptr<class Some>' 'void (pointer)'
  `-CXXNewExpr 'class Some *'
or:
  `-ImplicitCastExpr 'x *':'x *' <DerivedToBase>
    `-CXXNewExpr 'X *'
*/
RewriteImplicitCastFromNew RewriteImplicitCastFromNewCallback(Replacements);
Finder.addMatcher(
	implicitCastExpr(
		allOf(
			isUtilityPointer(),
			isConstructorConversionCast(),
			has(
				constructExpr(
					anyOf(
						has(
							newExpr()
						),
						has(
							implicitCastExpr( 
								has(
									newExpr()
								)
							)
						)
					)
				)
			)
		)
	).bind("cast"),
	&RewriteImplicitCastFromNewCallback);
