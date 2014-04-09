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
		if(type == "value_type") {
			// Use desugared type since we don't have a better info
			type = QualType::getAsString( cast->getType().getSplitDesugaredType() );
			// owning_ptr -> shared_ptr didn't get rewritten here yet (template?),
			// so do it if it's in the desugared type
			replace(type, "owning_ptr", "shared_ptr");
			replace(type, "access_ptr", "weak_ptr");
		}
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
RewriteImplicitCastFromNew RewriteImplicitCastFromNewCallback1(Replacements,
	"RewriteImplicitCastFromNew:implicitCastExpr>constructExpr>newExpr+implicitCastExpr>newExpr");
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
	&RewriteImplicitCastFromNewCallback1);

// For templates -- no implicit casts found there
/*
MaterializeTemporaryExpr 'value_type':'class utility::pointer::owning_ptr<struct X>' xvalue
`-CXXBindTemporaryExpr 'value_type':'class utility::pointer::owning_ptr<struct X>' (CXXTemporary)
  `-CXXConstructExpr 'value_type':'class utility::pointer::owning_ptr<struct X>' 'void (pointer)'
    `-CXXNewExpr 'struct X *'
*/

RewriteImplicitCastFromNew RewriteImplicitCastFromNewCallback2(Replacements,
	"RewriteImplicitCastFromNew:constructExpr>newExpr+bindTemporaryExpr");
Finder.addMatcher(
	constructExpr(
		allOf(
			isUtilityPointer(),
			has(
				newExpr()
			),
			hasParent(
				bindTemporaryExpr()
			)
		)
	).bind("cast"),
	&RewriteImplicitCastFromNewCallback2);
