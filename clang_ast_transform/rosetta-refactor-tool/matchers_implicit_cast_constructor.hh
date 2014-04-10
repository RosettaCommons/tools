////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in constructors
//   given X(SomeAP), using X(new Some) ==> X(SomeAP(new Some))

class RewriteImplicitCastInConstructor : public ReplaceMatchCallback {
public:
	RewriteImplicitCastInConstructor(
			tooling::Replacements *Replace,
			const char *tag = "RewriteImplicitCastInConstructor") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *construct = Result.Nodes.getStmtAs<Expr>("construct");
		const Expr *cast = Result.Nodes.getStmtAs<Expr>("cast") ?
			Result.Nodes.getStmtAs<Expr>("cast") :
			construct;
		
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
	
		doRewrite(sm, cast, origCode, newCode);
	}
};


// Implicit casts in constructs
RewriteImplicitCastInConstructor RewriteImplicitCastInConstructorCallback1(
	Replacements, "RewriteImplicitCastInConstructor:constructExpr>implicitCast.ConstructorConversion");
Finder.addMatcher(
	constructExpr(
		allOf(
			has(
				implicitCastExpr( isUtilityPointer() )
			),
			has(
				implicitCastExpr( isConstructorConversionCast() ).bind("cast")
			)
		)
	).bind("construct"),
	&RewriteImplicitCastInConstructorCallback1);

/*
CXXBindTemporaryExpr 'XOP':'class utility::pointer::owning_ptr<class X>' (CXXTemporary 0x6454ef0)
`-ImplicitCastExpr 'XOP':'class utility::pointer::owning_ptr<class X>' <ConstructorConversion>
  `-CXXConstructExpr 'XOP':'class utility::pointer::owning_ptr<class X>' 'void (pointer)'
    `-ImplicitCastExpr 'X *' <LValueToRValue>
*/

RewriteImplicitCastInConstructor RewriteImplicitCastInConstructorCallback2(
	Replacements, "RewriteImplicitCastInConstructor:constructExpr<implicitCast.ConstructorConversion");
Finder.addMatcher(
	constructExpr(
		allOf(
			hasParent(
				implicitCastExpr( isConstructorConversionCast() ).bind("cast")
			),
			has(
				implicitCastExpr( isLValueToRValueCast() )
			),
			isUtilityPointer()
		)
	).bind("construct"),
	&RewriteImplicitCastInConstructorCallback2);
	
/*
CXXConstructExpr 'xOP':'class utility::pointer::owning_ptr<class X>' 'void (const class utility::pointer::owning_ptr<class X> &)' elidable
`-MaterializeTemporaryExpr 'const class utility::pointer::owning_ptr<class X>' lvalue
  `-ImplicitCastExpr 'const class utility::pointer::owning_ptr<class X>' <NoOp>
    `-CXXBindTemporaryExpr 'owning_ptr<X>':'class utility::pointer::owning_ptr<class X>' (CXXTemporary 0x7f12889b8b90)
      `-CallExpr 'owning_ptr<class X>':'class utility::pointer::owning_ptr<class X>'
        |-ImplicitCastExpr 'owning_ptr<class X> (*)(const ResourceDescription &)' <FunctionToPointerDecay>
        | `-DeclRefExpr 'owning_ptr<class X> (const ResourceDescription &)' lvalue Function 0x7f12889b8950 'get_resource' 'owning_ptr<class X> (const ResourceDescription &)' (FunctionTemplate 0x7f1288a84800 'get_resource')
*/

RewriteImplicitCastInConstructor RewriteImplicitCastInConstructorCallback3(
	Replacements, "RewriteImplicitCastInConstructor:constructExpr>temporaryExpr>callExpr");

Finder.addMatcher(
	constructExpr(
		allOf(
			has(
				materializeTemporaryExpr(
					has(
						implicitCastExpr(
							has(	
								bindTemporaryExpr(
									has(
										callExpr(
											has(
												implicitCastExpr( isFunctionToPointerDecayCast() )
											)
										)
									)
								)
							)	
						)
					)
				)
			),
			isUtilityPointer()
		)
	).bind("construct"),
	&RewriteImplicitCastInConstructorCallback3);


