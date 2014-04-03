////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in constructors
//   given X(SomeAP), using X(new Some) ==> X(SomeAP(new Some))

class RewriteImplicitCastInConstructor : public ReplaceMatchCallback {
public:
	RewriteImplicitCastInConstructor(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

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
	
		doRewrite("RewriteImplicitCastInConstructor", sm, cast, origCode, newCode);
	}
};


// Implicit casts in constructs
RewriteImplicitCastInConstructor RewriteImplicitCastInConstructorCallback(Replacements);
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
	&RewriteImplicitCastInConstructorCallback);

/*
CXXBindTemporaryExpr 'xOP':'class utility::pointer::owning_ptr<class X>'
`-ImplicitCastExpr 'xOP':'class utility::pointer::owning_ptr<class X>' <ConstructorConversion>
  `-CXXConstructExpr 'xOP':'class utility::pointer::owning_ptr<class X>' 'void (pointer)'
*/

Finder.addMatcher(
	constructExpr(
		allOf(
			hasParent(
				implicitCastExpr( isConstructorConversionCast() ).bind("cast")
			),
			isUtilityPointer()
		)
	).bind("construct"),
	&RewriteImplicitCastInConstructorCallback);


/*
CXXConstructExpr 0x7f12889b8c30 <col:10, col:76> 'utility::sql_database::sessionOP':'class utility::pointer::owning_ptr<class utility::sql_database::session>' 'void (const class utility::pointer::owning_ptr<class utility::sql_database::session> &)' elidable
`-MaterializeTemporaryExpr 0x7f12889b8c10 <col:10, col:76> 'const class utility::pointer::owning_ptr<class utility::sql_database::session>' lvalue
  `-ImplicitCastExpr 0x7f12889b8bf8 <col:10, col:76> 'const class utility::pointer::owning_ptr<class utility::sql_database::session>' <NoOp>
    `-CXXBindTemporaryExpr 0x7f12889b8b98 <col:10, col:76> 'owning_ptr<class utility::sql_database::session>':'class utility::pointer::owning_ptr<class utility::sql_database::session>' (CXXTemporary 0x7f12889b8b90)
      `-CallExpr 0x7f12889b8b40 <col:10, col:76> 'owning_ptr<class utility::sql_database::session>':'class utility::pointer::owning_ptr<class utility::sql_database::session>'
        |-ImplicitCastExpr 0x7f12889b8b28 <col:10, col:55> 'owning_ptr<class utility::sql_database::session> (*)(const ResourceDescription &)' <FunctionToPointerDecay>
        | `-DeclRefExpr 0x7f12889b8a50 <col:10, col:55> 'owning_ptr<class utility::sql_database::session> (const ResourceDescription &)' lvalue Function 0x7f12889b8950 'get_resource' 'owning_ptr<class utility::sql_database::session> (const ResourceDescription &)' (FunctionTemplate 0x7f1288a84800 'get_resource')
*/

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
	&RewriteImplicitCastInConstructorCallback);

