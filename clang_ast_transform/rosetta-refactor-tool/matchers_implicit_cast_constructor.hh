////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in constructors

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
		
		if(!rewriteThisFile(cast, sm))
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
			hasDirect(
				implicitCastExpr( isUtilityPointer() )
			),
			hasDirect(
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
			hasDirect(
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
					hasDirect(
						implicitCastExpr(
							has(	
								bindTemporaryExpr(
									has(
										callExpr(
											hasDirect(
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

#if 0
/*
	ClassAOP aa = NULL;
	
  `-CXXConstructExpr 0x2b706a8 </data/rosetta/tools/clang_ast_transform/test-access_ptr.cc:72:12, /data/rosetta/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (const class utility::pointer::owning_ptr<class ClassA> &)' elidable
    `-MaterializeTemporaryExpr 0x2b70688 <col:18> 'const class utility::pointer::owning_ptr<class ClassA>' lvalue
      `-ImplicitCastExpr 0x2b70670 <col:18> 'const class utility::pointer::owning_ptr<class ClassA>' <NoOp>
        `-CXXBindTemporaryExpr 0x2b70618 <col:18> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x2b70610)
          `-ImplicitCastExpr 0x2b705f0 <col:18> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' <ConstructorConversion>
            `-CXXConstructExpr 0x2b705b8 <col:18> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
               `-ImplicitCastExpr 0x2b705a0 <col:18> 'pointer':'class ClassA *' <NullToPointer>

	ClassAOP aa(NULL);

  `-CXXConstructExpr 0x2a083b8 <col:12, col:19> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
    `-ImplicitCastExpr 0x2a083a0 </data/rosetta/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> 'pointer':'class ClassA *' <NullToPointer>
      `-GNUNullExpr 0x2a082f8 <col:18> 'long'
*/

RewriteImplicitCastInConstructor RewriteImplicitCastInConstructorCallback1(
	Replacements, "RewriteImplicitCastInConstructor:nullPtrLiteralExpr");

Finder.addMatcher(
	constructExpr(
		isUtilityPointer()
	).bind("construct"),
	&RewriteImplicitCastInConstructorCallback1);
#endif
