////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in variable (others?) declataions
//   SomeCAP x(0); ==> SomeCAP x();

class RewriteImplicitCastInDecl : public ReplaceMatchCallback {
public:
	RewriteImplicitCastInDecl(
			tooling::Replacements *Replace,
			const char *replacement ="",
			const char *tag = "RewriteImplicitCastInDecl") :
		ReplaceMatchCallback(Replace, tag), replacement(replacement) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *construct = Result.Nodes.getStmtAs<Expr>("construct");

		if(!rewriteThisFile(construct, sm))
			return;

		const std::string origCode = getText(sm, construct);
		size_t splitPos = origCode.find('(');
		if(splitPos == std::string::npos)
			return;

		const std::string leftSide = std::string(origCode, 0, splitPos);
		const std::string inside = trim(std::string(origCode, splitPos+1), "()\t ");
		
		if(inside != "0" && inside != "NULL")
			return;
		
		std::string newCode = trim(leftSide, "\t\r\n ") + replacement;
		doRewrite(sm, construct, origCode, newCode);
	}
private:
	const char *replacement;
};

/*
// Method variable declaration of an AP/OR with a NULL initializer
DeclStmt
`-VarDecl 'xCAP':'class utility::pointer::access_ptr<const X>'
  `-CXXConstructExpr 'xCAP':'class utility::pointer::access_ptr<const class X>' 'void (pointer)'
    `-ImplicitCastExpr 'pointer':'const class X *' <NullToPointer>
      `-IntegerLiteral 'int' 0
*/

RewriteImplicitCastInDecl RewriteImplicitCastInDeclCallback1(
	Replacements, 
	"",
	"RewriteImplicitCastInDecl:constructExpr>implitiCastExpt.NullToPointer+varDecl");
Finder.addMatcher(
	constructExpr(
		allOf(
			hasDirect(
				implicitCastExpr( isNullToPointerCast() )
			),
			hasParent(
				varDecl()
			),
			isUtilityPointer()
		)
	).bind("construct"),
	&RewriteImplicitCastInDeclCallback1);


/*
// c'tor initializer of an AP/OP with a NULL
CXXCtorInitializer Field 0x7fa4666ab230 'this_weak_ptr_' 'AtomAP':'class utility::pointer::access_ptr<class core::kinematics::tree::Atom>'
`-CXXConstructExpr 'AtomAP':'class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' 'void (pointer)'
  `-ImplicitCastExpr 'pointer':'class core::kinematics::tree::Atom *' <NullToPointer>
    `-IntegerLiteral 'int' 0
*/

RewriteImplicitCastInDecl RewriteImplicitCastInDeclCallback2(
	Replacements,
	"()",
	"RewriteImplicitCastInDecl:constructExpr>implitiCastExpt.NullToPointer+ctorInitializer");
Finder.addMatcher(
	constructExpr(
		allOf(
			hasDirect(
				implicitCastExpr( isNullToPointerCast() )
			),
			isUtilityPointer()
		)
	).bind("construct"),
	&RewriteImplicitCastInDeclCallback2);
