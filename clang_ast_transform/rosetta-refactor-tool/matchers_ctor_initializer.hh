////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace 0 and NULL in constructor initializers  

class RewriteCtorInitializer : public ReplaceMatchCallback {
	
public:
	RewriteCtorInitializer(
			tooling::Replacements *Replace,
			const char *tag ="RewriteCtorInitializer") :
		ReplaceMatchCallback(Replace, tag) {}
	
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *node = Result.Nodes.getStmtAs<Expr>("literal");

		if(!rewriteThisFile(node, sm))
			return;

		const std::string origCode( getText(sm, node) );
		if(origCode != "0" && origCode != "NULL")
			return;

		std::string newCode = "/* " + origCode + " */";
		doRewrite(sm, node, origCode, newCode);
	}
};

/*
CXXCtorInitializer Field 0x4f43710 'this_weak_ptr_' 'AtomTreeCAP':'class utility::pointer::access_ptr<const class core::kinematics::AtomTree>'
|-CXXConstructExpr 0x6779c38 <line:62:2, col:18> 'AtomTreeCAP':'class utility::pointer::access_ptr<const class core::kinematics::AtomTree>' 'void (pointer)'
| `-ImplicitCastExpr 0x6779c20 <col:17> 'pointer':'const class core::kinematics::AtomTree *' <NullToPointer>
|   `-IntegerLiteral 0x6779b68 <col:17> 'int' 0
*/

RewriteCtorInitializer RewriteCtorInitializerCallback1(Replacements, "CtorInitializer:null");
Finder.addMatcher(
	constructorDecl(
		forEachConstructorInitializer(
			ctorInitializer(
				withInitializer(
					constructExpr(
						has(
							integerLiteral().bind("literal")
						)
					).bind("construct")
				)
			)
		)
	),	&RewriteCtorInitializerCallback1);
