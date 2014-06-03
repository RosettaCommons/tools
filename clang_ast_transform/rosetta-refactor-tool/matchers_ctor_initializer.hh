////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace 0 in constructor initializers  

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

		std::string origCode( getText(sm, node) );
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

RewriteCtorInitializer RewriteCtorInitializerCallback1(Replacements, "CtorInitializer:0");
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




////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace NULL in constructor initializers  

/*
 * NOTE: This is a separate callback because NULL is actually __null in stddef.h,
 * and that's in a file that we don't rewrite here.
 */

class RewriteCtorInitializerHacky : public ReplaceMatchCallback {
	
public:
	RewriteCtorInitializerHacky(
			tooling::Replacements *Replace,
			const char *tag ="RewriteCtorInitializerHacky") :
		ReplaceMatchCallback(Replace, tag) {}
	
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *node = Result.Nodes.getStmtAs<Expr>("literal");
		const Expr *construct = Result.Nodes.getStmtAs<Expr>("construct");

		std::string matchedCode( getText(sm, node) );
		if(matchedCode != "__null")
			return;
		
		std::string origCode( getText(sm, construct) );
		std::string newCode = origCode;
		replace(newCode, "NULL", "/* NULL */");
		doRewrite(sm, construct, origCode, newCode);
	}
};

/*
|-CXXCtorInitializer Field 0x4ea45f0 'pep_prev_' 'NodeAP':'class utility::pointer::access_ptr<class core::environment::FoldTreeSketch::Node>'
| |-CXXConstructExpr 0x5b0dcd8 <line:325:3, col:19> 'NodeAP':'class utility::pointer::access_ptr<class core::environment::FoldTreeSketch::Node>' 'void (pointer)'
| | `-ImplicitCastExpr 0x5b0dcc0 </local/luki/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> 'pointer':'class core::environment::FoldTreeSketch::Node *' <NullToPointer>
| |   `-GNUNullExpr 0x5b0dbe8 <col:18> 'long'
*/

RewriteCtorInitializerHacky RewriteCtorInitializerHackyCallback1(Replacements, "CtorInitializer:nullExpr");
Finder.addMatcher(
	constructorDecl(
		forEachConstructorInitializer(
			ctorInitializer(
				withInitializer(
					constructExpr(
						has(
							expr(isNullExpr()).bind("literal")
						)
					).bind("construct")
				)
			)
		)
	),	&RewriteCtorInitializerHackyCallback1);
