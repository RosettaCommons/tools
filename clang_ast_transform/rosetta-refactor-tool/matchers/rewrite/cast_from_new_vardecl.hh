/*
	Replace implicit casts in "new" instantiation in variable declarations:

	FooOP foo( new Foo );
	FooOP foo = new Foo;
*/

class RewriteCastFromNewVarDecl : public ReplaceMatchCallback {
public:
	RewriteCastFromNewVarDecl(
		tooling::Replacements *Replace,
		const char *tag = "CastFromNewVarDecl") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
		const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");
		const Decl *vardecl = Result.Nodes.getStmtAs<Decl>("vardecl");

		if(!rewriteThisFile(castFrom, sm))
			return;

		// Get castFrom and castTo variable types
		const std::string castFromType(
			castFrom ? QualType::getAsString( castFrom->getType().split() ) : ""
		);
		const std::string castToType(
			castTo ? QualType::getAsString( castTo->getType().split() ) : ""
		);

		const std::string castFromTypeD(
			castFrom ? QualType::getAsString( castFrom->getType().getSplitDesugaredType() ) : ""
		);
		const std::string castToTypeD(
			castTo ? QualType::getAsString( castTo->getType().getSplitDesugaredType() ) : ""
		);

		const std::string origCode = getText(sm, vardecl);
		const std::string newConstructCode = getText(sm, castFrom);
		std::string prefix;

		size_t n = origCode.find('=');
		if(n == std::string::npos)
			n = origCode.find('(');
		if(n != std::string::npos)
			prefix = trim(origCode.substr(0, n));

		if(Debug) {
			const std::string locStr( vardecl->getSourceRange().getBegin().printToString(sm) );

			llvm::errs() << locStr << "\n";
			llvm::errs() << color("red") << origCode << color("") << "\n";

			llvm::errs() << "castFrom: " << color("purple") << castFromType << color("");
			if(castFromType != castFromTypeD)
				llvm::errs() << "          " << color("purple") << castFromTypeD << color("");
			llvm::errs() << "\n";

			llvm::errs() << "castTo:   " << color("brown") << castToType << color("");
			if(castToType != castToTypeD)
				llvm::errs() << "          " << color("brown") << castToTypeD << color("");
			llvm::errs() << "\n\n";
		}

		// Get original code
		if(origCode.empty())
			return;

		// Same thing, do nothing
		if(castFromTypeD == castToTypeD)
			return;

		std::string newCode;
		std::string type( castToType );
		if(type == "value_type" || type == "mapped_type") {
			// Use desugared type since we don't have a better info
			type = castToTypeD;
			// owning_ptr -> shared_ptr didn't get rewritten here yet (template?),
			// so do it if it's in the desugared type
			replace(type, "owning_ptr", "shared_ptr");
			replace(type, "access_ptr", "weak_ptr");
		}
		type = stripQualifiers(type);

		if(origCode == "0" || origCode == "NULL")
			newCode = prefix;
		else
			newCode = prefix + "( " + newConstructCode + " )";

		doRewrite(sm, vardecl, origCode, newCode);
		checkAndMarkSourceLocation(castFrom, sm);
	}
};

/*
	FuncOP func( new func::LinearPenaltyFunction( ... ) );

DeclStmt 0x9ef95e8 <line:84:5, col:123>
`-VarDecl 0x9ef9190 <col:5, col:122> col:33 func 'core::scoring::func::FuncOP':'class utility::pointer::owning_ptr<class core::scoring::func::Func>'
  `-CXXConstructExpr 0x9ef95b0 <col:33, col:122> 'core::scoring::func::FuncOP':'class utility::pointer::owning_ptr<class core::scoring::func::Func>' 'void (pointer)'
    `-ImplicitCastExpr 0x9ef9590 <col:39, col:120> 'pointer':'class core::scoring::func::Func *' <DerivedToBase (Func)>
      `-CXXNewExpr 0x9ef9448 <col:39, col:120> 'core::scoring::func::LinearPenaltyFunction *'
        `-CXXConstructExpr 0x9ef93f8 <col:43, col:120> 'core::scoring::func::LinearPenaltyFunction':'class core::scoring::func::LinearPenaltyFunction' 'void (const Real, const Real, const Real, const Real)'
...
*/

Finder.addMatcher(
	varDecl(
		has(
			constructExpr(
				isUtilityOwningPointer(),
				has(
					newExpr().bind("castFrom")
				)
			).bind("castTo")
		)
	).bind("vardecl"),
	new RewriteCastFromNewVarDecl(Replacements, "CastFromNewVarDecl:constructExpr"));

/*
  DeclStmt 0x5a27b28 <line:124:2, col:103>
  `-VarDecl 0x5a27610 <col:2, col:102> col:33 monteCarlo_ 'protocols::moves::MonteCarloOP':'class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>'
    `-ExprWithCleanups 0x5a27b10 <col:33, col:102> 'protocols::moves::MonteCarloOP':'class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>'
      `-CXXConstructExpr 0x5a27ad8 <col:33, col:102> 'protocols::moves::MonteCarloOP':'class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>' 'void (const class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo> &)' elidable
        `-MaterializeTemporaryExpr 0x5a27ab8 <col:47, col:102> 'const class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>' lvalue
          `-ImplicitCastExpr 0x5a27aa0 <col:47, col:102> 'const class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>' <NoOp>
            `-CXXBindTemporaryExpr 0x5a27a48 <col:47, col:102> 'protocols::moves::MonteCarloOP':'class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>' (CXXTemporary 0x5a27a40)
              `-ImplicitCastExpr 0x5a27a28 <col:47, col:102> 'protocols::moves::MonteCarloOP':'class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>' <ConstructorConversion>
                `-CXXConstructExpr 0x5a279f0 <col:47, col:102> 'protocols::moves::MonteCarloOP':'class utility::pointer::owning_ptr<class protocols::moves::MonteCarlo>' 'void (pointer)'
                  `-CXXNewExpr 0x5a27968 <col:47, col:102> 'protocols::moves::MonteCarlo *'
                    `-CXXConstructExpr 0x5a27920 <col:51, col:102> 'protocols::moves::MonteCarlo':'class protocols::moves::MonteCarlo' 'void (const Pose &, const ScoreFunction &, const Real)'
                      |-ImplicitCastExpr 0x5a27908 <col:80> 'const Pose':'const class core::pose::Pose' lvalue <NoOp>
                      | `-DeclRefExpr 0x5a276d0 <col:80> 'core::pose::Pose':'class core::pose::Pose' lvalue ParmVar 0x5a22700 'pose' 'core::pose::Pose &'
                      |-CXXOperatorCallExpr 0x5a27820 <col:86, col:87> 'const class core::scoring::ScoreFunction':'const class core::scoring::ScoreFunction' lvalue
                      | |-ImplicitCastExpr 0x5a27808 <col:86> 'reference (*)(void) const' <FunctionToPointerDecay>
                      | | `-DeclRefExpr 0x5a27780 <col:86> 'reference (void) const' lvalue CXXMethod 0x3a14c30 'operator*' 'reference (void) const'
                      | `-ImplicitCastExpr 0x5a27768 <col:87> 'const class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction>' lvalue <NoOp>
                      |   `-MemberExpr 0x5a27710 <col:87> 'core::scoring::ScoreFunctionCOP':'class utility::pointer::owning_ptr<const class core::scoring::ScoreFunction>' lvalue ->scorefxn_ 0x3a32fa0
                      |     `-CXXThisExpr 0x5a276f8 <col:87> 'class protocols::simple_moves::asym_fold_and_dock::AsymFoldandDockRbTrialMover *' this
                      `-FloatingLiteral 0x5a27860 <col:98> 'double' 2.000000e+00
*/

Finder.addMatcher(
	varDecl( has(
		exprWithCleanups( has(
			constructExpr(
				argumentCountIs(1),
				has(
					bindTemporaryExpr( has(
						constructExpr(
							isUtilityOwningPointer(),
							has(
								newExpr().bind("castFrom")
							)
						).bind("castTo")
					) )
				)

			)
		) ) 
	) ).bind("vardecl"),
	new RewriteCastFromNewVarDecl(Replacements, "CastFromNewVarDecl:constructExpr-constructExpr"));
