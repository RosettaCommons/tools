/*
	Replace implicit casts in "new" instantiation in expressions:

	FooOP foo( new Foo );
	FooOP foo = new Foo;
*/

class RewriteCastFromNewExpr : public ReplaceMatchCallback {
public:
	RewriteCastFromNewExpr(
		tooling::Replacements *Replace,
		const char *tag = "CastFromNewExpr") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
		const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");

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

		const std::string origCode = getText(sm, expr);
		const std::string newConstructCode = getText(sm, castFrom);
		std::string prefix;

		size_t n = origCode.find('=');
		if(n == std::string::npos)
			n = origCode.find('(');
		if(n != std::string::npos)
			prefix = trim(origCode.substr(0, n));

		if(Debug) {
			const std::string locStr( expr->getSourceRange().getBegin().printToString(sm) );

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

		doRewrite(sm, expr, origCode, newCode);
		checkAndMarkSourceLocation(castFrom, sm);
	}
};

/*
ExprWithCleanups 0x9e7cb58 <line:45:3, col:59> '<dependent type>'
`-BinaryOperator 0x9e7cb30 <col:3, col:59> '<dependent type>' '='
  |-MemberExpr 0x9e7a5e8 <col:3> 'core::pose::PoseOP':'class utility::pointer::owning_ptr<class core::pose::Pose>' lvalue ->start_pose 0x9e76990
  | `-CXXThisExpr 0x9e7a5d0 <col:3> 'DockingBenchmark<dock, TScale> *' this
  `-CXXFunctionalCastExpr 0x9e7cb08 <col:16, col:59> 'core::pose::PoseOP':'class utility::pointer::owning_ptr<class core::pose::Pose>' functional cast to core::pose::PoseOP <ConstructorConversion>
    `-CXXBindTemporaryExpr 0x9e7cae8 <col:16, col:57> 'core::pose::PoseOP':'class utility::pointer::owning_ptr<class core::pose::Pose>' (CXXTemporary 0x9e7cae0)
      `-CXXConstructExpr 0x9e7caa0 <col:16, col:57> 'core::pose::PoseOP':'class utility::pointer::owning_ptr<class core::pose::Pose>' 'void (pointer)'
        `-CXXNewExpr 0x9e7c9b8 <col:36, col:57> 'core::pose::Pose *'
          `-CXXConstructExpr 0x9e7c988 <col:40, col:57> 'core::pose::Pose':'class core::pose::Pose' 'void (void)'
*/


Finder.addMatcher(
	functionalCastExpr( has(
		bindTemporaryExpr( has(
			constructExpr(
				argumentCountIs(1),
				isUtilityPointer(),
				has(
					newExpr().bind("castFrom")
				)
			).bind("castTo")

		) )
	) ).bind("expr"),
	new RewriteCastFromNewExpr(Replacements, "CastFromNewExpr:functionalCastExpr"));
