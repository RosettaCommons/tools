/*
	Testing matcher -- just reports matching location and code snippet
*/

class MatchTester : public ReplaceMatchCallback {
public:
	MatchTester(
		tooling::Replacements *Replace,
		const char *tag = "MatchTester") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
		const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");

		if(!rewriteThisFile(expr, sm))
			return;

		const std::string castFromType(
			castFrom ? QualType::getAsString( castFrom->getType().split() ) : ""
		);
		const std::string castToType( 
			castTo ? QualType::getAsString( castTo  ->getType().split() ) : ""
		);

		const std::string castFromTypeD(
			castFrom ? QualType::getAsString( castFrom->getType().getSplitDesugaredType() ) : ""
		);
		const std::string castToTypeD( 
			castTo ? QualType::getAsString( castTo  ->getType().getSplitDesugaredType() ) : ""
		);

		const std::string origCode = getText(sm, expr);
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
};

MatchTester *MatchTesterCallback = new MatchTester(Replacements);

/*
	void set_a_vector1(ClassA *a) {
		as_[0] = a;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				implicitCastExpr(
					isLValueToRValueCast()
				)
			),
			isUtilityPointer()
		)
	).bind("expr"),
	MatchTesterCallback);

/*
	void set_aref_vector1(ClassA & a) {
		as_[0] = &a;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				unaryOperator()
			),
			isUtilityPointer()
		)
	).bind("expr"),
	MatchTesterCallback);

/*
	void new_a() {
		a_ = new ClassA;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				newExpr()
			),
			isUtilityPointer()
		)
	).bind("expr"),
	MatchTesterCallback);

/*
	void new_a() {
		a_ = new ClassA;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				thisExpr()
			),
			isUtilityPointer()
		)
	).bind("expr"),
	MatchTesterCallback);
	
/*
	void new_a_local() {
		ClassAOP aa = new ClassA;
	}
*/
Finder.addMatcher(
	constructExpr(
		allOf(
			has(
				implicitCastExpr( isConstructorConversionCast() )
			),
			isUtilityPointer()
		)
	).bind("expr"),
	MatchTesterCallback);

/*
	void foo() {
		ClassB b;
		ClassA *a = b(); // should match
	}
*/
// expr(hasType(recordDecl(hasName("utility::pointer::owning_ptr"))))
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				implicitCastExpr( has( declRefExpr( isVoidPtrOperator() ) ) )
			),
			has(
				implicitCastExpr( isUtilityPointer() )
			)
		)
	).bind("expr"),
	MatchTesterCallback);
