////////////////////////////////////////////////////////////////////////////////////////////////////
// Testing matcher -- just reports matching location and code snippet

class MatchTester : public ReplaceMatchCallback {
public:
	MatchTester(tooling::Replacements *Replace) : ReplaceMatchCallback(Replace) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		const FullSourceLoc FullLocation = FullSourceLoc(expr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;

		const std::string origCode = getText(sm, expr);
		const std::string locStr( expr->getSourceRange().getBegin().printToString(sm) );
			
		llvm::errs() 
			<< locStr << "\n" 
			<< "\t" << origCode << "\n"
			<< "\n";
	}
};

MatchTester MatchTesterCallback(Replacements);

/*
	void set_a_vector1(ClassA *a) {
		as_[0] = a;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDescendant(
				implicitCastExpr(
					isLValueToRValueCast()
				)
			),
			isUtilityPointer()
		)
	).bind("expr"),
	&MatchTesterCallback);

/*
	void set_aref_vector1(ClassA & a) {
		as_[0] = &a;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDescendant(
				unaryOperator()
			),
			isUtilityPointer()
		)
	).bind("expr"),
	&MatchTesterCallback);

/*
	void new_a() {
		a_ = new ClassA;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDescendant(
				newExpr()
			),
			isUtilityPointer()
		)
	).bind("expr"),
	&MatchTesterCallback);

/*
	void new_a() {
		a_ = new ClassA;
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDescendant(
				thisExpr()
			),
			isUtilityPointer()
		)
	).bind("expr"),
	&MatchTesterCallback);
	
/*
	void new_a_local() {
		ClassAOP aa = new ClassA;
	}
*/
Finder.addMatcher(
	constructExpr(
		allOf(
			hasDescendant(
				implicitCastExpr( isConstructorConversionCast() )
			),
			isUtilityPointer()
		)
	).bind("expr"),
	&MatchTesterCallback);

/*
	void foo() {
		ClassB b;
		ClassA *a = b(); // should match
	}
*/
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDescendant(
				implicitCastExpr( has( declRefExpr( isVoidPtrOperator() ) ) )
			),
			hasDescendant(
				implicitCastExpr( isUtilityPointer() )
			)
		)
	).bind("expr"),
	&MatchTesterCallback);
