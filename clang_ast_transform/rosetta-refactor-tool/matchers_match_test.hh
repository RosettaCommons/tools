////////////////////////////////////////////////////////////////////////////////////////////////////
// Testing matcher -- just reports matching location and code snippet

class MatchTester : public ReplaceMatchCallback {
public:
	MatchTester(
			tooling::Replacements *Replace,
			const char *tag = "MatchTester") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		const FullSourceLoc FullLocation = FullSourceLoc(expr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;

		const std::string origCodeStr = getText(sm, expr);
		const std::string locStr( expr->getSourceRange().getBegin().printToString(sm) );
			
		llvm::errs() 
			<< "@ " << locStr << " \033[36m(" << tag << ")\033[0m" 
			<< "- \033[31m" << origCodeStr << "\033[0m\n"
			<< "\n";
	}
};

MatchTester MatchTesterCallback(Replacements);

/************************************************************************
NOTE NOTE NOTE

In llvm/tools/clang/include/clang/ASTMatchers/ASTMatchersInternal.h:
class HasMatcher : public MatcherInterface<T> {
  bool matches(const T &Node, ASTMatchFinder *Finder,
               BoundNodesTreeBuilder *Builder) const override {
    return Finder->matchesChildOf(
        Node, ChildMatcher, Builder,
        //ASTMatchFinder::TK_IgnoreImplicitCastsAndParentheses,
        ASTMatchFinder::TK_AsIs,
        ASTMatchFinder::BK_First);
  }
};

Default traversal mode is to TK_IgnoreImplicitCastsAndParentheses.
It's not really documented and this doesn't help us.
Change to TK_AsIs. This may break something else... ?!
************************************************************************/

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
	&MatchTesterCallback);

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
	&MatchTesterCallback);

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
	&MatchTesterCallback);

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
	&MatchTesterCallback);
	
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
	&MatchTesterCallback);

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
	&MatchTesterCallback);
