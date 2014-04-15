////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in "new" instantiation

class RewriteImplicitCastFromNew : public ReplaceMatchCallback {
public:
	RewriteImplicitCastFromNew(
			tooling::Replacements *Replace,
			const char *tag = "RewriteImplicitCastFromNew") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
		const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");

		if(!rewriteThisFile(expr, sm))
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
		const std::string locStr( expr->getSourceRange().getBegin().printToString(sm) );
			
#ifdef DEBUG
		llvm::errs() << locStr << "\n";
		llvm::errs() << "\033[31m" << origCode << "\033[0m\n";

		llvm::errs() << "\033[32mcastFrom: " << castFromType << "\033[0m";
		if(castFromType != castFromTypeD)
			llvm::errs() << " : \033[32m" << castFromTypeD << "\033[0m";
		llvm::errs() << "\n";

		llvm::errs() << "\033[33mcastTo:   " << castToType << "\033[0m";
		if(castToType != castToTypeD)
			llvm::errs() << " : \033[33m" << castToTypeD << "\033[0m";
		llvm::errs() << "\n\n";
#endif

		// Get original code
		if(origCode.empty())
			return;
		
		// Same thing, do nothing
		if(castFromTypeD == castToTypeD)
			return;
			
		std::string newCode;
		std::string type( castToType );
		if(type == "value_type") {
			// Use desugared type since we don't have a better info
			type = castToTypeD;
			// owning_ptr -> shared_ptr didn't get rewritten here yet (template?),
			// so do it if it's in the desugared type
			replace(type, "owning_ptr", "shared_ptr");
			replace(type, "access_ptr", "weak_ptr");
		}
		if(beginsWith(type, "class "))
			type = std::string(type, 6);

		if(origCode == "0" || origCode == "NULL")
			newCode = type + "()";
		else
			newCode = type + "( " + origCode + " )";
	
		doRewrite(sm, expr, origCode, newCode);
	}
};

/*
	ClassAOP new_aap() {
		return new ClassA;
	}
	 
CXXMethodDecl 0x3bf6d90 </data/rosetta/tools/clang_ast_transform/test-access_ptr.cc:88:2, line:90:2> line:88:11 new_aap 'ClassAOP (void)'
`-CompoundStmt 0x3bfc3b8 <col:21, line:90:2>
  `-ReturnStmt 0x3bfc398 <line:89:3, col:14>
    `-ExprWithCleanups 0x3bfc380 <col:3, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>'
      `-CXXConstructExpr 0x3bfc348 <col:3, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (const class utility::pointer::owning_ptr<class ClassA> &)' elidable
        `-MaterializeTemporaryExpr 0x3bfc328 <col:10, col:14> 'const class utility::pointer::owning_ptr<class ClassA>' lvalue
          `-ImplicitCastExpr 0x3bfc310 <col:10, col:14> 'const class utility::pointer::owning_ptr<class ClassA>' <NoOp>
            `-CXXBindTemporaryExpr 0x3bfc2b8 <col:10, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x3bfc2b0)
              `-ImplicitCastExpr 0x3bfc298 <col:10, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' <ConstructorConversion>
                `-CXXConstructExpr 0x3bfc260 <col:10, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
                  `-CXXNewExpr 0x3bfc1d8 <col:10, col:14> 'class ClassA *'
                    `-CXXConstructExpr 0x3bfc1a8 <col:14> 'class ClassA' 'void (void)'

*/

RewriteImplicitCastFromNew RewriteImplicitCastFromNewCallback1(Replacements,
	"RewriteImplicitCastFromNew");
Finder.addMatcher(
	constructExpr(
		allOf(
			has(
				bindTemporaryExpr(
					has(
						constructExpr(
							has(
								newExpr().bind("castFrom")
							)
						).bind("expr")
					)
				).bind("castTo")
			),
			isUtilityPointer()			
		)
	),
	&RewriteImplicitCastFromNewCallback1);
