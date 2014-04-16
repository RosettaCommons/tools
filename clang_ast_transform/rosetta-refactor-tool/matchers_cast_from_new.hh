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
		const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
		const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");
		const Expr *expr = castFrom;

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
		if(type == "value_type" || type == "mapped_type") {
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
	 
  CXXConstructExpr 0x3bfc348 <col:3, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (const class utility::pointer::owning_ptr<class ClassA> &)' elidable
  `-MaterializeTemporaryExpr 0x3bfc328 <col:10, col:14> 'const class utility::pointer::owning_ptr<class ClassA>' lvalue
    `-ImplicitCastExpr 0x3bfc310 <col:10, col:14> 'const class utility::pointer::owning_ptr<class ClassA>' <NoOp>
      `-CXXBindTemporaryExpr 0x3bfc2b8 <col:10, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x3bfc2b0)
        `-ImplicitCastExpr 0x3bfc298 <col:10, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' <ConstructorConversion>
          `-CXXConstructExpr 0x3bfc260 <col:10, col:14> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
            `-CXXNewExpr 0x3bfc1d8 <col:10, col:14> 'class ClassA *'
              `-CXXConstructExpr 0x3bfc1a8 <col:14> 'class ClassA' 'void (void)'
*/

RewriteImplicitCastFromNew RewriteImplicitCastFromNewCallback1(Replacements,
	"RewriteImplicitCastFromNew:constructExpr");
Finder.addMatcher(
	newExpr(
		hasParent(
			constructExpr().bind("castTo")
		)
	).bind("castFrom"),
	&RewriteImplicitCastFromNewCallback1);


/*
	variables_[ varname ] = new VariableExpression( varname );

  CXXOperatorCallExpr 0x47d5c80 <line:1411:2, col:58> 'class utility::pointer::owning_ptr<class numeric::expression_parser::VariableExpression>' lvalue
  |-ImplicitCastExpr 0x47d5c68 <col:24> 'class utility::pointer::owning_ptr<class numeric::expression_parser::VariableExpression> &(*)(pointer)' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x47d5be8 <col:24> 'class utility::pointer::owning_ptr<class numeric::expression_parser::VariableExpression> &(pointer)' lvalue CXXMethod 0x47c4480 'operator=' 'class utility::pointer::owning_ptr<class numeric::expression_parser::VariableExpression> &(pointer)'
  |-CXXOperatorCallExpr 0x47d5600 <col:2, col:22> 'mapped_type':'class utility::pointer::owning_ptr<class numeric::expression_parser::VariableExpression>' lvalue
  | |-ImplicitCastExpr 0x47d55e8 <col:12, col:22> 'mapped_type &(*)(const key_type &)' <FunctionToPointerDecay>
  | | `-DeclRefExpr 0x47d5568 <col:12, col:22> 'mapped_type &(const key_type &)' lvalue CXXMethod 0x3859930 'operator[]' 'mapped_type &(const key_type &)'
  | |-MemberExpr 0x47d5510 <col:2> 'std::map<std::string, VariableExpressionOP>':'class std::map<class std::basic_string<char>, class utility::pointer::owning_ptr<class numeric::expression_parser::VariableExpression>, struct std::less<class std::basic_string<char> >, class std::allocator<struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class numeric::expression_parser::VariableExpression> > > >' lvalue ->variables_ 0x385ef10
  | | `-CXXThisExpr 0x47d54f8 <col:2> 'class numeric::expression_parser::SimpleExpressionCreator *' this
  | `-DeclRefExpr 0x47d5540 <col:14> 'const std::string':'const class std::basic_string<char>' lvalue ParmVar 0x47c0d80 'varname' 'const std::string &'
  `-CXXNewExpr 0x47d57a0 <col:26, col:58> 'class numeric::expression_parser::VariableExpression *'
    `-CXXConstructExpr 0x47d5768 <col:30, col:58> 'class numeric::expression_parser::VariableExpression' 'void (const std::string &)'
      `-DeclRefExpr 0x47d5648 <col:50> 'const std::string':'const class std::basic_string<char>' lvalue ParmVar 0x47c0d80 'varname' 'const std::string &'
*/

RewriteImplicitCastFromNew RewriteImplicitCastFromNewCallback2(Replacements,
	"RewriteImplicitCastFromNew:operatorCallExpr");
Finder.addMatcher(
	newExpr(
		hasParent(
			operatorCallExpr(
				has(
					declRefExpr(
						allOf(
							isNotClassOperator(),
							isUtilityPointer()
						)
					).bind("castTo")
				)
			)
		)
	).bind("castFrom"),
	&RewriteImplicitCastFromNewCallback2);
