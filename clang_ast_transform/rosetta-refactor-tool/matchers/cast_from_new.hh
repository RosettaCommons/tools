/*
	Replace implicit casts in "new" instantiation:

		ClassAOP new_aap() { return new ClassA; }
		AtomCOP foo() { return this; }
		foo(AtomCOP a); foo(0);

	variables_[ varname ] = new VariableExpression( varname );
*/

class RewriteCastFromNew : public ReplaceMatchCallback {
public:
	RewriteCastFromNew(
		tooling::Replacements *Replace,
		const char *tag = "CastFromNew") :
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

          ...
          `-CXXConstructExpr 0x47754b0 <col:41, col:64> 'InterpolatorOP':'class utility::pointer::owning_ptr<class numeric::interpolation::spline::Interpolator>' 'void (pointer)'
            `-ImplicitCastExpr 0x4775490 <col:41, col:64> 'pointer':'class numeric::interpolation::spline::Interpolator *' <DerivedToBase (Interpolator)>
              `-CXXNewExpr 0x4775400 <col:41, col:64> 'class numeric::interpolation::spline::SimpleInterpolator *'
                `-CXXConstructExpr 0x47753d0 <col:45, col:64> 'class numeric::interpolation::spline::SimpleInterpolator' 'void (void)'

*/

Finder.addMatcher(
	newExpr(
		anyOf(
			hasParent(
				constructExpr(isUtilityPointer()).bind("castTo")
			),
			hasParent(
				implicitCastExpr(
					hasParent(
						constructExpr(isUtilityPointer()).bind("castTo")
					)
				)
			)
		)
	).bind("castFrom"),
	new RewriteCastFromNew(Replacements, "CastFromNew:constructExpr"));


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
	new RewriteCastFromNew(Replacements, "CastFromNew:operatorCallExpr"));


/*
AtomCOP foo() { return this; }

`-ReturnStmt 0x43e7538 <line:228:3, col:10>
  `-ExprWithCleanups 0x43e7520 <col:3, col:10> 'AtomCOP':'class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>'
    `-CXXConstructExpr 0x43e74e8 <col:3, col:10> 'AtomCOP':'class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' 'void (const class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom> &)' elidable
      `-MaterializeTemporaryExpr 0x43e74c8 <col:10> 'const class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' lvalue
        `-ImplicitCastExpr 0x43e74b0 <col:10> 'const class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' <NoOp>
          `-CXXBindTemporaryExpr 0x43e7458 <col:10> 'AtomCOP':'class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' (CXXTemporary 0x43e7450)
            `-ImplicitCastExpr 0x43e7430 <col:10> 'AtomCOP':'class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' <ConstructorConversion>
              `-CXXConstructExpr 0x43e73f8 <col:10> 'AtomCOP':'class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' 'void (pointer)'
                `-ImplicitCastExpr 0x43e73d0 <col:10> 'pointer':'const class core::kinematics::tree::Atom *' <DerivedToBase (Atom_ -> Atom)>
                  `-CXXThisExpr 0x43e7380 <col:10> 'const class core::kinematics::tree::BondedAtom *' this
*/

// Let's not rewrite those automatically -- don't put this into an OP automatically
/*
Finder.addMatcher(
	thisExpr(
		hasParent(
			implicitCastExpr(
				hasParent(
					constructExpr( isUtilityPointer() ).bind("castTo")
				)
			)
		)
	).bind("castFrom"),
	new RewriteCastFromNew(Replacements, "CastFromNew:thisExpr"));
*/


/*
AtomTree.cc:1152
	root_ = src.root_->clone( 0, atom_pointer_ );
first argument is an AtomAP

`-CXXConstructExpr 0x6b8b3b8 <col:25> 'AtomAP':'class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' 'void (pointer)'
  `-ImplicitCastExpr 0x6b8b3a0 <col:25> 'pointer':'class core::kinematics::tree::Atom *' <NullToPointer>
    `-IntegerLiteral 0x6b8b2c0 <col:25> 'int' 0
*/


Finder.addMatcher(
	constructExpr(
		allOf(
			isUtilityAccessPointer(),
			has(
				integerLiteral().bind("castFrom")
			)
		)
	).bind("castTo"),
	new RewriteCastFromNew(Replacements, "CastFromNew:integerLiteral"));

