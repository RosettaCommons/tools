////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in assignments to class member variables:
//   class X {
//     SomeOP a_;
//   };
//   SomeAP a_ = new Some;
//   SomeAP a_ = 0;
//   SomeAP a_ = NULL;
//
// Replace implicit casts in assignments local varialbles variables:
//   x(SomeOP a) {
//     a = new Some;
//     a = 0;
//     a = NULL;
//   }
// or:
//   x() {
//     SomeOP a;
//     a = new Some;
//     a = 0;
//     a = NULL;
//   }
// OK

class RewriteImplicitCastInAssignment : public ReplaceMatchCallback {
public:
	RewriteImplicitCastInAssignment(
			tooling::Replacements *Replace, 
			const char *tag = "RewriteImplicitCastInAssignment") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *opercallexpr = Result.Nodes.getStmtAs<Expr>("opercallexpr");
		const CastExpr *castexpr = Result.Nodes.getStmtAs<CastExpr>("castexpr");
		const Expr *declrefexpr = Result.Nodes.getStmtAs<Expr>("expr");
		
		if(!rewriteThisFile(opercallexpr, sm))
			return;

		//llvm::errs() << QualType::getAsString(opercallexpr->getType().split()) << "\n";
		//llvm::errs() << QualType::getAsString(opercallexpr->getType().getSplitDesugaredType()) << "\n";
		
		// Get original code
		const std::string origCode = getText(sm, opercallexpr);
		if(origCode.empty())
			return;
			
		// Get cast type
		std::string type = QualType::getAsString(declrefexpr->getType().split());

		// Rewrite assignment	
		std::string leftSideCode = getTextToDelim(sm, opercallexpr, castexpr);
		if(leftSideCode.find('=') != std::string::npos)
			leftSideCode = std::string(leftSideCode, 0, leftSideCode.find('=')+1);
		std::string rightSideCode(origCode, leftSideCode.length());
		std::string codeOperType = std::string(origCode, leftSideCode.length()-1, 1);
		
		leftSideCode = trim(leftSideCode, " \t\n=");
		rightSideCode = trim(rightSideCode, " \t\n=");

		if(beginsWith(type, "std::map<") || beginsWith(type, "std::set<")) {
			// Hack to use second argument of map definition as type
			size_t j = type.find(',');
			type = trim(std::string(type, j+1), "<>, ");
		}

		//llvm::errs() << "LEFT:  " << leftSideCode << "\n";
		//llvm::errs() << "RIGHT: " << rightSideCode << "\n";
		//llvm::errs() << "OTYPE: " << codeOperType << "\n";
		//llvm::errs() << "TYPE:  " << type << "\n";
		
		if(codeOperType != "=") // we only care about = operators, not [], etc.
			return;
			
		// Skip these casts
		if(leftSideCode.length() <= 0)
			// not an assignment, would produce code with syntax error
			return;
		if(beginsWith(type, "const std::") || beginsWith(type, "std::"))
			// don't cast std:: classes
			return;
		if(type.find("::iterator") != std::string::npos)
			// ignore iterators
			return;
		if(type == "<bound member function type>")
			// not idea what these are, skip 'em
			return;

		if(!endsWith(type, "OP") && !endsWith(type, "AP") && !checkIsUtilityPointer(type)) {
			// Dirty rosetta-specific hack: only rewrite AP and OP;
			// without it, OP to Maps get also rewritten, i.e.:
			// resource_locators_[ locator_tag ] = ResourceLocatorsMap( resource_locator );
			return;
		}
		
		// Full type definition not yet rewritten in original code, so do it here
		replace(type, "owning_ptr", "shared_ptr");
		replace(type, "access_ptr", "weak_ptr");

		std::string newCode;
		if(rightSideCode == "0" || rightSideCode == "NULL") {
			//newCode = leftCode + " " + typedef_name + "()";
			newCode = leftSideCode + ".reset()";
		} else {
			newCode = leftSideCode + " = " + type + "( " + rightSideCode + " )";
			// don't rewrite if same type already (explicit cast)
			if(beginsWith(rightSideCode, type)) {
				return;
			}
		}

		doRewrite(sm, opercallexpr, origCode, newCode);
	}
};


// Assignment to class member OP/APs
/*
CXXOperatorCallExpr 'class utility::pointer::owning_ptr<...>' lvalue
|-ImplicitCastExpr 'class utility::pointer::owning_ptr<...>' <FunctionToPointerDecay>
| `-DeclRefExpr 'class utility::pointer::owning_ptr<...>'
|-MemberExpr 'xCOP':'class utility::pointer::owning_ptr<x>' lvalue ->child_
| `-CXXThisExpr 'class X *' this
		`-ImplicitCastExpr 'const owning_ptr<X>':'const class utility::pointer::owning_ptr<X>' lvalue <NoOp>
		 `-DeclRefExpr 'xOP':'class utility::pointer::owning_ptr<X>' lvalue.
*/

RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback1(Replacements,
	"RewriteImplicitCastInAssignment:operatorCallExpr>memberExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDirect(
				implicitCastExpr( isFunctionToPointerDecayCast() ).bind("castexpr")
			),
			has(
				memberExpr().bind("expr")
			),
			isUtilityPointer()
		)
	).bind("opercallexpr"),
	&RewriteImplicitCastInAssignmentCallback1);

// Assignment to local variable OPs/APs
/*
CXXOperatorCallExpr
|-ImplicitCastExpr 'class utility::pointer::owning_ptr<...> &(*)(pointer)' <FunctionToPointerDecay>
| `-DeclRefExpr ...
|-DeclRefExpr 'xCOP':'class utility::pointer::owning_ptr<...>' lvalue ...
`-ImplicitCastExpr...
  `-CXXNewExpr 'x *'
    `-CXXConstructExpr 'x':'class x' 'void (void)'
*/

RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback2(Replacements,
	"RewriteImplicitCastInAssignment:operatorCallExpr>declRefExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDirect(
				implicitCastExpr( isFunctionToPointerDecayCast() ).bind("castexpr")
			),
			has(
				declRefExpr().bind("expr")
			),
			isUtilityPointer()
		)
	).bind("opercallexpr"),
	&RewriteImplicitCastInAssignmentCallback2);
	
// Assignment to variable with operator, i.e. v_[x] = new Y
/*
CXXOperatorCallExpr 'class utility::pointer::owning_ptr<X>' lvalue
|-ImplicitCastExpr 'class utility::pointer::owning_ptr<X>' <FunctionToPointerDecay>
| `-DeclRefExpr 'class utility::pointer::owning_ptr<X>'
|-CXXOperatorCallExpr 'mapped_type':'class utility::pointer::owning_ptr<X>' lvalue
| |-ImplicitCastExpr 'mapped_type &(*)(const key_type &)' <FunctionToPointerDecay>
| | `-DeclRefExpr ...
| |-MemberExpr ...
| | `-CXXThisExpr ...
| `-DeclRefExpr ...
`-CXXNewExpr 'class X *'
	`-CXXConstructExpr...
		`-DeclRefExpr...
*/

RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback3(Replacements,
	"RewriteImplicitCastInAssignment:operatorCallExpr>operatorCallExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			hasDirect(
				implicitCastExpr(isFunctionToPointerDecayCast()).bind("castexpr")
			),
			has(
				operatorCallExpr(
					has(
						memberExpr().bind("expr")
					)
				)
			),
			isUtilityPointer()
		)
	).bind("opercallexpr"),
	&RewriteImplicitCastInAssignmentCallback3);

/*
CXXConstructExpr 'xOP':'class utility::pointer::owning_ptr<class X>' 'void (const class utility::pointer::owning_ptr<class X> &)' elidable
`-MaterializeTemporaryExpr 'const class utility::pointer::owning_ptr<class X>' lvalue
  `-ImplicitCastExpr 'const class utility::pointer::owning_ptr<class X>' <NoOp>
    `-CXXBindTemporaryExpr 'owning_ptr<X>':'class utility::pointer::owning_ptr<class X>' (CXXTemporary 0x7f12889b8b90)
      `-ImplicitCastExpr 'utility::pointer::owning_ptr<X>':'class utility::pointer::owning_ptr<class X>' <ConstructorConversion>
        `-CXXConstructExpr 'utility::pointer::owning_ptr<X>':'class utility::pointer::owning_ptr<class X>' 'void (pointer)'
          `-CallExpr 'owning_ptr<class X>':'class utility::pointer::owning_ptr<class X>'
            |-ImplicitCastExpr 'owning_ptr<class X> (*)(const ResourceDescription &)' <FunctionToPointerDecay>
            | `-DeclRefExpr 'owning_ptr<class X> (const ResourceDescription &)' lvalue Function 0x7f12889b8950 'get_resource' 'owning_ptr<class X> (const ResourceDescription &)' (FunctionTemplate 0x7f1288a84800 'get_resource')
*/

RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback4(
	Replacements, "RewriteImplicitCastInAssignment:constructExpr>temporaryExpr>>>callExpr");

Finder.addMatcher(
	constructExpr(
		allOf(
			has(
				materializeTemporaryExpr(
					has(
						bindTemporaryExpr(
							has(
								constructExpr(
									has(
										callExpr(
											has(
												implicitCastExpr( isFunctionToPointerDecayCast() )
											)
										)
									)
								).bind("expr")
							)
						).bind("castexpr")
					)
				)
			),
			isUtilityPointer()
		)
	).bind("opercallexpr"),
	&RewriteImplicitCastInAssignmentCallback4);
