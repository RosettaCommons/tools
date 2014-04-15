////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in owning/access pointers

class RewriteAssignmentsOper : public ReplaceMatchCallback {
public:
	RewriteAssignmentsOper(
			tooling::Replacements *Replace, 
			const char *tag = "RewriteAssignmentsOper") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		const Expr *castexpr = Result.Nodes.getStmtAs<Expr>("castexpr");
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

		// Rewrite assignment	
		std::string type(castToType);
		std::string leftSideCode = getTextToDelim(sm, expr, castexpr);
		if(leftSideCode.find('=') != std::string::npos)
			leftSideCode = std::string(leftSideCode, 0, leftSideCode.find('=')+1);
		std::string rightSideCode(origCode, leftSideCode.length());
		std::string codeOperType = std::string(origCode, leftSideCode.length()-1, 1);
		
		leftSideCode = trim(leftSideCode, " \t\n=");
		rightSideCode = trim(rightSideCode, " \t\n=");

#ifdef DEBUG
		llvm::errs() << "LEFT:  " << leftSideCode << "\n";
		llvm::errs() << "RIGHT: " << rightSideCode << "\n";
		llvm::errs() << "OTYPE: " << codeOperType << "\n";
		llvm::errs() << "TYPE:  " << type << "\n";
#endif
		
		// Skip these casts
		if(leftSideCode.length() <= 0)
			// not an assignment, would produce code with syntax error
			return;
		if(type.find("::iterator") != std::string::npos)
			// ignore iterators
			return;

		// Handle OPs contained in vectors, lists, maps, sets
		if(
			beginsWith(castToTypeD, "class utility::vector0<") ||
			beginsWith(castToTypeD, "class utility::vector1<") ||
			beginsWith(castToTypeD, "class std::map<") ||
			beginsWith(castToTypeD, "class std::vector<") ||
			beginsWith(castToTypeD, "class std::list<") ||
			beginsWith(castToTypeD, "class std::set<")
		) {
			const std::string castToTypeD_contained = extractTypeFromContainer(castToTypeD);
			if(castFromTypeD == castToTypeD_contained)
				return;
			type = extractTypeFromContainer(type);
		}

		if(beginsWith(type, "class "))
			type = std::string(type, strlen("class "));
			
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

		doRewrite(sm, expr, origCode, newCode);	}

protected:
	std::string extractTypeFromContainer(const std::string & container) {
		// Strip utility::vector1<utility::pointer::ownining_ptr<class ClassX>, allocator ... >
		// TODO: handle std::map<>
		size_t p = container.find('<');
		size_t q = container.find_last_of('>');
		if(p == std::string::npos)
			return container;
			
		std::string container_type = std::string(container, p + 1, q - p - 1); 
		q = container_type.find(',');
		if(q != std::string::npos)
			container_type = trim(std::string(container_type, 0, q));
		return container_type;
	}

};

/*
  ClassAOP a_;
  ClassA *a;
  a_ = a;

  `-CXXOperatorCallExpr 0x37bbff0 <line:47:3, col:8> 'class utility::pointer::owning_ptr<class ClassA>' lvalue
    |-ImplicitCastExpr 0x37bbfd8 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(*)(pointer)' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x37bbfb0 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(pointer)' lvalue CXXMethod 0x3769bc0 'operator=' 'class utility::pointer::owning_ptr<class ClassA> &(pointer)'
    |-MemberExpr 0x37bbb80 <col:3> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' lvalue ->a_ 0x376a820
    | `-CXXThisExpr 0x37bbb68 <col:3> 'class ClassB *' this
    `-ImplicitCastExpr 0x37bbf98 <col:8> 'class ClassA *' <LValueToRValue>
      `-DeclRefExpr 0x37bbbb0 <col:8> 'class ClassA *' lvalue ParmVar 0x37b7de0 'a' 'class ClassA *'

	utility::vector1<ClassAOP> as_;
	void set_a_vector1(ClassA *a) {
		as_[0] = a;
	}

  `-CXXOperatorCallExpr 0x2978818 <line:103:3, col:12> 'class utility::pointer::owning_ptr<class ClassA>' lvalue
    |-ImplicitCastExpr 0x2978800 <col:10> 'class utility::pointer::owning_ptr<class ClassA> &(*)(pointer)' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x29787d8 <col:10> 'class utility::pointer::owning_ptr<class ClassA> &(pointer)' lvalue CXXMethod 0x29231f0 'operator=' 'class utility::pointer::owning_ptr<class ClassA> &(pointer)'
    |-CXXOperatorCallExpr 0x29786f0 <col:3, col:8> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' lvalue
    | |-ImplicitCastExpr 0x29786d8 <col:6, col:8> 'reference (*)(const index_type)' <FunctionToPointerDecay>
    | | `-DeclRefExpr 0x2978650 <col:6, col:8> 'reference (const index_type)' lvalue CXXMethod 0x2962800 'operator[]' 'reference (const index_type)'
    | |-ImplicitCastExpr 0x2978618 <col:3> 'class utility::vectorL<1, class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue <UncheckedDerivedToBase (vectorL)>
    | | `-MemberExpr 0x29785a0 <col:3> 'utility::vector1<ClassAOP>':'class utility::vector1<class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue ->as_ 0x2970c60
    | |   `-CXXThisExpr 0x2978588 <col:3> 'class ClassB *' this
    | `-ImplicitCastExpr 0x2978638 <col:7> 'index_type':'unsigned long' <IntegralCast>
    |   `-IntegerLiteral 0x29785d0 <col:7> 'int' 0
    `-ImplicitCastExpr 0x29787c0 <col:12> 'class ClassA *' <LValueToRValue>
      `-DeclRefExpr 0x2978738 <col:12> 'class ClassA *' lvalue ParmVar 0x29724a0 'a' 'class ClassA *'
*/

RewriteAssignmentsOper RewriteAssignmentsOperCallback1(
	Replacements,
	"RewriteAssignmentsOper"
);
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				declRefExpr( isClassOperator() ).bind("castexpr")
			),
			anyOf(
				// a_ = a;
				has(
					memberExpr( isUtilityPointer() ).bind("castTo")
				),
				// as_[0] = a;
				has(
					operatorCallExpr(
						has(
							memberExpr( containsUtilityPointer() ).bind("castTo")
						)
					)
				)
			),
			has(
				declRefExpr( isNotClassOperator() ).bind("castFrom")
			)
		)
	).bind("expr"),
	&RewriteAssignmentsOperCallback1);

/*
  ClassAOP a_;
  ClassA & a;
  a_ = &a;

  `-CXXOperatorCallExpr 0x3f01580 <line:60:3, col:9> 'class utility::pointer::owning_ptr<class ClassA>' lvalue
    |-ImplicitCastExpr 0x3f01568 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(*)(pointer)' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x3f01540 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(pointer)' lvalue CXXMethod 0x3eaec20 'operator=' 'class utility::pointer::owning_ptr<class ClassA> &(pointer)'
    |-MemberExpr 0x3f01468 <col:3> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' lvalue ->a_ 0x3eaf880
    | `-CXXThisExpr 0x3f01450 <col:3> 'class ClassB *' this
    `-UnaryOperator 0x3f014c0 <col:8, col:9> 'class ClassA *' prefix '&'
      `-DeclRefExpr 0x3f01498 <col:9> 'class ClassA' lvalue ParmVar 0x3efd1c0 'a' 'class ClassA &'

	utility::vector1<ClassAOP> as_;
	void set_aref_vector1(ClassA & a) {
		as_[0] = &a;
	}

  `-CXXOperatorCallExpr 0x41ec4d0 <line:31:3, col:13> 'class utility::pointer::owning_ptr<class ClassA>' lvalue
    |-ImplicitCastExpr 0x41ec4b8 <col:10> 'class utility::pointer::owning_ptr<class ClassA> &(*)(pointer)' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x41ec490 <col:10> 'class utility::pointer::owning_ptr<class ClassA> &(pointer)' lvalue CXXMethod 0x41eac50 'operator=' 'class utility::pointer::owning_ptr<class ClassA> &(pointer)'
    |-CXXOperatorCallExpr 0x41ec398 <col:3, col:8> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' lvalue
    | |-ImplicitCastExpr 0x41ec380 <col:6, col:8> 'reference (*)(const index_type)' <FunctionToPointerDecay>
    | | `-DeclRefExpr 0x41ec358 <col:6, col:8> 'reference (const index_type)' lvalue CXXMethod 0x41d8f10 'operator[]' 'reference (const index_type)'
    | |-ImplicitCastExpr 0x41ec320 <col:3> 'class utility::vectorL<1, class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue <UncheckedDerivedToBase (vectorL)>
    | | `-MemberExpr 0x41ec2d0 <col:3> 'utility::vector1<ClassAOP>':'class utility::vector1<class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue ->as_ 0x41e22c0
    | |   `-CXXThisExpr 0x41ec2b8 <col:3> 'class ClassA *' this
    | `-ImplicitCastExpr 0x41ec340 <col:7> 'index_type':'unsigned long' <IntegralCast>
    |   `-IntegerLiteral 0x41ec300 <col:7> 'int' 0
    `-UnaryOperator 0x41ec408 <col:12, col:13> 'class ClassA *' prefix '&'
      `-DeclRefExpr 0x41ec3e0 <col:13> 'class ClassA' lvalue ParmVar 0x41e2c40 'a' 'class ClassA &'
*/

RewriteAssignmentsOper RewriteAssignmentsOperCallback2(
	Replacements,
	"RewriteAssignmentsOper:unaryOperator"
);
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				declRefExpr( isClassOperator() ).bind("castexpr")
			),
			anyOf(
				// a_ = a;
				has(
					memberExpr( isUtilityPointer() ).bind("castTo")
				),
				// as_[0] = a;
				has(
					operatorCallExpr(
						has(
							memberExpr( containsUtilityPointer() ).bind("castTo")
						)
					)
				)
			),
			has(
				unaryOperator(
					has(
						declRefExpr( isNotClassOperator() ).bind("castFrom")
					)
				)
			)
		)
	).bind("expr"),
	&RewriteAssignmentsOperCallback2);

/*
  ClassAOP a_;
  a_ = new ClassA;

  `-CXXOperatorCallExpr 0x26f1220 <line:65:3, col:12> 'class utility::pointer::owning_ptr<class ClassA>' lvalue
    |-ImplicitCastExpr 0x26f1208 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(*)(pointer)' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x26f11e0 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(pointer)' lvalue CXXMethod 0x269dba0 'operator=' 'class utility::pointer::owning_ptr<class ClassA> &(pointer)'
    |-MemberExpr 0x26f0c90 <col:3> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' lvalue ->a_ 0x269e800
    | `-CXXThisExpr 0x26f0c78 <col:3> 'class ClassB *' this
    `-CXXNewExpr 0x26f1120 <col:8, col:12> 'class ClassA *'
      `-CXXConstructExpr 0x26f10f0 <col:12> 'class ClassA' 'void (void)'
*/

RewriteAssignmentsOper RewriteAssignmentsOperCallback3(
	Replacements,
	"RewriteAssignmentsOper:newExpr"
);
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				declRefExpr( isClassOperator() ).bind("castexpr")
			),
			has(
				memberExpr( isUtilityPointer() ).bind("castTo")
			),
			has(
				newExpr().bind("castFrom")
			)
		)
	).bind("expr"),
	&RewriteAssignmentsOperCallback3);


/*
class X {
  SomeOP a_;
	void set_a_null() {
		a_ = NULL;
		a_ = 0;
	}
};

  |-CXXOperatorCallExpr 0x403dfa0 <line:42:3, /data/rosetta/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> 'class utility::pointer::owning_ptr<class ClassA>' lvalue
  | |-ImplicitCastExpr 0x403df88 </data/rosetta/tools/clang_ast_transform/test-access_ptr.cc:42:6> 'class utility::pointer::owning_ptr<class ClassA> &(*)(pointer)' <FunctionToPointerDecay>
  | | `-DeclRefExpr 0x403df08 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(pointer)' lvalue CXXMethod 0x3febb50 'operator=' 'class utility::pointer::owning_ptr<class ClassA> &(pointer)'
  | |-MemberExpr 0x403de40 <col:3> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' lvalue ->a_ 0x3fec7b0
  | | `-CXXThisExpr 0x403de28 <col:3> 'class ClassB *' this
  | `-ImplicitCastExpr 0x403def0 </data/rosetta/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> 'pointer':'class ClassA *' <NullToPointer>
  |   `-GNUNullExpr 0x403de70 <col:18> 'long'

  `-CXXOperatorCallExpr 0x403e108 </data/rosetta/tools/clang_ast_transform/test-access_ptr.cc:43:3, col:8> 'class utility::pointer::owning_ptr<class ClassA>' lvalue
    |-ImplicitCastExpr 0x403e0f0 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(*)(pointer)' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x403e0c8 <col:6> 'class utility::pointer::owning_ptr<class ClassA> &(pointer)' lvalue CXXMethod 0x3febb50 'operator=' 'class utility::pointer::owning_ptr<class ClassA> &(pointer)'
    |-MemberExpr 0x403e000 <col:3> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' lvalue ->a_ 0x3fec7b0
    | `-CXXThisExpr 0x403dfe8 <col:3> 'class ClassB *' this
    `-ImplicitCastExpr 0x403e0b0 <col:8> 'pointer':'class ClassA *' <NullToPointer>
      `-IntegerLiteral 0x403e030 <col:8> 'int' 0
*/

RewriteAssignmentsOper RewriteAssignmentsOperCallback4(Replacements,
	"RewriteAssignmentsOper:nullPtrLiteral/integerLiteral");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				declRefExpr( isUtilityPointer() ).bind("castexpr")
			),
			has(
				memberExpr().bind("castTo")
			),
			anyOf(
				has(
					integerLiteral(equals(0)).bind("castFrom")
				),
				has(
					// This doesn't work for some reason yet
					expr(nullPtrLiteralExpr()).bind("castFrom")
				)
			),
			isUtilityPointer()
		)
	).bind("expr"),
	&RewriteAssignmentsOperCallback4);


#if 0
/*
		ClassAOP aa = new ClassA;

  `-DeclStmt 0x2d417b8 <line:69:3, col:27>
    `-VarDecl 0x2d41560 <col:3, col:21> col:12 aa 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>'
      `-ExprWithCleanups 0x2d417a0 <col:12, col:21> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>'
        `-CXXConstructExpr 0x2d41768 <col:12, col:21> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (const class utility::pointer::owning_ptr<class ClassA> &)' elidable 
          `-MaterializeTemporaryExpr 0x2d41748 <col:17, col:21> 'const class utility::pointer::owning_ptr<class ClassA>' lvalue
            `-ImplicitCastExpr 0x2d41730 <col:17, col:21> 'const class utility::pointer::owning_ptr<class ClassA>' <NoOp>
              `-CXXBindTemporaryExpr 0x2d416d8 <col:17, col:21> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x2d416d0)
                `-ImplicitCastExpr 0x2d416b8 <col:17, col:21> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' <ConstructorConversion>
                  `-CXXConstructExpr 0x2d41680 <col:17, col:21> 'ClassAOP':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
                    `-CXXNewExpr 0x2d415f8 <col:17, col:21> 'class ClassA *'
                      `-CXXConstructExpr 0x2d415c8 <col:21> 'class ClassA' 'void (void)'
*/

RewriteAssignmentsOper RewriteAssignmentsOperCallback5(Replacements,
	"RewriteAssignmentsOper:varDecl+newExpr");
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
						).bind("castexpr")
					)
				).bind("castTo")
			),
			hasParent(varDecl()),
			isUtilityPointer()			
		)
	).bind("expr"),
	&RewriteAssignmentsOperCallback5);

#endif

