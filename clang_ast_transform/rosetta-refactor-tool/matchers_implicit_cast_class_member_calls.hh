////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace implicit casts in class member calls

class RewriteClassMemberCalls : public ReplaceMatchCallback {
public:
	RewriteClassMemberCalls(
			tooling::Replacements *Replace, 
			const char *tag = "RewriteClassMemberCalls") :
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

		std::string type(castToType);
		
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
		if(origCode == "0" || origCode == "NULL")
			newCode = type + "()";
		else
			newCode = type + "( " + origCode + " )";
		
		doRewrite(sm, expr, origCode, newCode);
	}
};


/*
	std::vector<ClassAOP> as_vector_;
	void set_a_vector1(ClassA *a) {
		as_vector_.push_back(a);
	}

  | `-CXXMemberCallExpr 0x2ef07d0 <col:3, col:25> 'void'
  |   |-MemberExpr 0x2ef07a0 <col:3, col:14> '<bound member function type>' .push_back 0x2e616a0
  |   | `-MemberExpr 0x2ef0670 <col:3> 'std::vector<ClassAOP>':'class std::vector<class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue ->as_vector_ 0x2e78910
  |   |   `-CXXThisExpr 0x2ef0658 <col:3> 'class ClassB *' this
  |   `-MaterializeTemporaryExpr 0x2ef08a8 <col:24> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' xvalue
  |     `-CXXBindTemporaryExpr 0x2ef0888 <col:24> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x2ef0880)
  |       `-CXXConstructExpr 0x2ef0848 <col:24> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
  |         `-ImplicitCastExpr 0x2ef0830 <col:24> 'class ClassA *' <LValueToRValue>
  |           `-DeclRefExpr 0x2ef0710 <col:24> 'class ClassA *' lvalue ParmVar 0x2ee9c70 'a' 'class ClassA *'
*/

RewriteClassMemberCalls RewriteClassMemberCallsCallback1(
	Replacements,
	"RewriteClassMemberCalls"
);

Finder.addMatcher(
	memberCallExpr(
		allOf(
			has(
				bindTemporaryExpr(
					has(
						constructExpr(
							has(
								declRefExpr().bind("castFrom")
							)
						)
					)
				).bind("expr")
			),
			has(
				memberExpr(
					has(
						memberExpr( containsUtilityPointer() ).bind("castTo")
					)
				)
			)
		)
	),
	&RewriteClassMemberCallsCallback1);


/*
	std::vector<ClassAOP> as_vector_;  
	void set_aref_vector1(ClassA & a) {
		as_vector_.push_back(&a);
	}

  | `-CXXMemberCallExpr 0x3e22090 <col:3, col:26> 'void'
  |   |-MemberExpr 0x3e22060 <col:3, col:14> '<bound member function type>' .push_back 0x3d7e850
  |   | `-MemberExpr 0x3e21f18 <col:3> 'std::vector<ClassAOP>':'class std::vector<class utility::pointer::owning_ptr<class ClassA>, class std::allocator<class utility::pointer::owning_ptr<class ClassA> > >' lvalue ->as_vector_ 0x3d95ac0
  |   |   `-CXXThisExpr 0x3e21f00 <col:3> 'class ClassB *' this
  |   `-MaterializeTemporaryExpr 0x3e22158 <col:24, col:25> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' xvalue
  |     `-CXXBindTemporaryExpr 0x3e22138 <col:24, col:25> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' (CXXTemporary 0x3e22130)
  |       `-CXXConstructExpr 0x3e220f0 <col:24, col:25> 'value_type':'class utility::pointer::owning_ptr<class ClassA>' 'void (pointer)'
  |         `-UnaryOperator 0x3e21fe0 <col:24, col:25> 'class ClassA *' prefix '&'
  |           `-DeclRefExpr 0x3e21fb8 <col:25> 'class ClassA' lvalue ParmVar 0x3e071a0 'a' 'class ClassA &'
*/

RewriteClassMemberCalls RewriteClassMemberCallsCallback2(
	Replacements,
	"RewriteClassMemberCalls:unaryOperator"
);
Finder.addMatcher(
	memberCallExpr(
		allOf(
			has(
				bindTemporaryExpr(
					has(
						constructExpr(
							has(
								unaryOperator(
									has(
										declRefExpr().bind("castFrom")
									)
								)
							)
						)
					)
				).bind("expr")
			),
			has(
				memberExpr(
					has(
						memberExpr( containsUtilityPointer() ).bind("castTo")
					)
				)
			)
		)
	),
	&RewriteClassMemberCallsCallback2);
