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
		const CXXOperatorCallExpr *opercallexpr = Result.Nodes.getStmtAs<CXXOperatorCallExpr>("opercallexpr");
		const CastExpr *castexpr = Result.Nodes.getStmtAs<CastExpr>("castexpr");
		const Expr *declrefexpr = Result.Nodes.getStmtAs<Expr>("expr");
		
		const FullSourceLoc FullLocation = FullSourceLoc(opercallexpr->getLocStart(), sm);
		if(FullLocation.getFileID() != sm.getMainFileID())
			return;

		// Get original code
		const std::string origCode = getText(sm, opercallexpr);
		if(origCode.empty())
			return;
			
		// Get cast type
		std::string type = QualType::getAsString(declrefexpr->getType().split());
		
		// Rewrite assignment	
		std::string leftSideCode = getTextToDelim(sm, opercallexpr, castexpr);
		std::string rightSideCode(origCode, leftSideCode.length());
		
		leftSideCode = trim(leftSideCode, " \t\n=");
		rightSideCode = trim(rightSideCode, " \t\n=");

		if(type == "<bound member function type>") {
			// Not sure what this is... internal type?
			return;
		}
		
		if(beginsWith(type, "std::map<")) {
			// Hack to use second argument of map definition as type
			size_t j = type.find(',');
			type = trim(std::string(type, j+1), "<>, ");
		}

		std::string newCode;
		if(rightSideCode == "0" || rightSideCode == "NULL") {
			//newCode = leftCode + " " + typedef_name + "()";
			newCode = leftSideCode + ".reset()";
		} else {
			newCode = leftSideCode + " = " + type + "( " + rightSideCode + " )";
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

RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback1(Replacements, "operCallExpr>memberExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				implicitCastExpr(isFunctionToPointerDecayCast()).bind("castexpr")
			),
			has(
				memberExpr(hasParent(operatorCallExpr())).bind("expr")
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

RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback2(Replacements, "operCallExpr>declRefExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				implicitCastExpr(isFunctionToPointerDecayCast()).bind("castexpr")
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

RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback3(Replacements, "operCallExpr>operCallExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
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
CXXOperatorCallExpr 0x2255200 <col:4, col:72> 'pointer':'class utility::pointer::ReferenceCount *'
|-ImplicitCastExpr 0x22551e8 <col:71, col:72> 'pointer (*)(void) const' <FunctionToPointerDecay>
| `-DeclRefExpr 0x2255168 <col:71, col:72> 'pointer (void) const' lvalue CXXMethod 0x2254b20 'operator()' 'pointer (void) const'
`-ImplicitCastExpr 0x2255268 <col:4, col:70> 'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' <NoOp>
  `-CXXBindTemporaryExpr 0x2255148 <col:4, col:70> 'ResourceOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' (CXXTemporary 0x2255140)
     `-CXXMemberCallExpr 0x2252900 <col:4, col:70> 'ResourceOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>'
       |-MemberExpr 0x22528a8 <col:4, col:37> '<bound member function type>' ->get_resource 0x1fc1d80
       | `-CallExpr 0x2252880 <col:4, col:34> 'class basic::resource_manager::ResourceManager *'
       |   `-ImplicitCastExpr 0x2252868 <col:4, col:21> 'class basic::resource_manager::ResourceManager *(*)(void)' <FunctionToPointerDecay>
       |     `-DeclRefExpr 0x2252808 <col:4, col:21> 'class basic::resource_manager::ResourceManager *(void)' lvalue CXXMethod 0x1fc0770 'get_instance' 'class basic::resource_manager::ResourceManager *(void)'
               `-DeclRefExpr 0x22528d8 <col:50> 'const ResourceDescription':'const class std::basic_string<char>' lvalue ParmVar 0x2252440 'resource_description' 'const ResourceDescription &'
*/

/*
// INCORRECT/NOT NEEDED
RewriteImplicitCastInAssignment RewriteImplicitCastInAssignmentCallback4(Replacements, "operCallExpr>>bindTemporaryExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				implicitCastExpr(isFunctionToPointerDecayCast()).bind("castexpr")
			),
			has(
				implicitCastExpr(
					has(
						bindTemporaryExpr( isUtilityPointer() ).bind("expr")
					)
				)
			)
		)
	).bind("opercallexpr"),
	&RewriteImplicitCastInAssignmentCallback4);
*/
