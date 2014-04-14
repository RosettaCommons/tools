////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace calls to operator() on owning_ptrs
/*
CXXOperatorCallExpr 'pointer':'const class X *'
|-ImplicitCastExpr 'pointer (*)(void) const' <FunctionToPointerDecay>
| `-DeclRefExpr 'pointer (void) const' lvalue
`-ImplicitCastExpr 'const class utility::pointer::owning_ptr<X>' lvalue <NoOp>
	`-DeclRefExpr 'TokenCOP':'class utility::pointer::owning_ptr<X>' 
*/

class RewriteVoidPtrOperator : public ReplaceMatchCallback {
public:
	RewriteVoidPtrOperator(
			tooling::Replacements *Replace,
			const char *tag = "RewriteVoidPtrOperator") :
		ReplaceMatchCallback(Replace, tag) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;

		const Expr *expr = Result.Nodes.getStmtAs<Expr>("expr");
		const Expr *castFrom = Result.Nodes.getStmtAs<Expr>("castFrom");
		const Expr *castTo = Result.Nodes.getStmtAs<Expr>("castTo");
		const Expr *nullNode = Result.Nodes.getStmtAs<Expr>("null");
		
		if(!rewriteThisFile(expr, sm))
			return;

		// Get original code
		const std::string origCode = getText(sm, expr);
		if(origCode.empty())
			return;
		
		if(origCode.find('*') != std::string::npos) {
			// Hack: try to avoid *e1() etc.
			return;
		}

		bool deref = true;
		const std::string castFromType(
			castFrom ? QualType::getAsString( castFrom->getType().getSplitDesugaredType() ) : ""
		);
		const std::string castToType( 
			castTo ? QualType::getAsString( castTo  ->getType().getSplitDesugaredType() ) : ""
		);
		const std::string nullStr( nullNode ? getText(sm, nullNode) : ""); 
		
		llvm::errs() << origCode << "\n";	
		llvm::errs() << "castFrom: " << castFromType << "\n";
		llvm::errs() << "castTo: " << castToType << "\n";
		
		if(castFromType == castToType)
			return;
			
		if(
			// i.e.: data() == NULL
			(nullStr == "NULL" || nullStr == "__null" || nullStr == "0") || 
			// Cast between AP and OP, pass pointer directly
			(checkIsUtilityPointer(castFromType) && checkIsUtilityPointer(castToType))
		) {
			deref = false;
		}
						
		// origCode should end with (), so strip that call operator
		std::string newCode = std::string(origCode, 0, origCode.length() -2);
		if(deref)
			newCode = "&(*" + newCode + ")";
			
		doRewrite(sm, expr, origCode, newCode);
	}
};

// Round 1
// (anyOf handles max 5 statements)

// The below matcher set covers these AST's

/*
CXXConstructExpr 0x7fcd48e1a460 <col:10, col:26> 'AtomOP':'class utility::pointer::owning_ptr<class core::kinematics::tree::Atom>' 'void (pointer)'
`-CXXOperatorCallExpr 0x7fcd48e1a3a0 <col:18, col:26> 'pointer':'class core::kinematics::tree::Atom *'
  |-ImplicitCastExpr 0x7fcd48e1a388 <col:25, col:26> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x7fcd48e1a360 <col:25, col:26> 'pointer (void) const' lvalue CXXMethod 0x7fcd48defcc0 'operator()' 'pointer (void) const'
  `-ImplicitCastExpr 0x7fcd48e1a3e0 <col:18> 'const class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' lvalue <NoOp>
    `-MemberExpr 0x7fcd48e1a330 <col:18> 'AtomAP':'class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' lvalue ->parent_ 0x7fcd48df0380

CXXConstructExpr 0x7fcd48e19a58 <col:10, col:27> 'AtomCOP':'class utility::pointer::owning_ptr<const class core::kinematics::tree::Atom>' 'void (pointer)'
`-ImplicitCastExpr 0x7fcd48e19a40 <col:19, col:27> 'pointer':'const class core::kinematics::tree::Atom *' <NoOp>
  `-CXXOperatorCallExpr 0x7fcd48e19940 <col:19, col:27> 'pointer':'class core::kinematics::tree::Atom *'
    |-ImplicitCastExpr 0x7fcd48e19928 <col:26, col:27> 'pointer (*)(void) const' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x7fcd48e198a0 <col:26, col:27> 'pointer (void) const' lvalue CXXMethod 0x7fcd48defcc0 'operator()' 'pointer (void) const'
    `-MemberExpr 0x7fcd48e19870 <col:19> 'const AtomAP':'const class utility::pointer::access_ptr<class core::kinematics::tree::Atom>' lvalue ->parent_ 0x7fcd48df0380
*/


RewriteVoidPtrOperator RewriteVoidPtrOperatorCallback1(Replacements,
	"RewriteVoidPtrOperator:constructExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			// CHILD EXPR: operator()
			has(
				declRefExpr( isVoidPtrOperator() )
			),
			// CHILD EXPR: castFrom
			anyOf(
				has(
					memberExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					bindTemporaryExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					declRefExpr( isUtilityPointer() ).bind("castFrom")
				)
			),
			// PARENT Expr/Stmt: castTo
			// could have used hasAncestor() here to cover both cases
			// (with and without implicit cast) but this is more strict/specific
			anyOf(
				hasParent(
					constructExpr( isUtilityPointer() ).bind("castTo")
				),
				hasParent(
					implicitCastExpr(
						hasParent(
							constructExpr( isUtilityPointer() ).bind("castTo")
						)
					)
				)
			)
		)
	).bind("expr"),
	&RewriteVoidPtrOperatorCallback1);


// Round 2

/*
CXXDynamicCastExpr 0x7f5caef5af28 <col:22, col:63> 'WrappedReal *' dynamic_cast<WrappedReal *> <Dynamic>
`-CXXOperatorCallExpr 0x7f5caef59de0 <col:54, col:61> 'pointer':'class utility::pointer::ReferenceCount *'
  |-ImplicitCastExpr 0x7f5caef59dc8 <col:60, col:61> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x7f5caef59d48 <col:60, col:61> 'pointer (void) const' lvalue CXXMethod 0x7f5caf1f97e0 'operator()' 'pointer (void) const'
  `-ImplicitCastExpr 0x7f5caef59e20 <col:54, col:59> 'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' <NoOp>
    `-CXXBindTemporaryExpr 0x7f5caef59d28 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' (CXXTemporary 0x7f5caef59d20)
      `-CXXMemberCallExpr 0x7f5caef59cf0 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>'
        `-MemberExpr 0x7f5caef59cc0 <col:54> '<bound member function type>' ->data 0x7f5caef8a990

ConditionalOperator 0x48c49e0 <col:8, col:31> 'const char *'
`-ImplicitCastExpr 0x48c4998 <col:8, col:18> '_Bool' <PointerToBoolean>
  `-CXXOperatorCallExpr 0x48c48d0 <col:8, col:18> 'pointer':'class utility::pointer::
    |-ImplicitCastExpr 0x48c48b8 <col:17, col:18> 'pointer (*)(void) const' <FunctionToPointerDecay>
    | `-DeclRefExpr 0x48c4838 <col:17, col:18> 'pointer (void) const' lvalue CXXMethod 0x48990f0 'operator()' 'pointer (void) const'
    `-MemberExpr 0x48c4808 <col:8, col:11> 'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>':'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' lvalue ->second 0x489db20
      `-CXXOperatorCallExpr 0x48c47c8 <col:8> 'pointer':'const struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > *'
        |-ImplicitCastExpr 0x48c47b0 <col:9> 'pointer (*)(void) const' <FunctionToPointerDecay>
        | `-DeclRefExpr 0x48c4788 <col:9> 'pointer (void) const' lvalue CXXMethod 0x489afc0 'operator->' 'pointer (void) const'
        `-ImplicitCastExpr 0x48c4770 <col:8> 'const struct std::_Rb_tree_const_iterator<struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > >' lvalue <NoOp>
          `-DeclRefExpr 0x48c4748 <col:8> 'class ResourceManager::ResourcesMap::const_iterator':'struct std::_Rb_tree_const_iterator<struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > >' lvalue Var 0x48c0dc0 'r' 'class ResourceManager::ResourcesMap::const_iterator':'struct std::_Rb_tree_const_iterator<struct std::pair<const class std::basic_string<char>, class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount> > >'


BinaryOperator 0x1ee6ec0 <col:12, /data/rosetta/clang/build/bin/../lib/clang/3.5.0/include/stddef.h:72:18> '_Bool' '=='
`-CXXOperatorCallExpr 0x1ee6e50 </data/rosetta/main/source/src/utility/signals/Link.hh:108:12, col:18> 'pointer':'struct utility::signals::LinkUnit *'
  |-ImplicitCastExpr 0x1ee6e38 <col:17, col:18> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x1ee6e10 <col:17, col:18> 'pointer (void) const' lvalue CXXMethod 0x1ee38d0 'operator()' 'pointer (void) const'
  `-MemberExpr 0x1ee6de0 <col:12> 'const LinkUnitOP':'const class utility::pointer::owning_ptr<struct utility::signals::LinkUnit>' lvalue ->unit_ 0x1ee5be0
    `-CXXThisExpr 0x1ee6dc8 <col:12> 'const class utility::signals::Link *' this
*/

RewriteVoidPtrOperator RewriteVoidPtrOperatorCallback2(Replacements,
	"RewriteVoidPtrOperator:misc");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			// CHILD EXPR: operator()
			has(
				declRefExpr( isVoidPtrOperator() )
			),
			// CHILD EXPR: castFrom
			anyOf(
				has(
					memberExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					bindTemporaryExpr( isUtilityPointer() ).bind("castFrom")
				)
			),
			// PARENT Expr/Stmt: castTo
			anyOf(
				hasParent(
					dynamicCastExpr().bind("castTo")
				),
				hasParent(
					implicitCastExpr(
						hasParent(
							conditionalOperator()
						)
					).bind("castTo")
				),
				hasParent(
					binaryOperator(
						hasDirect(
							implicitCastExpr().bind("null")
						)
					).bind("castTo")
				),
				hasParent(
					binaryOperator().bind("castTo")
				)
			)
		)
	).bind("expr"),
	&RewriteVoidPtrOperatorCallback2);
