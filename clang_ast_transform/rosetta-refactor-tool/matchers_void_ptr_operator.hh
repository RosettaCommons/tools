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

		const std::string castFromType( QualType::getAsString( castFrom->getType().getSplitDesugaredType() ) );
		const std::string castToType  ( QualType::getAsString( castTo  ->getType().getSplitDesugaredType() ) );
		bool deref = true;
		
		llvm::errs() << origCode << "\n";	
		llvm::errs() << "castFrom: " << castFromType << "\n";
		llvm::errs() << "castTo: " << castToType << "\n";
		
		if(castFromType == castToType)
			return;
			
		if(checkIsUtilityPointer(castFromType) && checkIsUtilityPointer(castToType)) {
			// Cast between AP and OP (whatever direction), can pass directly
			deref = false;
		}
				
		// origCode should end with (), so strip that call operator
		std::string newCode = std::string(origCode, 0, origCode.length() -2);
		if(deref)
			newCode = "&(*" + newCode + ")";
			
		doRewrite(sm, expr, origCode, newCode);
	}
};

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

CXXDynamicCastExpr 0x7f5caef5af28 <col:22, col:63> 'WrappedReal *' dynamic_cast<WrappedReal *> <Dynamic>
`-CXXOperatorCallExpr 0x7f5caef59de0 <col:54, col:61> 'pointer':'class utility::pointer::ReferenceCount *'
  |-ImplicitCastExpr 0x7f5caef59dc8 <col:60, col:61> 'pointer (*)(void) const' <FunctionToPointerDecay>
  | `-DeclRefExpr 0x7f5caef59d48 <col:60, col:61> 'pointer (void) const' lvalue CXXMethod 0x7f5caf1f97e0 'operator()' 'pointer (void) const'
  `-ImplicitCastExpr 0x7f5caef59e20 <col:54, col:59> 'const class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' <NoOp>
    `-CXXBindTemporaryExpr 0x7f5caef59d28 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>' (CXXTemporary 0x7f5caef59d20)
      `-CXXMemberCallExpr 0x7f5caef59cf0 <col:54, col:59> 'utility::pointer::ReferenceCountOP':'class utility::pointer::owning_ptr<class utility::pointer::ReferenceCount>'
        `-MemberExpr 0x7f5caef59cc0 <col:54> '<bound member function type>' ->data 0x7f5caef8a990
*/

RewriteVoidPtrOperator RewriteVoidPtrOperatorCallback1(Replacements,
	"RewriteVoidPtrOperator:constructExpr");
Finder.addMatcher(
	operatorCallExpr(
		allOf(
			has(
				declRefExpr( isVoidPtrOperator() )
			),
			anyOf(
				has(
					memberExpr( isUtilityPointer() ).bind("castFrom")
				),
				has(
					bindTemporaryExpr( isUtilityPointer() ).bind("castFrom")
				)
			),
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
				),
				hasParent(
					dynamicCastExpr( isUtilityPointer() ).bind("castTo")
				)
			)
		)
	).bind("expr"),
	&RewriteVoidPtrOperatorCallback1);
