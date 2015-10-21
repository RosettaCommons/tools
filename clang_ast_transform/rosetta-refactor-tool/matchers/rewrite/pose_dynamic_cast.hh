/*
	Replace dynamic casts in pose/conformation special cases:

	From:
		ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &( pose.conformation() ) );

	To:
		ProtectedConformationCOP conf = utility::pointer::dynamic_pointer_cast< ProtectedConformation const >( pose.conformation_ptr() );
*/

class RewritePoseDynamicCast : public ReplaceMatchCallback {
public:
	RewritePoseDynamicCast(
		tooling::Replacements *Replace,
		const char *replacementCastCode,
		const char *tag = "PoseDynamicCast") :
		ReplaceMatchCallback(Replace, tag),
		replacementCastCode(replacementCastCode) {}

	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;

		const Expr *dyncastexpr = Result.Nodes.getStmtAs<Expr>("dyncastexpr");
		const Expr *declrefexpr = Result.Nodes.getStmtAs<Expr>("declrefexpr");
		const Stmt *unaryoperator = Result.Nodes.getStmtAs<Stmt>("unaryoperator");
		const Expr *membercallexpr = Result.Nodes.getStmtAs<Expr>("membercallexpr");
		
		if(!rewriteThisFile(dyncastexpr, sm))
			return;

		// Get original code and cast type
		const std::string origCode = getText(sm, dyncastexpr);
		const std::string origCastType = QualType::getAsString( dyncastexpr->getType().split() );
		const std::string declRefType = QualType::getAsString( declrefexpr->getType().split() );
		const std::string origCastCode = getText(sm, unaryoperator);
		const std::string origCastCodeInside = getText(sm, membercallexpr);

		if(origCode.empty() || origCastCode.empty() || origCastType.empty())
			return;

		// We only are about Pose.conformation() here
		if(!endsWith(declRefType, "core::pose::Pose") || !endsWith(origCastCodeInside, ".conformation"))
			return;

		std::string castType = stripQualifiers(origCastType);
		if(endsWith(castType, "*"))
			castType = trim(std::string(castType, 0, castType.length() -1));
		if(origCastType.find("const ") != std::string::npos)
			castType += " const";	

		std::string newCastCode = origCastCodeInside;
		replace(newCastCode, ".conformation", ".conformation_ptr");

		std::string newCode = 
			replacementCastCode + "< " + castType + " > "
				+ "( " + newCastCode + "() )";
	
		doRewrite(sm, dyncastexpr, origCode, newCode);
	}
private:
	std::string replacementCastCode;
};


/*
  ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &( pose.conformation() ) );

  |-DeclStmt 0xa2675a0 <line:61:3, col:107>
  | `-VarDecl 0xa2643b0 <col:3, col:106> col:28 conf 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>'
  |   `-ExprWithCleanups 0xa267588 <col:28, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>'
  |     `-CXXConstructExpr 0xa267550 <col:28, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' 'void (const class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation> &)' elidable
  |       `-MaterializeTemporaryExpr 0xa267530 <col:35, col:106> 'const class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' lvalue
  |         `-ImplicitCastExpr 0xa267518 <col:35, col:106> 'const class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' <NoOp>
  |           `-CXXBindTemporaryExpr 0xa2671f8 <col:35, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' (CXXTemporary 0xa2671f0)
  |             `-ImplicitCastExpr 0xa2671d8 <col:35, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' <ConstructorConversion>
  |               `-CXXConstructExpr 0xa2671a0 <col:35, col:106> 'ProtectedConformationCOP':'class utility::pointer::owning_ptr<const class protocols::environment::ProtectedConformation>' 'void (pointer)'
  |                 `-CXXDynamicCastExpr 0xa264578 <col:35, col:106> 'const class protocols::environment::ProtectedConformation *' dynamic_cast<const class protocols::environment::ProtectedConformation *> <Dynamic>
  |                   `-UnaryOperator 0xa264548 <col:81, col:104> 'const Conformation *' prefix '&'
  |                     `-ParenExpr 0xa2644f8 <col:82, col:104> 'const Conformation':'const class core::conformation::Conformation' lvalue
  |                       `-CXXMemberCallExpr 0xa2644d0 <col:84, col:102> 'const Conformation':'const class core::conformation::Conformation' lvalue
  |                         `-MemberExpr 0xa2644a0 <col:84, col:89> '<bound member function type>' .conformation 0x97efd70
  |                           `-DeclRefExpr 0xa264408 <col:84> 'const core::pose::Pose':'const class core::pose::Pose' lvalue ParmVar 0xa263b80 'pose' 'const core::pose::Pose &'
*/

Finder.addMatcher(
	cxxDynamicCastExpr(
		hasParent(
			cxxConstructExpr( isUtilityPointer() ).bind("constructexpr")
		),
		has(
			unaryOperator(
				has(
					cxxMemberCallExpr(
						has(
							memberExpr(
								has(
									declRefExpr().bind("declrefexpr")
								)
							).bind("membercallexpr")
						)
					)
				)
			).bind("unaryoperator")
		)
	).bind("dyncastexpr"),
	new RewritePoseDynamicCast(Replacements, "utility::pointer::dynamic_pointer_cast", "DynamicCast:constructExpr"));
