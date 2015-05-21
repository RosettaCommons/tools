////////////////////////////////////////////////////////////////////////////////////////////////////
// Our custom clang AST matchers

namespace clang {
namespace ast_matchers {

AST_MATCHER(Expr, trueExpr) {
  return true;
}

AST_MATCHER(Decl, trueDecl) {
  return true;
}

AST_MATCHER(CastExpr, isNullToPointer) {
  return Node.getCastKind() == CK_NullToPointer ||
    Node.getCastKind() == CK_NullToMemberPointer;
}

AST_MATCHER(Expr, isNullExpr) {
  return Node.getStmtClass() == AbstractConditionalOperator::GNUNullExprClass;
}

AST_MATCHER(CastExpr, isNonNoopCast) {
  return Node.getCastKind() != CK_NoOp;
}

AST_MATCHER(CastExpr, isLValueToRValueCast) {
  return Node.getCastKind() == CK_LValueToRValue;
}

AST_MATCHER(CastExpr, isFunctionToPointerDecayCast) {
  return Node.getCastKind() == CK_FunctionToPointerDecay;
}

AST_MATCHER(CastExpr, isConstructorConversionCast) {
  return Node.getCastKind() == CK_ConstructorConversion;
}

AST_MATCHER(CastExpr, isNullToPointerCast) {
  return Node.getCastKind() == CK_NullToPointer;
}

AST_MATCHER(DeclRefExpr, isCallOperator) {
	return Node.getNameInfo().getName().getAsString() == "operator()";
}

AST_MATCHER(DeclRefExpr, isNotOperator) {
	llvm::errs() << Node.getNameInfo().getName().getAsString() << "\n";
	return Node.getNameInfo().getName().getAsString() == "operator!";
}

AST_MATCHER(DeclRefExpr, isClassOperator) {
	return checkIsClassOperator( Node.getNameInfo().getName().getAsString() );
}

AST_MATCHER(DeclRefExpr, isNotClassOperator) {
	return !checkIsClassOperator( Node.getNameInfo().getName().getAsString() );
}

AST_MATCHER(Decl, isTypedefDecl) {
  //return (Node.getKind() == NK_Typedef);
  return !strcmp(Node.getDeclKindName(), "Typedef");
}

AST_MATCHER(Expr, isUtilityPointer) {

	// Sugared type
	QualType T = Node.getType();
	SplitQualType T_split = T.split();
	if(checkIsUtilityPointer(QualType::getAsString(T_split)))
		return true;
			
	if(!T.isNull()) {
		// Desugared type
		SplitQualType D_split = T.getSplitDesugaredType();
		if (T_split != D_split) {
			if(checkIsUtilityPointer(QualType::getAsString(D_split)))
				return true;
		}
	}

	return false;
}

AST_MATCHER(Expr, isUtilityAccessPointer) {

	// Sugared type
	QualType T = Node.getType();
	SplitQualType T_split = T.split();
	if(checkIsUtilityAccessPointer(QualType::getAsString(T_split)))
		return true;
			
	if(!T.isNull()) {
		// Desugared type
		SplitQualType D_split = T.getSplitDesugaredType();
		if (T_split != D_split) {
			if(checkIsUtilityAccessPointer(QualType::getAsString(D_split)))
				return true;
		}
	}

	return false;
}

AST_MATCHER(Expr, isUtilityOwningPointer) {

	// Sugared type
	QualType T = Node.getType();
	SplitQualType T_split = T.split();
	if(checkIsUtilityOwningPointer(QualType::getAsString(T_split)))
		return true;

	if(!T.isNull()) {
		// Desugared type
		SplitQualType D_split = T.getSplitDesugaredType();
		if (T_split != D_split) {
			if(checkIsUtilityOwningPointer(QualType::getAsString(D_split)))
				return true;
		}
	}

	return false;
}

AST_MATCHER(Expr, containsUtilityPointer) {

	// Sugared type
	QualType T = Node.getType();
	SplitQualType T_split = T.split();
	if(checkContainsUtilityPointer(QualType::getAsString(T_split)))
		return true;
			
	if(!T.isNull()) {
		// Desugared type
		SplitQualType D_split = T.getSplitDesugaredType();
		if (T_split != D_split) {
			if(checkContainsUtilityPointer(QualType::getAsString(D_split)))
				return true;
		}
	}

	return false;
}

AST_MATCHER(Type, sugaredNullptrType) {
	const Type *DesugaredType = Node.getUnqualifiedDesugaredType();
	if (const BuiltinType *BT = dyn_cast<BuiltinType>(DesugaredType))
		return BT->getKind() == BuiltinType::NullPtr;
	return false;
}

AST_MATCHER(QualType, isReal) {
	return Node->isRealFloatingType();
}

} // end namespace ast_matchers
} // end namespace clang
