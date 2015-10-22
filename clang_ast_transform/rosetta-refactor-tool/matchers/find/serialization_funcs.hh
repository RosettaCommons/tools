/*
	Find records with save and load functions and also the member variables
	that are mentioned inside of these functions.
*/

class SerializationFuncFinder : public ReplaceMatchCallback {

public:
	SerializationFuncFinder(tooling::Replacements *Replace) :
		ReplaceMatchCallback(Replace, "RecordDeclFinder")
		{}

	// Main callback for all matches
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {

		SourceManager *sm = Result.SourceManager;
		const CXXOperatorCallExpr * opcall = Result.Nodes.getStmtAs<CXXOperatorCallExpr>("op_call");
		const MemberExpr * member_var = Result.Nodes.getStmtAs<MemberExpr>("member");
		//const CXXMethodDecl * save_method = Result.Nodes.getStmtAs<CXXMethodDecl>("savemethod");
		const FunctionTemplateDecl * save_method = Result.Nodes.getStmtAs<FunctionTemplateDecl>("savemethod");

		std::cout << "opcall " << opcall << " member_var " << member_var << " save_method " << save_method << std::endl;

		if(!rewriteThisFile(opcall, *sm))
			return;
		//std::cout << "found one" << std::endl; opcall->dump(); std::cout << "\n"; member_var->dump(); std::cout << "\n" << std::endl;
		std::cout << "Save method:" << save_method << " " << member_var << " " << opcall << " ";
		save_method->dump();
		
		//std::cout << std::endl;
		//std::cout << save_method->getThisType( *Result.Context )->getPointeeType().getCanonicalType().getUnqualifiedType().getAsString();
		// 
		std::cout << " " << member_var->getMemberNameInfo().getAsString();
		std::cout << std::endl;
		
	}
};

//   | |-CXXMethodDecl 0x4aab890 <line:11:3, line:14:3> line:11:8 save 'void (class cereal::BinaryOutputArchive &) const'
//   | | |-TemplateArgument type 'class cereal::BinaryOutputArchive'
//   | | |-ParmVarDecl 0x4aab7d0 <col:14, col:24> col:24 arc 'class cereal::BinaryOutputArchive &'
//   | | `-CompoundStmt 0x4aac1b8 <col:36, line:14:3>
//   | |   |-CXXOperatorCallExpr 0x4aabd40 <line:12:5, col:17> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue
//   | |   | |-ImplicitCastExpr 0x4aabd28 <col:5, col:17> 'class cereal::BinaryOutputArchive &(*)(const int &&)' <FunctionToPointerDecay>
//   | |   | | `-DeclRefExpr 0x4aabca8 <col:5, col:17> 'class cereal::BinaryOutputArchive &(const int &&)' lvalue CXXMethod 0x4aabba0 'operator()' 'class cereal::BinaryOutputArchive &(const int &&)'
//   | |   | |-ImplicitCastExpr 0x4aabd88 <col:5> 'class cereal::OutputArchive<class cereal::BinaryOutputArchive, 1>' lvalue <UncheckedDerivedToBase (OutputArchive)>
//   | |   | | `-DeclRefExpr 0x4aab998 <col:5> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue ParmVar 0x4aab7d0 'arc' 'class cereal::BinaryOutputArchive &'
//   | |   | `-MemberExpr 0x4aab220 <col:10> 'const int' lvalue ->myint_ 0x4aaae80
//   | |   |   `-CXXThisExpr 0x4aab208 <col:10> 'const class MyClass *' this
//   | |   `-CXXOperatorCallExpr 0x4aac150 <line:13:5, col:19> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue
//   | |     |-ImplicitCastExpr 0x4aac138 <col:5, col:19> 'class cereal::BinaryOutputArchive &(*)(const float &&)' <FunctionToPointerDecay>
//   | |     | `-DeclRefExpr 0x4aac0b8 <col:5, col:19> 'class cereal::BinaryOutputArchive &(const float &&)' lvalue CXXMethod 0x4aabfb0 'operator()' 'class cereal::BinaryOutputArchive &(const float &&)'
//   | |     |-ImplicitCastExpr 0x4aac198 <col:5> 'class cereal::OutputArchive<class cereal::BinaryOutputArchive, 1>' lvalue <UncheckedDerivedToBase (OutputArchive)>
//   | |     | `-DeclRefExpr 0x4aabda8 <col:5> 'class cereal::BinaryOutputArchive':'class cereal::BinaryOutputArchive' lvalue ParmVar 0x4aab7d0 'arc' 'class cereal::BinaryOutputArchive &'
//   | |     `-MemberExpr 0x4aab2c0 <col:10> 'const float' lvalue ->myfloat_ 0x4aaaee0
//   | |       `-CXXThisExpr 0x4aab2a8 <col:10> 'const class MyClass *' this

// This works fine, but cannot report on multiple data members serialized within a single call operator, e.g. arc( myfloat_, mychar_ );
// Finder.addMatcher(
// 	operatorCallExpr(
// 		hasDescendant( memberExpr().bind("member")),
// 		hasAncestor( methodDecl(
// 			hasName( "save" )
// 			//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::BinaryOutputArchive")))))
// 			//hasTemplateArgument(0,refersToType(asString("class cereal::BinaryOutputArchive")))
// 		).bind("savemethod") )
// 	).bind( "op_call" ), new SerializationFuncFinder(Replacements) );

Finder.addMatcher(
	cxxOperatorCallExpr(
		hasDescendant( memberExpr().bind("member")),
		hasAncestor( functionTemplateDecl(
			hasName( "save" )
			//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::BinaryOutputArchive")))))
			,has( templateTypeParmDecl())
		).bind("savemethod") )
	).bind( "op_call" ), new SerializationFuncFinder(Replacements) );



// This for some reason reports each match multiple times! Grrr
//
// Finder.addMatcher(
// 	memberExpr( hasAncestor( operatorCallExpr(
// 		hasAncestor( methodDecl(
// 			hasName("save")
// 			//,hasParameter(0,hasType(referenceType(pointee(asString("class cereal::JSONOutputArchive")))))
// 		).bind("savemethod") )).bind("op_call"))
// 	).bind("member"),
// 	new SerializationFuncFinder(Replacements) );

//Finder.addMatcher(
//	recordDecl(
//		hasParent(
//			decl().bind("parent")
//		)
//	).bind("recorddecl"),
//	new RecordDeclFinder(Replacements));
	
