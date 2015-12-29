// Declares clang::SyntaxOnlyAction.
#include "clang/Frontend/FrontendActions.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"

#include <clang/AST/RecursiveASTVisitor.h>
#include <clang/Basic/SourceLocation.h>
#include <clang/Frontend/CompilerInstance.h>
#include <clang/AST/Comment.h>

// Declares llvm::cl::extrahelp.
#include "llvm/Support/CommandLine.h"


#include <context.hpp>
#include <function.hpp>
#include <class.hpp>

using namespace clang::tooling;
using namespace llvm;

// Apply a custom category to all command-line options so that they are the
// only ones displayed.
static llvm::cl::OptionCategory BinderToolCategory("binder options");

// CommonOptionsParser declares HelpMessage with a description of the common
// command-line options related to the compilation database and input files.
// It's nice to have this help message in all tools.
static cl::extrahelp CommonHelp(CommonOptionsParser::HelpMessage);

// A help message for this specific tool can be added afterwards.
static cl::extrahelp MoreHelp("\nMore help text...\n");


using namespace clang;

using std::string;

// using namespace clang::driver;
// using namespace clang::tooling;
// using namespace llvm;





class ClassVisitor : public RecursiveASTVisitor<ClassVisitor>
{
public:
    explicit ClassVisitor(DeclContext *dc) : decl_context(dc) {}


	virtual bool VisitEnumDecl(EnumDecl *record) {
		errs() << "ClassVisitor EnumDecl: " << record->getQualifiedNameAsString() << "\n";
		record->dump();
        return true;
	}

private:
    DeclContext *decl_context;
};

string wrap_CXXRecordDecl(CXXRecordDecl *R)
{
	// for(auto c = R->captures_begin(); c != R->captures_end(); ++c) {
	// 	c->dump();
	// }

	// for(auto c = R->decls_begin(); c != R->decls_end(); ++c) {
	// 	c->dump();
	// }

	ClassVisitor v{R};
	v.TraverseDecl( R );

	//R->dump();

	return "";
}


class BinderVisitor : public RecursiveASTVisitor<BinderVisitor>
{
public:
    explicit BinderVisitor(CompilerInstance *ci) : ast_context( &( ci->getASTContext() ) ) {}


	virtual bool VisitFunctionDecl(FunctionDecl *record)
	{
		if( record->isCXXInstanceMember() ) return true;
		if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;


		binder::Item I { binder::bind_function(binder::_module_variable_name_, record) };
		//outs() << I;
		context.add(I);

        return true;
    }

	virtual bool VisitCXXRecordDecl(CXXRecordDecl *record) {
		if( record->isCXXInstanceMember() ) return true;
		if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;

		binder::Item I { binder::bind_class(record) };
		context.add(I);

        return true;
    }

	// virtual bool VisitClassTemplateSpecializationDecl(ClassTemplateSpecializationDecl *record) {
	// 	if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;

	// 	errs() << "Visit ClassTemplateSpecializationDecl:" << record->getQualifiedNameAsString() << "\n";
	// 	record->dump();
    //     return true;
	// }

	// virtual bool VisitTemplateDecl(TemplateDecl *record) {
	// 	//if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;
	// 	errs() << "Visit TemplateDecl: " << record->getQualifiedNameAsString() << "\n";
	// 	record->dump();
    //     return true;
	// }


	// virtual bool VisitTypedefDecl(TypedefDecl *record) {
	// 	if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;

	// 	errs() << "Visit TypedefDecl: " << record->getQualifiedNameAsString() << "\n";
	// 	// record->dump();
    //     return true;
	// }


	// virtual bool VisitFieldDecl(FieldDecl *record) {
	// 	errs() << "Visit FieldDecl: " << record->getQualifiedNameAsString() << "\n";
	// 	record->dump();
    //     return true;
	// }

	// virtual bool VisitEnumDecl(EnumDecl *record) {
	// 	errs() << "GlobalVisitor EnumDecl: " << record->getQualifiedNameAsString() << " isCXXClassMember:" << record->isCXXClassMember() << " isCXXInstanceMember:" << record->isCXXInstanceMember() << " isExternallyVisible:" << record->isExternallyVisible() << "\n";
	// 	record->dump();
    //     return true;
	// }

	// virtual bool VisitNamedDecl(NamedDecl *record) {
	// 	errs() << "Visit NamedRecord: " << record->getQualifiedNameAsString() << "\n";
    //     return true;
	// }


	void generate(void) { context.generate(); }
private:
    ASTContext *ast_context;

	binder::Context context;
};


class BinderASTConsumer : public ASTConsumer
{
private:
    std::unique_ptr<BinderVisitor> visitor;

public:
    // override the constructor in order to pass CI
    explicit BinderASTConsumer(CompilerInstance *ci) : visitor(new BinderVisitor(ci)) {}

    // override this to call our ExampleVisitor on the entire source file
    virtual void HandleTranslationUnit(ASTContext &context)
	{
        visitor->TraverseDecl( context.getTranslationUnitDecl() );
		visitor->generate();
    }
};


class BinderFrontendAction : public ASTFrontendAction {
public:
    virtual std::unique_ptr<clang::ASTConsumer> CreateASTConsumer(CompilerInstance &ci, StringRef file) {
        return std::unique_ptr<ASTConsumer>( new BinderASTConsumer(&ci) );
    }
};


int main(int argc, const char **argv)
{
  CommonOptionsParser op(argc, argv, BinderToolCategory);

  ClangTool tool(op.getCompilations(), op.getSourcePathList());

  return tool.run(newFrontendActionFactory<BinderFrontendAction>().get());
}



/*
VisitFunctionDecl
        //string funcName = func->getNameInfo().getName().getAsString();
        //string funcName = func->getNameInfo().getAsString();

		// QualType qt = record->getReturnType();
		// string rqtype = qt.getAsString();
		// //string rqtype = QualType::getAsString( qt.split() );

		// string funcName = record->getQualifiedNameAsString();


        // //if (funcName == "do_math") {
        // //    rewriter.ReplaceText(func->getLocation(), funcName.length(), "add5");
		// errs() << "Visit function def: " << rqtype << ' ' << funcName << " isCXXClassMember:" << record->isCXXClassMember() << " isCXXInstanceMember:" << record->isCXXInstanceMember() << " isExternallyVisible:" << record->isExternallyVisible() << " isDefined:" << record->isDefined() << "\n";

		// // errs() << "  name: " << record->getNameAsString() << "\n";
		// // errs() << "  Return Type (getCanonicalType): " << qt.getCanonicalType().getAsString() << "\n";
		// // errs() << "  Return Type (getCanonicalTypeInternal): " << qt->getCanonicalTypeInternal().getAsString() << "\n";
		// // errs() << "  Return Type (getCanonicalTypeUnqualified): " << QualType( qt->getCanonicalTypeUnqualified() ).getAsString() << "\n";
		//errs() << "  Source location:" << record->getLocation().printToString( ast_context->getSourceManager() ) << "\n";
		//errs() << "  Source location: " << ast_context->getSourceManager().getFilename( record->getLocation() )  << "\n";
		// errs() << "  Source location: " << ast_context->getSourceManager().getFileEntryForID( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).getFileID() )->getName() << "\n";

		// if( auto comment = ast_context->getLocalCommentForDeclUncached(record) ) comment->dumpColor();


		// for(uint i=0; i<record->getNumParams(); ++i) {
		// 	errs() << "  param[" << i << "]: " << record->getParamDecl(i)->getName() << ' ' << record->getParamDecl(i)->getOriginalType().getAsString() << "\n";
		// }

		// //record->getReturnType().dump();
        // //}
		// errs() << "\n";



virtual bool VisitRecordDecl(RecordDecl *record) {
virtual bool VisitCXXRecordDecl(CXXRecordDecl *record) {
	// 	if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;

    //     //string record_name = record->getNameInfo().getName().getAsString();
	// 	//errs() << "Visit record def: " << record->getNameAsString() << "\n";
	// 	errs() << "Visit VisitCXXRecordDecl:" << record->getQualifiedNameAsString() << " isCXXClassMember:" << record->isCXXClassMember() << " isCXXInstanceMember:" << record->isCXXInstanceMember() << " isExternallyVisible:" << record->isExternallyVisible() << "\n";

	// 	//errs() << "                  " << record->getTypeSourceInfo()->getType()->getAsString() << "\n";
	// 	//errs() << "                  " << record->getTypeForDecl()->getTypeClassName() << "\n";
	// 	//errs() << "                  " << record->getDeclName().getAsString() << "\n";

	// 	errs() << wrap_CXXRecordDecl(record);

    //     return true;
	// }



		*/
