// Declares clang::SyntaxOnlyAction.
#include "clang/Frontend/FrontendActions.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"

#include <clang/AST/RecursiveASTVisitor.h>
#include <clang/Frontend/CompilerInstance.h>

// Declares llvm::cl::extrahelp.
#include "llvm/Support/CommandLine.h"

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


class BinderVisitor : public RecursiveASTVisitor<BinderVisitor>
{
public:
    explicit BinderVisitor(CompilerInstance *ci) : ast_context( &( ci->getASTContext() ) ) {}

	virtual bool VisitFunctionDecl(FunctionDecl *func)
	{
        string funcName = func->getNameInfo().getName().getAsString();
        //if (funcName == "do_math") {
        //    rewriter.ReplaceText(func->getLocation(), funcName.length(), "add5");
		errs() << "Visit function def: " << funcName << "\n";
        //}
        return true;
    }

	//virtual bool VisitRecordDecl(RecordDecl *record) {
	virtual bool VisitCXXRecordDecl(CXXRecordDecl *record) {

        //string record_name = record->getNameInfo().getName().getAsString();
		//errs() << "Visit record def: " << record->getNameAsString() << "\n";
		errs() << "Visit VisitCXXRecordDecl:" << record->getQualifiedNameAsString() << " isCXXClassMember:" << record->isCXXClassMember() << "\n";

		//errs() << "                  " << record->getTypeSourceInfo()->getType()->getAsString() << "\n";
		//errs() << "                  " << record->getTypeForDecl()->getTypeClassName() << "\n";
		//errs() << "                  " << record->getDeclName().getAsString() << "\n";
        return true;
	}

	virtual bool VisitClassTemplateSpecializationDecl(ClassTemplateSpecializationDecl *record) {
		errs() << "Visit ClassTemplateSpecializationDecl:" << record->getQualifiedNameAsString() << "\n";
        return true;
	}

	virtual bool VisitTemplateDecl(TemplateDecl *record) {
		errs() << "Visit TemplateDecl:" << record->getQualifiedNameAsString() << "\n";
		record->dump();
        return true;
	}



	// virtual bool VisitNamedDecl(NamedDecl *record) {
	// 	errs() << "Visit NamedRecord: " << record->getQualifiedNameAsString() << "\n";
    //     return true;
	// }



private:
    ASTContext *ast_context;
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
