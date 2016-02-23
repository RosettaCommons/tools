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
#include <util.hpp>

using namespace clang::tooling;
using namespace llvm;


// Apply a custom category to all command-line options so that they are the
// only ones displayed.
static llvm::cl::OptionCategory BinderToolCategory("Binder options");

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


cl::opt<std::string> O_root_module("root-module", cl::desc("Name of root module"), /*cl::init("example"),*/ cl::cat(BinderToolCategory));
//cl::opt<std::string> O_root_module("root-module", cl::desc("Name of root module"), cl::Required, cl::cat(BinderToolCategory));

cl::opt<int> O_max_file_size("max-file-size", cl::desc("Specify maximum length of generated source files"), cl::init(1024*4), cl::cat(BinderToolCategory));

cl::opt<std::string> O_prefix("prefix", cl::desc("Output prefix for all generated files. Might contain directories."), cl::init(""), cl::cat(BinderToolCategory));

cl::list<std::string> O_bind("bind", cl::desc("Namespace to bind, could be specified more then once. Specify \"\" to bind all namespaces.")); // , cl::OneOrMore
cl::list<std::string> O_skip("skip", cl::desc("Namespace to skip, could be specified more then once")); // , cl::OneOrMore


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


	bool shouldVisitTemplateInstantiations () const { return true; }


	virtual bool VisitFunctionDecl(FunctionDecl *F)
	{
		if( F->isCXXInstanceMember() or isa<CXXMethodDecl>(F) ) return true;

		if( binder::is_bindable(F) ) {
			binder::BinderOP b = std::make_shared<binder::FunctionBinder>(F);
			context.add(b);
		}

		// if( FullSourceLoc(F->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;
		// if( binder::is_bindable(F)  and  !binder::namespace_from_named_decl(F).compare(0, 7, "utility") ) {
		// 	outs() << "Adding function: " << F->getQualifiedNameAsString() << "\n";
		// 	binder::BinderOP b{ std::make_shared<binder::FunctionBinder>(F) };
		// 	context.add(b);
		// }
		// else {
		// 	//outs() << "Skipping " << record->getQualifiedNameAsString() << " because it is not bindable...\n";
		// }

        return true;
    }

	virtual bool VisitCXXRecordDecl(CXXRecordDecl *C) {
		if( C->isCXXInstanceMember()  or  C->isCXXClassMember() ) return true;

		if( binder::is_bindable(C) ) {
			binder::BinderOP b = std::make_shared<binder::ClassBinder>(C);
			context.add(b);
		}

		// if( FullSourceLoc(C->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;
		// if( isa<ClassTemplateSpecializationDecl>(C) ) {
		// 	errs() << "Visit VisitCXXRecordDecl->ClassTemplateSpecializationDecl:" << C->getQualifiedNameAsString() << binder::template_specialization(C) << "\n";
		// }
		// if( binder::is_bindable(C)  and  !binder::namespace_from_named_decl(C).compare(0, 7, "utility") ) {
		// 	outs() << "Adding class: " << C->getQualifiedNameAsString() << "\n";
		// 	binder::BinderOP b{ std::make_shared<binder::ClassBinder>(C) };
		// 	context.add(b);
		// }

        return true;
    }

	// virtual bool VisitClassTemplateSpecializationDecl(ClassTemplateSpecializationDecl *C) {
	// 	if( FullSourceLoc(C->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;
	// 	errs() << "Visit ClassTemplateSpecializationDecl:" << C->getQualifiedNameAsString() << binder::template_specialization(C) << "\n";
	// 	C->dump();
    //     return true;
	// }

	// virtual bool VisitTemplateDecl(TemplateDecl *record) {
	// 	//if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;
	// 	errs() << "Visit TemplateDecl: " << record->getQualifiedNameAsString() << "\n";
	// 	//record->dump();
    //     return true;
	// }

	// virtual bool VisitClassTemplateDecl(ClassTemplateDecl *record) {
	// 	//if( FullSourceLoc(record->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;
	// 	errs() << "Visit ClassTemplateDecl: " << record->getQualifiedNameAsString() << binder::template_specialization( record->getTemplatedDecl() ) << "\n";
	// 	//record->dump();
    //     return true;
	// }

	virtual bool VisitTypedefDecl(TypedefDecl *T) {
		if( FullSourceLoc(T->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;

		//errs() << "Visit TypedefDecl: " << T->getQualifiedNameAsString() << "  Type: " << T->getUnderlyingType()->getCanonicalTypeInternal()/*getCanonicalType()*/.getAsString() << "\n";
		// record->dump();
        return true;
	}


	// virtual bool VisitFieldDecl(FieldDecl *record) {
	// 	errs() << "Visit FieldDecl: " << record->getQualifiedNameAsString() << "\n";
	// 	record->dump();
    //     return true;
	// }

	virtual bool VisitEnumDecl(EnumDecl *E) {
		if( FullSourceLoc(E->getLocation(), ast_context->getSourceManager() ).isInSystemHeader() ) return true;
		//errs() << "Visitor EnumDecl: " << E->getQualifiedNameAsString() << " isCXXClassMember:" << E->isCXXClassMember() << " isCXXInstanceMember:" << E->isCXXInstanceMember() << " isExternallyVisible:" << E->isExternallyVisible() << "\n";
		//E->dump();
        return true;
	}

	// virtual bool VisitNamedDecl(NamedDecl *record) {
	// 	errs() << "Visit NamedRecord: " << record->getQualifiedNameAsString() << "\n";
    //     return true;
	// }


	void generate(void) {
		context.generate(O_root_module, O_bind, O_skip, O_prefix, O_max_file_size);
	}

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

	//outs() << "Root module: " << O_root_module << "\n";
	//for(auto &s : O_bind) outs() << "Binding: " << s << "\n";

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
