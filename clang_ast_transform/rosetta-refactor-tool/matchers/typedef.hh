/*
	Replace owning_ptr/access_ptr in declarations

	Type definitions:
		typedef utility::pointer::access_ptr< Some > SomeAP
		typedef utility::pointer::owning_ptr< Some > SomeOP
		typedef utility::pointer::access_ptr< Some const > SomeCAP
		typedef utility::pointer::owning_ptr< Some const > SomeCOP

	Parameter Declarations:
		void foo( utility::vector1< utility::pointer::owning_ptr<class ClassA> > & ops_ );
		void foo( utility::vector1< utility::pointer::access_ptr<class ClassA> > & aps_ );

	Field Declarations (in class header):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > ops_;
		utility::vector1< utility::pointer::access_ptr<class ClassA> > aps_;

	Function Declarations (return type):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > foo();
		utility::vector1< utility::pointer::access_ptr<class ClassA> > foo();

	Method Declarations (return type):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > foo();
		utility::vector1< utility::pointer::access_ptr<class ClassA> > foo();

	Variable Declarations (in method):
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > ops_;
		utility::vector1< utility::pointer::access_ptr<class ClassA> > aps_;
*/

class RewriteDecl : public ReplaceMatchCallback {
	
public:
	RewriteDecl(
		tooling::Replacements *Replace,
		const char *tag ="Decl",
		const char *delim ="") :
		ReplaceMatchCallback(Replace, tag),
		delim(delim) {}
	
	virtual void run(const ast_matchers::MatchFinder::MatchResult &Result) {
		SourceManager &sm = *Result.SourceManager;
		const Decl *node = Result.Nodes.getStmtAs<Decl>("decl");

		if(!rewriteThisFile(node, sm))
			return;
		
		const std::string origCode( getText(sm, node) );
		std::string newCode( origCode );
		std::string suffix;

		if(delim && *delim) {
			size_t p = newCode.find(delim);
			if(p != std::string::npos) {
				suffix = newCode.substr(p);
				newCode = newCode.substr(0, p);
			}
		}
		
		replace(newCode, "owning_ptr", "shared_ptr");
		replace(newCode, "access_ptr", "weak_ptr");

		doRewrite(sm, node, origCode, newCode + suffix);
	}
	
private:
	const char *delim;
};


// Typedefs for access_ptr and owning_ptr
RewriteDecl RewriteDeclCallback1(Replacements, "Decl:isTypedefDecl");
Finder.addMatcher(
	decl( isTypedefDecl() ).bind("decl"),
	&RewriteDeclCallback1);

// Parameter declaration in methods
RewriteDecl RewriteDeclCallback2(Replacements, "Decl:paramVarDecl");
Finder.addMatcher(
	parmVarDecl().bind("decl"),
	&RewriteDeclCallback2);

// Field (variable) declaration in classes
RewriteDecl RewriteDeclCallback3(Replacements, "Decl:fieldDecl");
Finder.addMatcher(
	fieldDecl().bind("decl"),
	&RewriteDeclCallback3);

// Method return type declaration
RewriteDecl RewriteDeclCallback4(Replacements, "Decl:methodDecl", "{");
Finder.addMatcher(
	methodDecl().bind("decl"),
	&RewriteDeclCallback4);

// Function return type declaration
RewriteDecl RewriteDeclCallback5(Replacements, "Decl:functionDecl", "{");
Finder.addMatcher(
	functionDecl().bind("decl"),
	&RewriteDeclCallback5);

// Local method varialble declaration
RewriteDecl RewriteDeclCallback6(Replacements, "Decl:varDecl");
Finder.addMatcher(
	varDecl().bind("decl"),
	&RewriteDeclCallback6);