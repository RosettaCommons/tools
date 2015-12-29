#include <string>
#include <iostream>

// std::string foo(int a) { return "foo1"; }
// std::string foo(int a, int b) { return "foo2"; }
// std::string foo(int a, int b, int c) { return "foo3"; }


// std::string ffoo(int ) { return "foo1"; }
// std::string ffoo(int , int ) { return "foo2"; }
// std::string ffoo(int , int , int ) { return "foo3"; }






enum class Global_Enum {G_E1, G_E2, G_E3};




namespace aaaa {

class A
{
public:

	A() {}
	A(int) {}
	A(int, int) {}

	int foo() { return 42; }
	std::string foo(int a) { return "foo:" + std::to_string(a); }
	std::string foo(int a, std::string const &b, int c) { return "foo: " + b + " "+ std::to_string(c); }
};

struct S
{
	S() {}

	int foo() { return 43; }
};


void hi_(void)
{
	std::cout << "Hi in aaaa..." << std::endl;
}


};

void hi(void)
{
	std::cout << "Hi..." << std::endl;
}
