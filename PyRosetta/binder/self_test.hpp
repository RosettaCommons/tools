#include <string>
#include <iostream>

class A {
public:
	int n;
	float f;
	double d;
};

//std::string foo1(std::string s = std::string("aaa") ) { return s; }

//std::string foo2(std::string s = "aaa" ) { return s; }
//void foo3(A a=A() ) {}


std::string foo(int a) { return "foo1"; }
std::string foo(int a, int b) { return "foo2"; }
std::string foo(int a, int b, int c) { return "foo3: a="+std::to_string(a)+" b="+std::to_string(b)+" c="+std::to_string(c); }




namespace aaaa {

namespace bbbb {

namespace cccc {
std::string foo_cccc(int a) { return "foo1"; }
}

}

namespace dd {
std::string foo_dd(int a) { return "foo1"; }
}

}



// enum class Global_Enum {G_E1, G_E2, G_E3};




// namespace aaaa {

// class A
// {
// public:

// 	A() {}
// 	A(int) {}
// 	A(int, int) {}

// 	int foo() { return 42; }
// 	std::string foo(int a) { return "foo:" + std::to_string(a); }
// 	std::string foo(int a, std::string const &b, int c) { return "foo: " + b + " "+ std::to_string(c); }

// 	void foo_a(A a=A()) {}
// };

// struct S
// {
// 	S() {}

// 	int foo() { return 43; }
// };


// void hi_(void)
// {
// 	std::cout << "Hi in aaaa..." << std::endl;
// }


// };

// void hi(void)
// {
// 	std::cout << "Hi..." << std::endl;
// }
