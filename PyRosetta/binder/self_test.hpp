#include <iostream>

enum class Global_Enum {G_E1, G_E2, G_E3};



namespace aaaa {

class A
{
	A() {}
public:
	A(int) {}
protected:
	A(int, int) {}
};

struct S
{
	S() {}
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
