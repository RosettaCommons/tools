#include <iostream>

enum class Global_Enum {G_E1, G_E2, G_E3};



namespace aaaa {

class A
{
	A() {}
	A(int) {}
	A(int, int) {}
};

};

void hi(void)
{
	std::cout << "Hi..." << std::endl;
}
