#include <iostream>

namespace foo {
	namespace bar {
		int func(int i) {
			return 2*i;
		}
	}
}

int main (int argc, char const *argv[])
{
	std::cout << 30 << " " << foo::bar::not_func(30) << std::endl;
	return 0;
}