#include <test_ref.fwd.hh>

struct MyStruct {
	int i;
};

int msfunc(MyStruct & ms) {
	return ms.i;
}