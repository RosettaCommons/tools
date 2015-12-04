

namespace aaaa {

class A
{
};

};


namespace bbbb {

class B : public aaaa::A
{

public:
	class B_A {};

	void foo() {}
};

template <class T>
class Q : public B
{
public:
	T value;
};

};



struct S
{
	int s1;
};


void global_function_foo() {}

void global_function_foo( bbbb::Q<int> ) {}


typedef bbbb::Q<int> T_Q_int;
using U_Q_int = bbbb::Q<int>;
using U_Q_float = bbbb::Q<float>;


void global_function_foo_t( T_Q_int );
void global_function_foo_u( U_Q_int );
void global_function_foo_u( U_Q_float _) { _.foo(); _.value += 1; }


class T_W : public T_Q_int {};
class U_W : public U_Q_int {};
//class U_W_float : public U_Q_float {};
