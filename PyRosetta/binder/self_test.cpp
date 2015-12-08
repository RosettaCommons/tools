
#include <vector>

enum class Global_Enum {G_E1, G_E2, G_E3};



namespace aaaa {

class A
{
	A() {}
	A(int) {}
	A(int, int) {}
};

};


namespace bbbb {

class B : public aaaa::A
{

public:
	class B_A {};

	enum class B_Enum {B_E1, B_E2, B_E3} ;

	void foo() {}

	B_Enum e;

	int a, b, c;

	static void B_static_foo() {}
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
void global_function_foo( bbbb::Q<int> )
{
	struct Inner_struct {};
	enum class Inner_Enum {I_E1, I_E2, I_E3};

	int fpppp(int a);

}

using Vector_A = std::vector<aaaa::A>;

typedef std::vector<aaaa::A> _Vector_A_typedef;

void global_function_foo_vector( std::vector<aaaa::A> ) {}


typedef bbbb::Q<int> T_Q_int;
using U_Q_int = bbbb::Q<int>;
using U_Q_float = bbbb::Q<float>;


void global_function_foo_t( T_Q_int );
void global_function_foo_u( U_Q_int );
void global_function_foo_u( U_Q_float _) { _.foo(); _.value += 1; }


class T_W : public T_Q_int {};
class U_W : public U_Q_int {};
//class U_W_float : public U_Q_float {};
