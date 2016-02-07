#include <string>
//#include <iostream>

/*
namespace aa {

class A {
public:
	int n;
	float f;
	double d;
};

typedef int aa_INT;

inline void foo_aa(A& a, aa_INT) {}

typedef void (* UtilityExitCallBack)(void);

void set_main_exit_callback( UtilityExitCallBack = 0 );

}

//std::string foo1(std::string s = std::string("aaa") ) { return s; }

//std::string foo2(std::string s = "aaa" ) { return s; }
//void foo3(A a=A() ) {}

inline bool foo(bool a, bool&b, aa::A) { return true; }
// std::string foo(int a) { return "foo1"; }
// std::string foo(int a, int b) { return "foo2"; }
// std::string foo(int a, int b, int c) { return "foo3: a="+std::to_string(a)+" b="+std::to_string(b)+" c="+std::to_string(c); }


inline void foo_A(aa::A& a, aa::aa_INT) {}

*/


namespace utility {

// class UA;

// class A {
// public:
// 	int n;

// 	void foo() {};
// };

template <typename T>
class TUA;

template <typename T>
class TUA {
public:
	T t;
};


typedef TUA<int> TUA_int;


void foo_TUA(TUA_int a) {};


template< bool >
struct vectorL_IndexSelector
{
};


/// @brief vectorL index type selector: Negative lower index specialization
template<>
struct vectorL_IndexSelector< false >
{
	inline
	static
	bool
	ge() { return true; }
};


class CSI_Sequence
{
public:
	/// @brief constructor
	CSI_Sequence(std::string sequence) {}

	operator std::string() const { return sequence_; }

	/// @brief operator to output our sequence so we can write: std::cout << CSI_SequenceObject
	//friend std::ostream & operator << (std::ostream & os, CSI_Sequence const &sq) { os << sq.sequence_; return os; }

	void foo(int) const {};
private:
	std::string sequence_;
};

template <class T>
void swap(TUA<T>) {}

void foo_i(int **a) {}
std::string foo_s(std::string const &a="aaa") { return a; }



} // namespace utility
