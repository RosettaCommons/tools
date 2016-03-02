#include <string>
#include <vector>
#include <memory>
//#include <iostream>

#include <self_test.incl.hpp>

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

enum E1 {E1A, E1B};


namespace utility {

// template <class T> class N
// {
// 	enum OptionTypes {
// 		UNKNOWN_OPTION,
// 	};

// 	void private_foo(std::vector<OptionTypes> );

// public:
// 	void foo() {}
// };

// void foo_N(N<int>) {}


//void foo() {}
//void foo_string_const(std::string) {}

// void foo_string_const(std::string const &) {}
// void foo_string(std::string  &) {}
// void foo_int_const(int const &) {}
// void foo_int(int &) {}

// class CA {};
// template<typename T>
// void foo_aaaa(int) {}

// //template<>
// //void foo_aaaa<aaaa::A>(int);

// void foo()
// {
// 	foo_aaaa<aaaa::A>(0);
// }



template<typename T>
class TA
{
public:
	T value;
};


template<int S>
class TI
{
public:
	int values[S];
};


class B : public TA<aaaa::A> {};

//void foo_A( std::shared_ptr< const aaaa::A > &) {}
//void foo_A( std::shared_ptr< TA<aaaa::A> > &) {}

typedef TA<aaaa::A> TA_A;

void foo_TA(TA<aaaa::A>) {}

template<int T>
void foo_T(TI<T> &) {}

// typedef TA<aaaa::A> TA_A;

// class TC : public TA<aaaa::A> {};

// class A {
// public:
// 	int n;

// 	void foo() {};
// };

// template <typename T>
// class TUA;

// template <typename T>
// class TUA {
// public:
// 	T t;
// };


// typedef TUA<int> TUA_int;


// void foo_TUA(TUA_int a) {};


// template< bool >
// struct vectorL_IndexSelector
// {
// };


// /// @brief vectorL index type selector: Negative lower index specialization
// template<>
// struct vectorL_IndexSelector< false >
// {
// 	inline
// 	static
// 	bool
// 	ge() { return true; }
// };


// template <class T>
// void swap(TUA<T>) {}

// void foo_i(int **a) {}
// std::string foo_s(std::string const &a="aaa") { return a; }


// class TUAC : public TUA_int
// {};

// void foo(TUA_int) {}
// void foo(TUA<float>) {}

// namespace csi {

// enum E2 {E2A, E2B};


// namespace inner1 {
// namespace inner2 {
// namespace inner3 {
// void foo() {}
// }
// }
// }

// class CSI_Sequence
// {
// public:
// 	/// @brief constructor
// 	CSI_Sequence(std::string sequence) {}

// 	operator std::string() const { return sequence_; }

// 	/// @brief operator to output our sequence so we can write: std::cout << CSI_SequenceObject
// 	//friend std::ostream & operator << (std::ostream & os, CSI_Sequence const &sq) { os << sq.sequence_; return os; }

// 	void foo(int) const {};
// private:
// 	std::string sequence_;
// };

// void foo_csi(char * const argv[]) {}

// class ESFT : public std::enable_shared_from_this<ESFT>
// {
// public:
// 	int a;
// };


// class ESFT2 : public ESFT
// {
// public:
// 	std::shared_ptr<ESFT2> shared_from_this() {
//         return shared_from_this();
//     }

// 	void * fn;
// 	int a;

// };

// void foo(std::string a = "abc") {}
// std:: string foo(int a) { return ""; }

//} // namespace csi

} // namespace utility
