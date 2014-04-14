#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <string>

class ClassA;
class ClassB;

typedef utility::pointer::access_ptr< ClassA > ClassAAP;
typedef utility::pointer::owning_ptr< ClassA > ClassAOP;
typedef utility::pointer::owning_ptr< ClassA const > ClassACOP;
typedef utility::pointer::owning_ptr< ClassB > ClassBOP;

// Dummy class
class ClassA {
public:
	ClassA() { }
	void owning_ptr_acquire(ClassA *){}
	void add_ref() {}
	void remove_ref() {}
	void show() {}
};

// Fancy class
class ClassB {
private:
	ClassAOP a_;
	ClassBOP b_;
	utility::vector1<ClassAOP> as_;
	std::string s_;
	
public:
	ClassB() { }
	void owning_ptr_acquire(ClassB *){}
	void add_ref() {}
	void remove_ref() {}
	void show() {}
	
	void set_a_null() {
		a_ = NULL;
		a_ = 0;
	}
	
	// Pointer to OP conversion
	void set_a_ptr(ClassA *a) {
		a_ = a;
	}

	void set_a_ptr_casted(ClassA *a) {
		a_ = ClassAOP( a ); // show not match
		// a_.reset( a ); // show not match
	}

	void set_a_op(ClassAOP a) {
		a_ = a;	// should not match
	}

	void set_a_ref(ClassA & a) {
		a_ = &a;
	}
	
	void new_a() {
		a_ = new ClassA;
	}
	
	void new_a_local() {
		ClassAOP aa = new ClassA;
	}
	
	void set_to_self(ClassBOP b) {
		b->b_ = this;
	}

	void test_op(ClassAOP a) {
		runtime_assert( a );
	}

	// Vector1	
	void set_a_vector1(ClassA *a) {
		as_[0] = a;
	}

	void set_aop_vector1(ClassAOP a) {
		as_[0] = a; // should not match
	}

	void set_aref_vector1(ClassA & a) {
		as_[0] = &a;
	}
	
	// operator() call
	ClassA * operator()() {
		return &(*a_); // should not match
	}
	
	// Should not match
	void set_str_char_ptr(const char *s) {
		s_ = s;
	}	
};

void foo() {
	ClassAAP a(new ClassA);
	void *p = a();
}
