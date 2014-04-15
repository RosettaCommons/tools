#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <string>
#include <set>
#include <list>
#include <vector>
#include <map>

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
	void *operator()() { return this; }
};

// Fancy class
class ClassB {
private:
	ClassAOP a_;
	ClassBOP b_;
	utility::vector1<ClassAOP> as_;
	utility::vector1< utility::pointer::owning_ptr<class ClassA> > as2_;
	std::vector<ClassAOP> as_vector_;
	std::list<ClassAOP> as_list_;
	std::set<ClassAOP> as_set_;
	std::map<std::string,ClassAOP> as_map_;
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
		a_ = ClassAOP( a ); // should not match
		// a_.reset( a ); // should not match
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

	void set_local_null1() {
		ClassAOP aa = NULL;
	}

	void set_local_null2() {
		ClassAOP aa(NULL);
	}
	
	void set_local_null3() {
		ClassAOP aa = 0;
	}

	void set_local_null4() {
		ClassAOP aa(0);
	}

	ClassAOP new_aap() {
		return new ClassA;
	}
	
	void set_to_self(ClassBOP b) {
		b->b_ = this; // should not match
	}

	void test_op(ClassAOP a) {
		runtime_assert( a );
	}

	// Vector1	
	void set_a_vector1(ClassA *a) {
		as_[0] = a;
		as_map_["some"] = a;
		as_vector_.push_back(a);
		as_set_.insert(a);
		as_list_.push_back(a);
	}

	void set_aop_vector1(ClassAOP a) {
		// should not match
		as_[0] = a;
		as_map_["some"] = a;
		as_vector_.push_back(a);
		as_set_.insert(a);
		as_list_.push_back(a);
	}

	void set_aref_vector1(ClassA & a) {
		as_[0] = &a;
		as_map_["some"] = &a;
		as_vector_.push_back(&a);
		as_set_.insert(&a);
		as_list_.push_back(&a);
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
	ClassAAP a(new ClassA); // should not match
	ClassAAP b = new ClassA;
	void *p = a(); // for op_void_ptr match tester
}
