#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <string>
#include <set>
#include <list>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////
// utility::pointer typedefs

class ClassA;
class ClassB;
class ClassA2;

typedef utility::pointer::access_ptr< ClassA > ClassAAP;
typedef utility::pointer::owning_ptr< ClassA > ClassAOP;
typedef utility::pointer::owning_ptr< ClassA const > ClassACOP;
typedef utility::pointer::owning_ptr< ClassA2 > ClassA2OP;
typedef utility::pointer::owning_ptr< ClassB > ClassBOP;

///////////////////////////////////////////////////////////////////////
// utility::pointer usage 

// Dummy class
class ClassA {
public:
	ClassA() { }
	void owning_ptr_acquire(ClassA *){}
	void add_ref() {}
	void remove_ref() {}
	void show() const {}
	void *operator()() { return this; }
};

class ClassA2 : public ClassA {
public:
	ClassA2() { }
	void method_in_a2() { }
};

// Fancy class
class ClassB {
private:
	ClassAOP a_;
	ClassBOP b_;
	ClassAAP aap_;
	utility::vector1<ClassAOP> as_;
	utility::vector1< utility::pointer::owning_ptr<class ClassA> > as2_;
	std::vector<ClassAOP> as_vector_;
	std::list<ClassAOP> as_list_;
	std::set<ClassAOP> as_set_;
	std::map<std::string,ClassAOP> as_map_;
	std::string s_;
	
public:
	ClassB() : aap_(0) { }
	void owning_ptr_acquire(ClassB *){}
	void add_ref() {}
	void remove_ref() {}
	void show() {}

	///////////////////////////////////////////////////////////////////////
	// Decl

	void bar() {
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > ops_;
		utility::vector1< utility::pointer::access_ptr<class ClassA> > aps_;
	}

	void bar(utility::vector1< utility::pointer::owning_ptr<class ClassA> > & ops_ ) {
		ops_.clear();
	}

	utility::vector1< utility::pointer::owning_ptr<class ClassA> > zzz() {
		// this comment owning_ptr should not get rewritten
		utility::vector1< utility::pointer::owning_ptr<class ClassA> > r;
		return r;
	}	

	static ClassBOP factory() {
		return new ClassB;
	}
	
	///////////////////////////////////////////////////////////////////////
	// shared ptr assignment
  
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
		ClassAOP aa = NULL; // missed
	}

	void set_local_null2() {
		ClassAOP aa(NULL); // missed
	}
	
	void set_local_null3() {
		ClassAOP aa = 0; // missed
	}

	void set_local_null4() {
		ClassAOP aa(0); // missed
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

	// Reference passing
	void set_aref_vector1(ClassA & a) {
		as_[0] = &a;
		as_map_["some"] = &a;
		as_vector_.push_back(&a);
		as_set_.insert(&a);
		as_list_.push_back(&a);
	}
	
	void taking_aref(ClassA const & a) {
		a.show();
	}

	void taking_aop(ClassAOP a) {
		a->show();
	}
	
	void calling_taking_op() {
		ClassAOP a( new ClassA );
		ClassA & aref = *a;
		taking_aref(aref);
		taking_aop(a);
	}
	
	///////////////////////////////////////////////////////////////////////
	// return pointer operator

	// operator() call
	ClassA * operator()() {
		return &(*a_); // should not match
	}
	
	///////////////////////////////////////////////////////////////////////
	// casts
  
#ifdef CASTS
	void dynamic_cast_test(ClassAOP a) {
		ClassA2OP a2 = dynamic_cast< ClassA2 * > ( a() );
		a2->method_in_a2();
	}
#endif

	///////////////////////////////////////////////////////////////////////
	// Calls

	void caller(ClassBOP b) {
		b->new_a();
		b->set_a_null();
		if(b) {
			b->new_a_local();
		}
		ClassB::factory();
		ClassBOP my_b = ClassB::factory();
	}

	///////////////////////////////////////////////////////////////////////
	// Should not match
	
	void set_str_char_ptr(const char *s) {
		s_ = s;
	}

};

///////////////////////////////////////////////////////////////////////
// Functions

void foo() {
	ClassAAP a(new ClassA); // should not match
	ClassAAP b = new ClassA;
	void *p = a(); // for op_void_ptr match tester
}

void bar() {
	utility::vector1< utility::pointer::owning_ptr<class ClassA> > ops_;
	utility::vector1< utility::pointer::access_ptr<class ClassA> > aps_;
}

void bar(utility::vector1< utility::pointer::owning_ptr<class ClassA> > & ops_ ) {
	ops_.clear();
}

utility::vector1< utility::pointer::owning_ptr<class ClassA> > zzz() {
	// this comment owning_ptr should not get rewritten
	utility::vector1< utility::pointer::owning_ptr<class ClassA> > r;
	return r;
}
