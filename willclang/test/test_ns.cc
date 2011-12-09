#include <include/test_ns.hh>
int func() {
	ns1::val1;
	ns1::ns2::val2;
	ns1::ns2::ns3::val3;
	ns1::ns2::ns3::ns4::val4;
	ns1::ns2::ns3::ns4::ns5::val5;

	using namespace ns1;
	using namespace ns2;
	using namespace ns3;

	val1;
	val2;
	val3;
	ns4::val4;
	ns4::ns5::val5;

	using namespace ns4;
	using namespace ns5;

	val1;
	val2;
	val3;
	val4;
	val5;
}