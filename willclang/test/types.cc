struct MyStruct {};

typedef MyStruct MSP;

void test_types() {
	int i;
	int *ip;
	int & ir(i);
	MyStruct m;
	MyStruct *mp;
	MyStruct &mr(m);
	MSP msp;
}