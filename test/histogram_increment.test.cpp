#include <iostream>
#include <PhysTools/histogram.h>

using namespace phys_tools::histograms;

void both_empty(histogram<1> h1, histogram<1> h2){
	std::cout << "Both empty:" << std::endl;
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void first_empty(histogram<1> h1, histogram<1> h2){
	std::cout << "First empty:" << std::endl;
	h2.add(1,amount(2));
	h2.add(4,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void second_empty(histogram<1> h1, histogram<1> h2){
	std::cout << "Second empty:" << std::endl;
	h1.add(1,amount(1));
	h1.add(2,amount(2));
	h1.add(3,amount(3));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_nn(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: neither, neither):" << std::endl;
	h1.add(1,amount(2));
	h1.add(3,amount(2));
	h2.add(1,amount(1));
	h2.add(2,amount(2));
	h2.add(3,amount(3));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_n1(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: neither, 1):" << std::endl;
	h1.add(1,amount(2));
	h1.add(4,amount(2));
	h2.add(1,amount(1));
	h2.add(2,amount(2));
	h2.add(3,amount(3));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_n2(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: neither, 2):" << std::endl;
	h1.add(1,amount(1));
	h1.add(2,amount(2));
	h1.add(3,amount(3));
	h2.add(1,amount(2));
	h2.add(4,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_1n(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: 1, neither):" << std::endl;
	h1.add(0,amount(2));
	h1.add(3,amount(2));
	h2.add(1,amount(1));
	h2.add(2,amount(2));
	h2.add(3,amount(3));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_11(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: 1, 1):" << std::endl;
	h1.add(0,amount(2));
	h1.add(3,amount(2));
	h2.add(1,amount(1));
	h2.add(2,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_12(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: 1, 2):" << std::endl;
	h1.add(0,amount(2));
	h1.add(2,amount(2));
	h2.add(1,amount(1));
	h2.add(3,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_2n(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: 2, neither):" << std::endl;
	h1.add(1,amount(1));
	h1.add(2,amount(2));
	h1.add(3,amount(3));
	h2.add(0,amount(2));
	h2.add(3,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_21(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: 2, 1):" << std::endl;
	h1.add(1,amount(1));
	h1.add(2,amount(2));
	h1.add(3,amount(3));
	h2.add(0,amount(2));
	h2.add(2,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_22(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: 2, 2):" << std::endl;
	h1.add(1,amount(1));
	h1.add(2,amount(2));
	h2.add(0,amount(2));
	h2.add(3,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_nonempty_no_overlap(histogram<1> h1, histogram<1> h2){
	std::cout << "Both nonempty (larger extent: 1, 2; no overlap):" << std::endl;
	h1.add(1,amount(1));
	h1.add(2,amount(2));
	h2.add(3,amount(2));
	h2.add(4,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void overflow_and_underflow(histogram<1> h1, histogram<1> h2){
	std::cout << "Both have underflow and overflow:" << std::endl;
	h1.getAxis(0)->setLowerLimit(0);
	h1.getAxis(0)->setUpperLimit(5);
	h2.getAxis(0)->setLowerLimit(0);
	h2.getAxis(0)->setUpperLimit(5);
	h1.add(-1,amount(1));
	h1.add(2,amount(1));
	h1.add(6,amount(1));
	h2.add(-3,amount(2));
	h2.add(3,amount(3));
	h2.add(8,amount(4));
	auto r = h1+h2;
	std::cout << r;
	std::cout << "U " << r.getUnderflow() << std::endl;
	std::cout << "O " << r.getOverflow() << std::endl;
	std::cout << std::endl;
}

void first_uninitialized(histogram<1> h2){
	std::cout << "First uninitialized:" << std::endl;
	histogram<1> h1; //not initialized with axes!
	h2.add(3,amount(1));
	h2.add(4,amount(2));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void second_uninitialized(histogram<1> h1){
	std::cout << "Second uninitialized:" << std::endl;
	histogram<1> h2; //not initialized with axes!
	h1.add(1,amount(2));
	h1.add(3,amount(4));
	auto r = h1+h2;
	std::cout << r << std::endl;
}

void both_uninitialized(){
	std::cout << "Both uninitialized:" << std::endl;
	histogram<1> h1,h2; //not initialized with axes!
	auto r = h1+h2;
	std::cout << r << std::endl;
}

int main(){
	histogram<1> h1(LinearAxis(0,1)), h2(LinearAxis(0,1));
	
	both_empty(h1,h2);
	first_empty(h1,h2);
	second_empty(h1,h2);
	both_nonempty_nn(h1,h2);
	both_nonempty_n1(h1,h2);
	both_nonempty_n2(h1,h2);
	both_nonempty_1n(h1,h2);
	both_nonempty_11(h1,h2);
	both_nonempty_12(h1,h2);
	both_nonempty_2n(h1,h2);
	both_nonempty_21(h1,h2);
	both_nonempty_22(h1,h2);
	both_nonempty_no_overlap(h1,h2);
	overflow_and_underflow(h1,h2);
	first_uninitialized(h2);
	second_uninitialized(h1);
	both_uninitialized();
}