#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	
	histogram<2> h1(LinearAxis(0,1),LinearAxis(0,1));
	h1.add(0,0,amount(1));
	h1.add(0,1,amount(2));
	h1.add(1,0,amount(3));
	h1.add(1,1,amount(4));
	std::cout << h1 << std::endl;
	
	histogram<Dynamic> h2(LinearAxis(0,1),LinearAxis(0,1));
	h2.add(0,0,amount(1));
	h2.add(0,1,amount(2));
	h2.add(1,0,amount(3));
	h2.add(1,1,amount(4));
	std::cout << h2 << std::endl;
}