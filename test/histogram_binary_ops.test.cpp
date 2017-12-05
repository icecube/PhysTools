#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	
	histogram<2> h1(LinearAxis(0,1),LinearAxis(0,1));
	histogram<2> h2(LinearAxis(0,1),LinearAxis(0,1));
	
	h1.add(1,1,amount(2));
	h1.add(2,2,amount(3));
	
	h2.add(0,0,amount(2));
	h2.add(0,1,amount(2));
	h2.add(1,0,amount(2));
	h2.add(1,1,amount(2));
	
	std::cout << h1+h2 << std::endl;
	std::cout << h1-h2 << std::endl;
	std::cout << h1*h2 << std::endl;
	std::cout << h1/h2 << std::endl;
}