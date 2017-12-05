#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	
	histogram<1> h1(FixedLinearAxis(1,5,4));
	h1.add(0,amount(2));
	h1.add(6,amount(3));
	
	std::cout << h1.getUnderflow() << ' ' << h1.getOverflow() << std::endl;
	
	histogram<Dynamic> h2;
	h2.setDimensions(1);
	h2.setAxis(0,new FixedLinearAxis(1,5,4));
	h2.add(0,amount(4));
	h2.add(6,amount(5));
	
	std::cout << h2.getUnderflow() << ' ' << h2.getOverflow() << std::endl;
}