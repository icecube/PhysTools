#include <iostream>
#include <PhysTools/histogram.h>

using namespace phys_tools::histograms;

int main(){
	histogram<2> h(LinearAxis(0,1),LinearAxis(0,1));
	
	//this has too few arguments and should not compile
#pragma phys_tools_testing failure_begin
	h.add(2);
#pragma phys_tools_testing failure_end
	std::cout << h;
}