#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	
	LinearAxis a(0,5);
	a.setUpperLimit(26);
	
	histogram<1> h(a);
	
	for(double x=0; x<30; x++)
		h.add(x);
	
	std::cout << h << std::endl;
	std::cout << "overflow " << h.getOverflow() << std::endl;
	
	for(auto it=h.begin(), end=h.end(); it!=end; it++){
		auto other=h.findBinIterator(it);
		if(other==end)
			std::cout << "OUT OF RANGE" << std::endl;
		else
			std::cout << other.getCoordinate(0) << std::endl;
	}
}