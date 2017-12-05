#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	
	//default construct
	histogram<Dynamic> h1;
	//then set dimensions
	h1.setDimensions(2);
	//then set axes
	h1.setAxis(0,new LinearAxis(0,1));
	h1.setAxis(1,new LinearAxis(0,1));
	
	//construct with dimensions
	histogram<Dynamic> h2(1);
	//then set axes
	h2.setAxis(0,new LinearAxis(0,1));
	
	//construct with axes
	histogram<Dynamic> h3(LinearAxis(0,1),LinearAxis(0,1));
	
	std::cout << h1.getDimensions() << std::endl;
	std::cout << h2.getDimensions() << std::endl;
	std::cout << h3.getDimensions() << std::endl;
}
