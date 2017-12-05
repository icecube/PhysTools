#include <PhysTools/histogram.h>

int main(){
	using namespace phys_tools::histograms;
	
	//default construct
	histogram<2> h1;
	//set axes after the fact
	h1.setAxes(LinearAxis(0,1),LinearAxis(0,1));
	
	//construct with axes
	histogram<2> h2(LinearAxis(0,1),LinearAxis(0,1));
}