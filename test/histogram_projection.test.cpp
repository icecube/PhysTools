#include <PhysTools/histogram.h>
#include <iostream>

template<typename HistType>
void fillHist(HistType& h){
	using phys_tools::histograms::amount;
	
	h.add(0,0,0,amount(1));
	h.add(1,0,0,amount(2));
	h.add(0,1,0,amount(3));
	h.add(1,1,0,amount(4));
	h.add(0,0,1,amount(5));
	h.add(1,0,1,amount(6));
	h.add(0,1,1,amount(7));
	h.add(1,1,1,amount(8));
}

template<typename HistType>
void testProjection(const HistType& h){
	for(unsigned int i=0; i<3; i++)
		std::cout << h.project(i) << std::endl;
}

int main(){
	using namespace phys_tools::histograms;
	
	histogram<3> h1(LinearAxis(0,1),LinearAxis(0,1),LinearAxis(0,1));
	fillHist(h1);
	testProjection(h1);
	
	histogram<Dynamic> h2(LinearAxis(0,1),LinearAxis(0,1),LinearAxis(0,1));
	fillHist(h2);
	testProjection(h2);
}
