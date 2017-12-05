#include <iostream>
#include <PhysTools/histogram.h>

int main(){
	using namespace phys_tools::histograms;

	histogram<2> h1(LinearAxis(0,1),LinearAxis(1,5));
	histogram<2>::internalCoordinate intC[2];
	double extC[2];

	extC[0]=.5;
	extC[1]=2;
	h1.getCoordinates(extC,intC);
	std::cout << '(' << extC[0] << ',' << extC[1] << ") -> (" 
		<< intC[0] << ',' << intC[1] << ')' << std::endl;
}
