#include "PhysTools/histogram.h"
#include "PhysTools/bin_types.h"

int main(){
	using namespace phys_tools::histograms;
	
	histogram<1,sqErrorValue<>> h(LinearAxis(0,1));
	//implicit errors
	h.add(0,amount(1));
	h.add(1,amount(1));
	h.add(2,amount(3));
	h.add(3);
	h.add(3);
	h.add(4,amount(1));
	h.add(4,amount(1));
	h.add(5,amount(3));
	h.add(5,amount(3));
	//explicit errors
	h.add(6,amount(sqErrorValue<>(2,1)));
	h.add(7,amount(sqErrorValue<>(2,1)));
	h.add(7,amount(sqErrorValue<>(2,1)));
	
	std::cout.precision(3);
	std::cout << h;
}