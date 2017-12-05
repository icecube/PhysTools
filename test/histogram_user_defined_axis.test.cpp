#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	
	histogram<1> h_base(FixedUserDefinedAxis({0,1.1,2.5,6.7,8.0}));
	h_base.setUseContentScaling(false);
	
	{ //lower edge of the first bin should land in the first bin
		histogram<1> h(h_base);
		h.add(0);
		std::cout << h; //0 1
    }
	
	{ //values inside first bin should land in the first bin
		histogram<1> h(h_base);
		h.add(1);
		std::cout << h; //0 1
    }
	
	{ //lower edge of the second bin should land in the second bin
		histogram<1> h(h_base);
		h.add(1.1);
		std::cout << h; //1.1 1
    }
	
	{ //values inside second bin should land in the second bin
		histogram<1> h(h_base);
		h.add(2);
		std::cout << h; //1.1 1
    }
	
	{ //lower edge of the last bin should land in the last bin
		histogram<1> h(h_base);
		h.add(6.7);
		std::cout << h; //6.7 1
    }
	
	{ //values inside last bin should land in the last bin
		histogram<1> h(h_base);
		h.add(7);
		std::cout << h; //6.7 1
    }
	
	{ //adding to a bin to the left, then one to the right should work
		histogram<1> h(h_base);
		h.add(2);
		h.add(4);
		std::cout << h; //1.1 1\n2.5 1
	}
	
	{ //adding to a bin to the right, then one to the left should also work
		histogram<1> h(h_base);
		h.add(4);
		h.add(2);
		std::cout << h; //1.1 1\n2.5 1
	}
}