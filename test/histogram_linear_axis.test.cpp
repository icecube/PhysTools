#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	
	{ //no limits
		histogram<1> h_base(LinearAxis(0,1));
		
		{ //bin edges should count toward the bin to the left
			histogram<1> h(h_base);
			h.add(0);
			std::cout << h; //0 1
		}
		
		{ //middles of bins should be obvious
			histogram<1> h(h_base);
			h.add(0.5);
			std::cout << h; //0 1
		}
		
		{ //bin edges should count toward the bin to the left
			histogram<1> h(h_base);
			h.add(17);
			std::cout << h; //17 1
		}
		
		{ //middles of bins should be obvious
			histogram<1> h(h_base);
			h.add(17.5);
			std::cout << h; //17 1
		}
		
		{ //adding to a bin to the left, then one to the right should work
			histogram<1> h(h_base);
			h.add(4);
			h.add(6);
			std::cout << h; //4 1\n5 0\n6 1
		}
		
		{ //adding to a bin to the right, then one to the left should also work
			histogram<1> h(h_base);
			h.add(6);
			h.add(4);
			std::cout << h; //4 1\n5 0\n6 1
		}
		
		{ //infinites should under/overflow
			histogram<1> h(h_base);
			h.add(-std::numeric_limits<double>::infinity());
			h.add(std::numeric_limits<double>::infinity());
			std::cout << "Underflow: " << h.getUnderflow() << std::endl;
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			h.add(std::numeric_limits<double>::quiet_NaN()); //NaNs should overflow
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			std::cout << std::endl;
		}
	}
	
	{ //limits
		histogram<1> h_base(LinearAxis(0,1));
		h_base.getAxis(0)->setLowerLimit(-7);
		h_base.getAxis(0)->setUpperLimit(22);
		
		{ //allowed values should be inserted normally
			histogram<1> h(h_base);
			h.add(-2);
			std::cout << h; //-2 1
		}
		
		{ //values below the lower limit should underflow
			histogram<1> h(h_base);
			h.add(-8);
			std::cout << "Underflow: " << h.getUnderflow() << std::endl;
			std::cout << h; //prints nothing!
			std::cout << std::endl;
		}
		
		{ //the lower limit itself should _not_ underflow
			histogram<1> h(h_base);
			h.add(-7);
			std::cout << "Underflow: " << h.getUnderflow() << std::endl;
			std::cout << h;
		}
		
		{ //values above the upper limit should overflow
			histogram<1> h(h_base);
			h.add(25);
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			std::cout << h; //prints nothing!
			std::cout << std::endl;
		}
		
		{ //the upper limit itself _should_ underflow
			histogram<1> h(h_base);
			h.add(22);
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			std::cout << h; //prints nothing!
			std::cout << std::endl;
		}
	}
}