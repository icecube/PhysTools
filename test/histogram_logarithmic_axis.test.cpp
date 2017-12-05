#include <PhysTools/histogram.h>
#include <iostream>

int main(){
	using namespace phys_tools::histograms;
	std::cout.precision(7);
	
	{ //no limits
		histogram<1> h_base(LogarithmicAxis(0,0.5));
		h_base.setUseContentScaling(false);
		
		{ //bin edges should count toward the bin to the left
			histogram<1> h(h_base);
			h.add(10);
			std::cout << h; //10 1
		}
		
		{ //middles of bins should be obvious
			histogram<1> h(h_base);
			h.add(20);
			std::cout << h; //10 1
		}
		
		{ //bin edges should count toward the bin to the left
			histogram<1> h(h_base);
			h.add(100);
			std::cout << h; //100 1
		}
		
		{ //middles of bins should be obvious
			histogram<1> h(h_base);
			h.add(200);
			std::cout << h; //100 1
		}
		
		{ //adding to a bin to the left, then one to the right should work
			histogram<1> h(h_base);
			h.add(20);
			h.add(200);
			std::cout << h; //4 1\n5 0\n6 1
		}
		
		{ //adding to a bin to the right, then one to the left should also work
			histogram<1> h(h_base);
			h.add(200);
			h.add(20);
			std::cout << h; //4 1\n5 0\n6 1
		}
		
		{ //zeros and infinites should under/overflow
			histogram<1> h(h_base);
			h.add(0);
			h.add(std::numeric_limits<double>::infinity());
			std::cout << "Underflow: " << h.getUnderflow() << std::endl;
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			h.add(-1); //negative number should underflow
			h.add(std::numeric_limits<double>::quiet_NaN()); //NaNs should overflow
			std::cout << "Underflow: " << h.getUnderflow() << std::endl;
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			std::cout << std::endl;
		}
	}
	
	{ //limits
		histogram<1> h_base(LogarithmicAxis(0,0.5));
		h_base.setUseContentScaling(false);
		h_base.getAxis(0)->setLowerLimit(20);
		h_base.getAxis(0)->setUpperLimit(2000);
		
		{ //allowed values should be inserted normally
			histogram<1> h(h_base);
			h.add(100);
			std::cout << h; //100 1
		}
		
		{ //values below the lower limit should underflow
			histogram<1> h(h_base);
			h.add(10);
			std::cout << "Underflow: " << h.getUnderflow() << std::endl;
			std::cout << h; //prints nothing!
			std::cout << std::endl;
		}
		
		{ //the lower limit itself should _not_ underflow
			histogram<1> h(h_base);
			h.add(20);
			std::cout << "Underflow: " << h.getUnderflow() << std::endl;
			std::cout << h;
		}
		
		{ //values above the upper limit should overflow
			histogram<1> h(h_base);
			h.add(2500);
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			std::cout << h; //prints nothing!
			std::cout << std::endl;
		}
		
		{ //the upper limit itself _should_ underflow
			histogram<1> h(h_base);
			h.add(2000);
			std::cout << "Overflow: " << h.getOverflow() << std::endl;
			std::cout << h; //prints nothing!
			std::cout << std::endl;
		}
	}
}