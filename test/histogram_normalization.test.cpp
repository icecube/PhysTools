#include <PhysTools/histogram.h>

template<typename HistType>
void renormalize_histogram(HistType& h){
	double n=h.integral();
	h/=n;
	assert(std::abs(h.integral()-1.0)<1e-6);
	h*=5;
	assert(std::abs(h.integral()-5.0)<1e-6);
}

int main(){
	using namespace phys_tools::histograms;
	
	histogram<2> h1(LinearAxis(0,1),LinearAxis(0,1));
	h1.add(0,1,amount(5));
	h1.add(1,2,amount(7));
	h1.add(2,0,amount(3));
	h1.add(2,1,amount(4));
	renormalize_histogram(h1);
	
	histogram<Dynamic> h2(LinearAxis(0,1),LinearAxis(0,1));
	h2.add(0,1,amount(5));
	h2.add(1,2,amount(7));
	h2.add(2,0,amount(3));
	h2.add(2,1,amount(4));
	renormalize_histogram(h2);
	
	histogram<2> h3(LogarithmicAxis(0,.5),LogarithmicAxis(0,.5));
	h3.add(5,10,amount(5));
	h3.add(10,50,amount(7));
	h3.add(50,5,amount(3));
	h3.add(50,10,amount(4));
	renormalize_histogram(h3);
	
	histogram<Dynamic> h4(LogarithmicAxis(0,.5),LogarithmicAxis(0,.5));
	h4.add(5,10,amount(5));
	h4.add(10,50,amount(7));
	h4.add(50,5,amount(3));
	h4.add(50,10,amount(4));
	renormalize_histogram(h4);
}