#include <PhysTools/histogram.h>
#include <iostream>

template<typename HistType>
void exerciseIterators(HistType& h){
	assert(h.getDimensions()==2);
	
	auto it=h.begin(), begin=h.begin(), end=h.end();
	//we can iterate over the histogram, moving most rapdily through the first dimension
	for(; it!=end; it++)
		std::cout << (double)*it << ' ' << it.getCoordinate(0) << ' ' << it.getCoordinate(1) << ' '
			<< it.getBinEdge(0) << ' ' << it.getBinEdge(1) << std::endl;
	//having reached the end of the histogram, we can go back again
	while(it!=begin){
		it--;
		std::cout << (double)*it << ' ' << it.getCoordinate(0) << ' ' << it.getCoordinate(1) << ' '
			<< it.getBinEdge(0) << ' ' << it.getBinEdge(1) << std::endl;
	}
	
	auto rit=h.rbegin(), rbegin=h.rbegin(), rend=h.rend();
	//we can also do the whole exercise in reverse
	for(; rit!=rend; rit++)
		std::cout << (double)*rit << ' ' << rit.getCoordinate(0) << ' ' << rit.getCoordinate(1) << ' '
			<< rit.getBinEdge(0) << ' ' << rit.getBinEdge(1) << std::endl;
	//having reached the end of the histogram, we can go back again
	while(rit!=rbegin){
		rit--;
		std::cout << (double)*rit << ' ' << rit.getCoordinate(0) << ' ' << rit.getCoordinate(1) << ' '
			<< rit.getBinEdge(0) << ' ' << rit.getBinEdge(1) << std::endl;
	}
}

int main(){
	using namespace phys_tools::histograms;
	
	histogram<2> h1(LinearAxis(0,1),LinearAxis(0,1));
	h1.add(1,1,amount(1));
	h1.add(1,2,amount(2));
	h1.add(2,2,amount(3));
	h1.add(2,1,amount(4));
	exerciseIterators(h1);
	
	histogram<Dynamic> h2(LinearAxis(0,1),LinearAxis(0,1));
	h2.add(1,1,amount(1));
	h2.add(1,2,amount(2));
	h2.add(2,2,amount(3));
	h2.add(2,1,amount(4));
	exerciseIterators(h2);
	
	//attempting to iterate over empty histograms should safely do nothing
	histogram<2> empty1(LinearAxis(0,1),LinearAxis(0,1));
	exerciseIterators(empty1);
	histogram<Dynamic> empty2(LinearAxis(0,1),LinearAxis(0,1));
	exerciseIterators(empty2);
}