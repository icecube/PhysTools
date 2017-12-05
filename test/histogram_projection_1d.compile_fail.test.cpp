#include <PhysTools/histogram.h>

int main(){
	using namespace phys_tools::histograms;
	histogram<1> h(LinearAxis(0,1));
	h.add(0,amount(1));
	h.add(1,amount(2));
#pragma phys_tools_testing failure_begin
	auto proj=h.project(0);
#pragma phys_tools_testing failure_end
}