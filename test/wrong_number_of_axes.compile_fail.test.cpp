#include <PhysTools/histogram.h>

int main(){
	using namespace phys_tools::histograms;
#pragma phys_tools_testing failure_begin
	histogram<2> h1(LinearAxis(0,1));
#pragma phys_tools_testing failure_end
}