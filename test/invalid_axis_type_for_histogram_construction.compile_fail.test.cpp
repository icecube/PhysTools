#include <PhysTools/histogram.h>

int main(){
	using namespace phys_tools::histograms;
#pragma phys_tools_testing failure_begin
	histogram<1>(5); //int is not a valid axis type!
#pragma phys_tools_testing failure_end
}
