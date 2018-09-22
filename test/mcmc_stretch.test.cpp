#include <iostream>
#include <cmath>
#include <random>
#include <cassert>

#include <PhysTools/mcmc.h>
#include <PhysTools/histogram.h>
#include <PhysTools/optimization/ParameterSet.h>

const double mu=2;
const double sigma=.7;

struct llh{
	double operator()(const std::vector<double>& c) const{
		assert(c.size()==1);
		double z=(c[0]-mu)/sigma;
		return(.5*z*z);
	}
};

int main(){
	using namespace phys_tools;
	std::mt19937 rng(137);
	std::vector<std::vector<double>> initialEnsemble;
	for(size_t i=0; i<10; i++)
		initialEnsemble.push_back(std::vector<double>{.1*i});

	ParameterSet params;
	params.addParameter("x");
	params.setParameterLowerLimit("x",-1);
	params.setParameterUpperLimit("x",5);
	const size_t nSamples=5000000;
	auto samples=markovSample(llh(),phys_tools::StretchMove(),params,
					 		 nSamples,10000,10,rng,std::move(initialEnsemble));
	phys_tools::histograms::histogram<1> h(phys_tools::histograms::LinearAxis(0,.2));
	h.setUseContentScaling(false);
	for(const auto& sample : samples)
		h.add(sample.coordinates[0]);
	//std::cout << h << std::endl;
	double expSum=0, obsSum=0;
	double chi2=0;
	size_t nBins=0;
	for(auto it=h.begin(); it!=h.end(); it++){
		double tMin=(it.getBinEdge(0)-mu)/(sqrt(2)*sigma);
		double tMax=(it.getBinEdge(0)+it.getBinWidth(0)-mu)/(sqrt(2)*sigma);
		double expected=nSamples*(erf(tMax)-erf(tMin))/2;
		double observed=*it;
		//std::cout << "Obs: " << observed << " Exp: " << expected <<
		//" [" << expected-sqrt(expected) << ',' << expected+sqrt(expected) << ']';
		//if(observed>expected-sqrt(expected) && observed<expected+sqrt(expected))
		//	std::cout << " *";
		//std::cout << std::endl;
		expSum+=expected;
		double diff=(observed-expected);
		chi2+=diff*diff/expected;
		nBins++;
	}
	//std::cout << "ExpSum: " << expSum << std::endl;
	//std::cout << "Chi2: " << chi2 << std::endl;
	//std::cout << "Chi2/ndof: " << chi2/nBins << std::endl;
	if(std::abs(chi2/nBins-1)>.1)
		std::cout << "Distribution not well sampled" << std::endl;
}
