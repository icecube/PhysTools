#include <iostream>
#include <PhysTools/likelihood/likelihood.h>

using namespace phys_tools;
using namespace phys_tools::likelihood;

//not a 'problem' at all, but an analytically defined function for which we know
//all true properties
struct problem{
	static constexpr double xm=3;
	static constexpr double ym=2;
	static constexpr double zm=1;
	static constexpr double xs=1;
	static constexpr double ys=2;
	static constexpr double zs=3;

	double llh(const std::vector<double>& params) const{
		double x=params[0];
		double y=params[1];
		double z=params[2];

		double tx=.5*(x-xm)*(x-xm)/(xs*xs);
		double ty=.5*(y-ym)*(y-ym)/(ys*ys);
		double tz=.5*(z-zm)*(z-zm)/(zs*zs);
		return(tx+ty+tz);
	}
};

likelihoodPoint fit(const problem& prob, const ParameterSet& ps){
	likelihoodPoint result;
	result.params={problem::xm,problem::ym,problem::zm};
	for(unsigned int i=0; i<3; i++){
		if(ps.isFixed(i))
			result.params[i]=ps.getParameterValue(i);
	}
	result.likelihood=prob.llh(result.params);
	return(result);
}

int main(){
	const double tol=.01;
	problem prob;
	ParameterSet ps;
	ps.addParameter("x");
	ps.setParameterProperty("x","UncertaintySearchPlan",
							UncertaintySearchPlan{[](double x){ return(x-1); },
								[](double x){ return(x+1); }, tol});
	ps.addParameter("y");
	ps.setParameterProperty("y","UncertaintySearchPlan",
							UncertaintySearchPlan{[](double x){ return(x-1); },
								[](double x){ return(x+1); }, tol});
	ps.addParameter("z");
	ps.setParameterProperty("z","UncertaintySearchPlan",
							UncertaintySearchPlan{[](double x){ return(x-1); },
								[](double x){ return(x+1); }, tol});
	ps.setParameterValue("x",1.5);
	ps.fixParameter("x");

	likelihoodPoint bestFit{{1.5,2,1},0};
	bestFit.likelihood=prob.llh(bestFit.params);

	auto f=std::function<likelihoodPoint(const problem&,const ParameterSet&)>(&fit);
	auto results=findParameterUncertainties(prob,ps,f,bestFit);
	//for(auto result : results)
	//	std::cout << '[' << result.first << ',' << result.second << "]\n";

	//first parameter is fixed
	if(results[0].first!=ps.getParameterValue(0))
		std::cout << "Left endpoint for parameter 0 should be fixed value\n";
	if(results[0].second!=ps.getParameterValue(0))
		std::cout << "Right endpoint for parameter 0 should be fixed value\n";
	//second parameter is free
	if(std::abs(results[1].first-(problem::ym-problem::ys))>tol)
		std::cout << "Left endpoint for parameter 1 not within tolerance of theoretical value: "
		<< results[1].first << " instead of " << (problem::ym-problem::ys) << std::endl;
	if(std::abs(results[1].second-(problem::ym+problem::ys))>tol)
		std::cout << "Right endpoint for parameter 1 not within tolerance of theoretical value: "
		<< results[1].first << " instead of " << (problem::ym+problem::ys) << std::endl;
	//third parameter is free
	if(std::abs(results[2].first-(problem::zm-problem::zs))>tol)
		std::cout << "Left endpoint for parameter 1 not within tolerance of theoretical value: "
		<< results[2].first << " instead of " << (problem::zm-problem::zs) << std::endl;
	if(std::abs(results[2].second-(problem::zm+problem::zs))>tol)
		std::cout << "Right endpoint for parameter 1 not within tolerance of theoretical value: "
		<< results[2].first << " instead of " << (problem::zm+problem::zs) << std::endl;
}
