#include <iostream>
#include "../PhysTools/autodiff.h"
#include "../PhysTools/optimization/lbfgsb/lbfgsb.h"

using namespace phys_tools::lbfgsb;

struct Rosenbrock{
	double a, b;
	Rosenbrock(double a, double b):a(a),b(b){}

	template<typename T>
	T operator()(T x, T y) const{
		T d1=a-x;
		T d2=y-x*x;
		return(d1*d1+b*d2*d2);
	}

	//related handy property: the minimum function value is always 0
	std::pair<double,double> minimum() const{
		return(std::make_pair(a,a*a));
	}
};

int main(){
	std::cout.precision(10);
	const double xtol=1e-6;
	const double ftol=1e-10;

	//test a classic
	{
		std::cout << "R1:\n";
		Rosenbrock r(1,100);
		phys_tools::ParameterSet params;
		params.addParameter("x");
		params.setParameterValue("x",-3);
		params.addParameter("y");
		params.setParameterValue("x",-4);
		LBFGSB_Driver driver(params);
		driver.setChangeTolerance(1e-10);
		driver.setGradientTolerance(1e-10);
		bool success=driver.minimize(makeSimpleBFGSFunction<2>(r));
		auto min=driver.minimumPosition();
		//std::cout << " found minimum at (" << min[0] << ',' << min[1] << ')' << std::endl;
		//std::cout << " minimum value " << driver.minimumValue() << std::endl;
		auto trueMin=r.minimum();
		if(std::abs(min[0]-trueMin.first)<xtol)
			std::cout << "x within tolerance of truth" << std::endl;
		else
			std::cout << "x not within tolerance of truth: "
			<< "true value=" << trueMin.first
			<< " obtained value=" << min[0] << std::endl;
		if(std::abs(min[1]-trueMin.second)<xtol)
			std::cout << "y within tolerance of truth" << std::endl;
		else
			std::cout << "y not within tolerance of truth: "
			<< "true value=" << trueMin.second
			<< " obtained value=" << min[1] << std::endl;
		if(driver.minimumValue()<ftol) //true value is zero, function is non-negative
			std::cout << "f within tolerance of truth" << std::endl;
		else
			std::cout << "f not within tolerance of truth: "
			<< "true value=0"
			<< " obtained value=" << driver.minimumValue() << std::endl;
	}
	{
		std::cout << "R2:\n";
		Rosenbrock r(2.2,14);
		phys_tools::ParameterSet params;
		params.addParameter("x");
		params.setParameterValue("x",7);
		params.addParameter("y");
		params.setParameterValue("x",-4);
		LBFGSB_Driver driver(params);
		driver.setChangeTolerance(1e-10);
		driver.setGradientTolerance(1e-10);
		bool success=driver.minimize(makeSimpleBFGSFunction<2>(r));
		auto min=driver.minimumPosition();
		//std::cout << " found minimum at (" << min[0] << ',' << min[1] << ')' << std::endl;
		//std::cout << " minimum value " << driver.minimumValue() << std::endl;
		auto trueMin=r.minimum();
		if(std::abs(min[0]-trueMin.first)<xtol)
			std::cout << "x within tolerance of truth" << std::endl;
		else
			std::cout << "x not within tolerance of truth: "
			<< "true value=" << trueMin.first
			<< " obtained value=" << min[0] << std::endl;
		if(std::abs(min[1]-trueMin.second)<xtol)
			std::cout << "y within tolerance of truth" << std::endl;
		else
			std::cout << "y not within tolerance of truth: "
			<< "true value=" << trueMin.second
			<< " obtained value=" << min[1] << std::endl;
		if(driver.minimumValue()<ftol) //true value is zero, function is non-negative
			std::cout << "f within tolerance of truth" << std::endl;
		else
			std::cout << "f not within tolerance of truth: "
			<< "true value=0"
			<< " obtained value=" << driver.minimumValue() << std::endl;
	}
	/*{
		std::cout << "R3:\n";
		const double a=1, b=100;
		phys_tools::ParameterSet params;
		params.addParameter("x");
		params.setParameterValue("x",-3);
		params.addParameter("y");
		params.setParameterValue("x",-4);
		LBFGSB_Driver driver(params);
		driver.setChangeTolerance(1e-10);
		driver.setGradientTolerance(1e-10);
		bool success=driver.minimize(makeSimpleBFGSFunction<2>([=](auto x, auto y){
			decltype(x) d1=a-x;
			decltype(x) d2=y-x*x;
			return(d1*d1+b*d2*d2);
		}));
		auto min=driver.minimumPosition();
		auto trueMin=std::make_pair(a,a*a);
		if(std::abs(min[0]-trueMin.first)<xtol)
			std::cout << "x within tolerance of truth" << std::endl;
		else
			std::cout << "x not within tolerance of truth: "
			<< "true value=" << trueMin.first
			<< " obtained value=" << min[0] << std::endl;
		if(std::abs(min[1]-trueMin.second)<xtol)
			std::cout << "y within tolerance of truth" << std::endl;
		else
			std::cout << "y not within tolerance of truth: "
			<< "true value=" << trueMin.second
			<< " obtained value=" << min[1] << std::endl;
		if(driver.minimumValue()<ftol) //true value is zero, function is non-negative
			std::cout << "f within tolerance of truth" << std::endl;
		else
			std::cout << "f not within tolerance of truth: "
			<< "true value=0"
			<< " obtained value=" << driver.minimumValue() << std::endl;
	}*/
}
