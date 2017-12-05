#include <iostream>
#include "../PhysTools/autodiff.h"
#include "../PhysTools/lbfgsb/lbfgsb.h"

using namespace phys_tools::lbfgsb;

struct absv{
	template<typename T>
	T operator()(T x) const{
		using std::abs;
		return(abs(x));
	}
	
	double minimumPosition() const{ return(0); }
	double minimumValue() const{ return(0); }
};

struct quadratic{
	double a,b,c;
	quadratic(double a, double b, double c):a(a),b(b),c(c){}

	template<typename T>
	T operator()(T x) const{
		return((a*x + b)*x + c);
	}
	
	double minimumPosition() const{
		//y=a*x^2+b*x+c
		//dy/dx=2*a*x+b
		//dy/dx=0 => -b/(2*a)=x
		return(-b/(2*a));
	}
	double minimumValue() const{ return(operator()(minimumPosition())); }
};

struct sine{
	template<typename T>
	T operator()(T x) const{
		using std::sin;
		return(sin(x));
	}
	
	double minimumPosition() const{
		return(6*atan(1.)); //3 pi/2
	}
	double minimumValue() const{ return(-1); }
};

struct mess{
	template<typename T>
	T operator()(T x) const{
		using std::sin;
		using std::exp;
		return(sin(x)*exp(x)+1/x);
	}
};

int main(){
	std::cout.precision(10);
	const double xtol=1e-6;
	const double ftol=1e-10;
	const double gtol=1e-10;
	//make sure we can handle a trivial case where the quadratic assumptions in
	//the minimzer will work perfectly
	{
		std::cout << "quadratic:\n";
		quadratic q{3.2,1.2,-7};
		
		LBFGSB_Driver driver;
		driver.addParameter(26.7);
		driver.setChangeTolerance(0); //use only gradient stopping condition
		driver.setGradientTolerance(gtol);
		bool success=driver.minimize(makeSimpleBFGSFunction<1>(q));
		if(!success){
			std::cout << " minimization failed!" << std::endl;
			std::cout << ' ' << driver.errorMessage() << std::endl;
		}
		auto min=driver.minimumPosition();
		auto trueMin=q.minimumPosition();
		if(std::abs(min[0]-trueMin)<xtol)
			std::cout << " x within tolerance of truth" << std::endl;
		else
			std::cout << " x not within tolerance of truth: "
			<< "true value=" << trueMin
			<< " obtained value=" << min[0] << std::endl;
		if(driver.minimumValue()<q.minimumValue()+ftol)
			std::cout << " f within tolerance of truth" << std::endl;
		else
			std::cout << " f not within tolerance of truth: "
			<< "true value=" << q.minimumValue()
			<< " obtained value=" << driver.minimumValue() << std::endl;
		auto result=q(phys_tools::autodiff::FD<1>(driver.minimumPosition().front(),0));
		if(std::abs(result.derivative(0))<gtol)
			std::cout << " f' is within tolerance" << std::endl;
		else
			std::cout << " f' is not within tolerance: " << result.derivative(0) << std::endl;
	}
	//make sure we can handle a not completely trivial case
	{
		std::cout << "sine:\n";
		sine s;
		
		LBFGSB_Driver driver;
		driver.addParameter(3,1,2*atan(1),10*atan(1)); //search within [pi/2,5 pi/2]
		driver.setChangeTolerance(0);
		driver.setGradientTolerance(gtol); //use only gradient stopping condition
		bool success=driver.minimize(makeSimpleBFGSFunction<1>(s));
		if(!success){
			std::cout << " minimization failed!" << std::endl;
			std::cout << ' ' << driver.errorMessage() << std::endl;
		}
		auto min=driver.minimumPosition();
		auto trueMin=s.minimumPosition();
		if(std::abs(min[0]-trueMin)<xtol)
			std::cout << " x within tolerance of truth" << std::endl;
		else
			std::cout << " x not within tolerance of truth: "
			<< "true value=" << trueMin
			<< " obtained value=" << min[0] << std::endl;
		if(driver.minimumValue()<s.minimumValue()+ftol)
			std::cout << " f within tolerance of truth" << std::endl;
		else
			std::cout << " f not within tolerance of truth: "
			<< "true value=" << s.minimumValue()
			<< " obtained value=" << driver.minimumValue() << std::endl;
		auto result=s(phys_tools::autodiff::FD<1>(driver.minimumPosition().front(),0));
		if(std::abs(result.derivative(0))<gtol)
			std::cout << " f' is within tolerance" << std::endl;
		else
			std::cout << " f' is not within tolerance: " << result.derivative(0) << std::endl;
	}
	//we should be able to blunder to some sort of result for a continuous, convex
	//function whose derivative is discontinuous
	{
		std::cout << "abs:\n";
		absv abs;
		
		LBFGSB_Driver driver;
		driver.addParameter(16);
		driver.setChangeTolerance(1e-10); //use only change stooping condition
		driver.setGradientTolerance(0);
		bool success=driver.minimize(makeSimpleBFGSFunction<1>(abs));
		//!!! Omit this check because LBFGSB does not really like this function
		//if(!success){
		//	std::cout << " minimization failed!" << std::endl;
		//	std::cout << ' ' << driver.errorMessage() << std::endl;
		//}
		auto min=driver.minimumPosition();
		auto trueMin=abs.minimumPosition();
		if(std::abs(min[0]-trueMin)<xtol)
			std::cout << " x within tolerance of truth" << std::endl;
		else
			std::cout << " x not within tolerance of truth: "
			<< "true value=" << trueMin
			<< " obtained value=" << min[0] << std::endl;
		if(driver.minimumValue()<abs.minimumValue()+ftol)
			std::cout << " f within tolerance of truth" << std::endl;
		else
			std::cout << " f not within tolerance of truth: "
			<< "true value=" << abs.minimumValue()
			<< " obtained value=" << driver.minimumValue() << std::endl;
		//!!! Omit gradient check, since this function's gradient never goes to zero
	}
	//make sure we can handle some ugly mess
	{
		std::cout << "mess:\n";
		mess m;
		
		LBFGSB_Driver driver;
		driver.addParameter(2,1,1e-4,2.3);
		driver.setChangeTolerance(0); //use only gradient stopping condition
		driver.setGradientTolerance(gtol);
		bool success=driver.minimize(makeSimpleBFGSFunction<1>(m));
		if(!success){
			std::cout << " minimization failed!" << std::endl;
			std::cout << ' ' << driver.errorMessage() << std::endl;
		}
		//settle for just checking the function's gradient
		auto result=m(phys_tools::autodiff::FD<1>(driver.minimumPosition().front(),0));
		if(std::abs(result.derivative(0))<gtol)
			std::cout << " f' is within tolerance" << std::endl;
		else
			std::cout << " f' is not within tolerance: " << result.derivative(0) << std::endl;
	}
}