#ifndef LBFGSB_INTERFACE_H
#define LBFGSB_INTERFACE_H

#include <cassert>
#include <cstdint>
#include <vector>
#include "PhysTools/autodiff.h"
#include "PhysTools/optimization/ParameterSet.h"

namespace phys_tools{
///The LBFGSB algorithm for function minimization with gradients
namespace lbfgsb{

///Base class for all functions which can be minimized with LBFGSB_Driver.
class BFGS_FunctionBase{
public:
	///Evaluate the value of the function at x
	virtual double evalF(std::vector<double> x) const=0;
	///Evaluate the value of the function and its gradient at x
	virtual std::pair<double,std::vector<double>> evalFG(std::vector<double> x) const=0;
};

///Base class for common cases of functions which can be minimized with LBFGSB_Driver.
///
/// This adapter class allows any subtype which implements
/// \code
///  template<typename T>
///  T operator()(T x1, T x2,..., T xN) const;
/// \endcode
/// with any number of parameters, which is assumed to be an automatic
/// differentiation-ready evaluation function, to be easily used with LBFGSB_Driver.
/// For example:
/// \code
///  struct Rosenbrock : public lbfgsb::SimpleBFGSFunction<Rosenbrock,2>{
///      template<typename T>
///      T operator()(T x, T y) const{
///          using std::sin;
///          using std::exp;
///          return(sin(x)*exp(x)+1/x);
///      }
///  };
/// \endcode
/// This is useful for simple functions with small numbers of parameters,
/// known at compile time.
/// Note, however, that using makeSimpleBFGSFunction may be simpler still.
///\tparam Derived the derived type which implements some function to be minimized
///\tparam Arity the number of parameters taken by Derived
template<typename Derived, unsigned Arity>
struct SimpleBFGSFunction : public BFGS_FunctionBase{
	static_assert(Arity>0,"There is no need to search for the minimum of something with no parameters");
private:
	//type for unpacking a runtime vector into compile-time function parameters
	template<unsigned I, typename T>
	struct computer{
		//wacky unused parameter is for overload selection
		static double convertArgType(double, double x, unsigned int index){ return(x); }
		//make sure each FD gets a assigned the correct index
		static phys_tools::autodiff::FD<Arity>
		convertArgType(phys_tools::autodiff::FD<Arity>, double x, unsigned int index){
			return(phys_tools::autodiff::FD<Arity>(x,index));
		}

		template<typename... Args>
		T invokeCompute(const Derived& f, const std::vector<double>& x, Args... args){
			return computer<I-1,T>().invokeCompute(f,x,convertArgType(T(),x[I-1],I-1),args...);
		}
	};
	template<typename T>
	struct computer<0,T>{
		template<typename... Args>
		T invokeCompute(const Derived& f, const std::vector<double>& x, Args... args){
			return f(args...);
		}
	};
public:
	virtual double evalF(std::vector<double> x) const{
		assert(x.size()==Arity);
		return computer<Arity,double>().invokeCompute(*(const Derived*)this,x);
	}
	virtual std::pair<double,std::vector<double>> evalFG(std::vector<double> x) const{
		assert(x.size()==Arity);
		using DType=phys_tools::autodiff::FD<Arity>;
		DType y=computer<Arity,DType>().invokeCompute(*(const Derived*)this,x);
		//reuse existing storage
		std::vector<double> g=std::move(x);
		for(unsigned int i=0; i<Arity; i++)
			g[i]=y.derivative(i);
		return(std::make_pair(y.value(),std::move(g)));
	}
};

namespace detail{
template<unsigned int Arity, typename FuncType>
struct simpleFunctionWrapper : public SimpleBFGSFunction<simpleFunctionWrapper<Arity,FuncType>,Arity>{
	FuncType f;
	simpleFunctionWrapper(FuncType&& f):f(f){}
	template<typename... Args, typename T=typename std::common_type<Args...>::type>
	T operator()(Args... args) const{
		return(f.template operator()<T>(args...));
	}
	operator FuncType&(){ return(f); }
	operator FuncType&() const{ return(f); }
};
} //namespace detail

///Adapts any suitably callable object to be minimized with LBFGSB_Driver.
///
/// Like SimpleBFGSFunction this function makes it easy to adapt an object which
/// implements
/// \code
///  template<typename T>
///  T operator()(T x1, T x2,..., T xN) const;
/// \endcode
/// with any number of parameters, which is assumed to be an automatic
/// differentiation-ready evaluation function, to be easily used with LBFGSB_Driver.
/// However, unlike SimpleBFGSFunction, it imposes no further requirements on that
/// object (it need not inherit from any particular base class). This can reduce
/// boilerplate code, and allows C++14 generic lambdas to be minimized.
/// Example of a plain function object type:
/// \code
///  struct Rosenbrock{
///      template<typename T>
///      T operator()(T x, T y) const{
///          using std::sin;
///          using std::exp;
///          return(sin(x)*exp(x)+1/x);
///      }
///  };
///  // ...
///  LBFGSB_Driver driver;
///  driver.addParameter(-3);
///  driver.addParameter(-4);
///  driver.minimize(makeSimpleBFGSFunction<2>(Rosenbrock()));
/// \endcode
/// Or, a generic lambda:
/// \code
///  const double a=1, b=100;
///  driver.minimize(makeSimpleBFGSFunction<2>([=](auto x, auto y){
///     decltype(x) d1=a-x;
///     decltype(x) d2=y-x*x;
///     return(d1*d1+b*d2*d2);
///  }));
/// \endcode
/// \tparam Arity the arity of the function being adapted
template<unsigned int Arity, typename FuncType>
detail::simpleFunctionWrapper<Arity,FuncType>
makeSimpleBFGSFunction(FuncType&& f){
	return(detail::simpleFunctionWrapper<Arity,FuncType>{std::forward<FuncType>(f)});
}

///Driver class for function minimization
class LBFGSB_Driver{
private:
	//variables related to the full set of function parameters specified by the user
	///Fit parameters
	const ParameterSet& parameters;
	///Current values of parameters
	std::vector<double> parameterValues;

	//variables related to the minimizer itself

	///number of previous gradients to store for the hessian approximation
	size_t historySize;
	///tolerance for absolute change in function values; if the change in
	///value between iterations is smaller than this value, the minimization
	///is considered to have converged
	double changeTol;
	///tolerance for magnitude of the gradient; if after an iteration the
	///largest component of the scaled function gradient (the gradient
	///divided element-wise by parameterScales) the minimization is considered
	///to have converged
	double gradTol;

	///records the value of the objective function at the best point found by minimizing
	///not valid until after minimize() has been called at least once!
	double finalFunctionValue;
	///Number of function evaluations performed during the last call to minimize()
	size_t nEvals;

    std::string lastStatus;

public:
	LBFGSB_Driver(const ParameterSet& ps):
	parameters(ps),
	parameterValues(ps.getParameterValues()),
	historySize(10),
	changeTol(1e-8),
	gradTol(1e-8),
	finalFunctionValue(std::numeric_limits<double>::quiet_NaN()),
	nEvals(0)
	{}

	///Set the number of past gradient values which should be used to estimate the hessian
	void setHistorySize(unsigned int size){
		historySize=size;
	}

	///Set the termination tolerance threshold for changes in successive values of the objective function
	void setChangeTolerance(double tol){
		changeTol=tol;
	}
	///Set the termination tolerance threshold for the maximum component of the scaled gradient
	void setGradientTolerance(double tol){
		gradTol=tol;
	}

	///Perform a minimization of a function with respect to the free parameters
	///\param func the objective function to be minimized
	///\returns true if minimization successfully converged
	///\post the values of all free parameters will be modified to match the location of the best point
	///the minimizer was able to find (even if the minimization was not deemed to have converged), and
	///can be retrieved by GetParameterValue() or minimumPosition(), and the function value at that point
	///can be retrieved using minimumValue(). Additionally, numberOfEvaluations() can report the number of
	///function evaluations performed during the minimization attempt.
	bool minimize(const BFGS_FunctionBase& func){
		//number of currently free variables
		int effectiveNVar=parameters.numberOfFreeParameters();

		if(effectiveNVar==0){ //the 0D case is so easy we don't have to call a library
			finalFunctionValue=func.evalF(parameterValues);
			nEvals=1;
			lastStatus.clear();
			return(true);
		}
		assert(effectiveNVar>0);

		///scaling factors for the parameters
		std::vector<double> parameterScales(effectiveNVar,1.);
		//get any scales specified by the user, keeping ones when not specified
		for(std::size_t i=0; i<effectiveNVar; i++){
			std::size_t idx=parameters.getFreeParameterIndex(i);
			if(parameters.parameterHasProperty(idx,"GralTolScale"))
				parameterScales[i]=parameters.getParameterProperty<double>(idx,"GralTolScale");
		}
		///the lower bounds for each variable
		std::vector<double> lowerBounds=parameters.getFreeParameterLowerBounds();
		///the upper bounds for each variable
		std::vector<double> upperBounds=parameters.getFreeParameterUpperBounds();
		///defines which bounds are used for each variable
		/// 0 : variable i is unbounded
		/// 1 : variable i is bounded below only
		/// 2 : variable i is bounded below and above
		/// 3 : variable i is bounded above only
		std::vector<int> boundsTypes(effectiveNVar,0);
		for(std::size_t i=0; i<effectiveNVar; i++){
			if(lowerBounds[i]>-std::numeric_limits<double>::infinity())
				boundsTypes[i]=1;
			if(upperBounds[i]<std::numeric_limits<double>::infinity())
				boundsTypes[i]=(boundsTypes[i]==1 ? 2 : 3);
		}

		///buffer for receiving instructions from setulb
		char task[60];
		///buffer for internal calculations (iwa)
		std::vector<int> internalBuf(3*effectiveNVar);
		///buffer for internal gradient history calculations (wa)
		std::vector<double> internalGradBuf((2*historySize + 5)*effectiveNVar + 11*historySize*historySize + 8*historySize);
		///buffer for internal string manipulation
		char textBuf[60];
		///buffer for integer diagnostics
		int diagnosticBuf[48];
		///buffer for floating point diagnostics
		double diagnosticBuf2[29];
		int iprint=0; //dummy var
		int histSize=historySize;

		//position known to the minimizer, contains only free variables
		std::vector<double> pos(effectiveNVar);
		double fval;
		//complete gradient with all variables
		std::vector<double> fullGrad;
		//gradient known to the minimizer, contains only free variables
		std::vector<double> grad(effectiveNVar);

		//prepare for first iteration
		memset(task, 0, 60);
		strncpy(task, "START", 5);
		parameters.extractFreeParameters(parameterValues,pos);
		nEvals=0;

		bool done=false;

		while(!done){
			setulb_(&effectiveNVar,     //n
					&histSize,          //m
					&pos[0],            //x
					&lowerBounds[0],    //l
					&upperBounds[0],    //u
					&boundsTypes[0],    //nbd
					&parameterScales[0],//scale
					&fval,              //f
					&grad[0],           //g
					&changeTol,         //factr
					&gradTol,           //pgtol
					&internalGradBuf[0],//wa
					&internalBuf[0],    //iwa
					&task[0],           //task
					&iprint,            //iprint, unused
					&textBuf[0],        //csave
					&diagnosticBuf[0],  //lsave
					&diagnosticBuf[4],  //isave
					&diagnosticBuf2[0]  //dsave
					);

			switch(task[0]){
				case 'A': //failure
					//copy back fit point and value
					parameters.insertFreeParameters(pos,parameterValues);
					finalFunctionValue=fval;
					done=true;
					//failure indicated to the user by return value below
					break;
				case 'C': //converged
					//copy back fit point and value
					parameters.insertFreeParameters(pos,parameterValues);
					finalFunctionValue=fval;
					done=true;
					break;
				case 'E':
					lastStatus=task;
					throw std::runtime_error("Invalid parameters passed to setulb");
				case 'F': //need to evaluate function
					parameters.insertFreeParameters(pos,parameterValues);
					std::tie(fval,fullGrad)=func.evalFG(parameterValues);
					parameters.extractFreeParameters(fullGrad,grad);
					nEvals++;
					break;
				case 'N': //iteration completed
					//check iteration count/number of evaluations?
					break;
			}
		}
        bool success=(task[0]=='C');
		lastStatus=task;
		return(success);
	}

    std::string getLastStatus() const{
        return(lastStatus);
    }

	///Get the smallest function value found by the minimizer
	///\pre minimize() must have been called previously
	double minimumValue() const{
		return(finalFunctionValue);
	}

	///Get the parametervalues at which the smallest function value
	///was found by the minimizer
	///\pre minimize() must have been called previously
	std::vector<double> minimumPosition() const{
		return(parameterValues);
	}

	///Get the number of function (and gradient) evaluations
	///performed during the last minimization attempt
	///\pre minimize() must have been called previously
	size_t numberOfEvaluations() const{
		return(nEvals);
	}
};

} //namespace lbfgsb
} //namespace phys_tools

#endif
