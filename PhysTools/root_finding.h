#ifndef PHYSTOOLS_ROOT_FINDING_H
#define PHYSTOOLS_ROOT_FINDING_H

#include <iostream>

namespace phys_tools{

//search for the the value x such that func(c)==target where c is in [a,b]
//and fa=func(a), fb=func(b), assuming that func is monotonic on this domain.
//This is some sort of simplified form of Brent's method.
template<typename FuncType>
double findRoot(FuncType func, double target,
                double a, double b,
                double fa, double fb,
                double precision=1e-6, bool verbose=false){
	bool increasing=fb>fa;
	fa-=target;
	fb-=target;

	double c=(a+b)/2; //begin with a bisection
	double fc, e;

	while((b-a)>precision){
		fc=func(c)-target;
		if(verbose)
			std::cout << "  f(" << c << ")=" << fc+target << std::endl;

		//calculate the quadratic interpolation
		double r=fc/fb;
		double s=fc/fa;
		double t=fa/fb;
		double p=s*(t*(r-t)*(b-c)+(r-1)*(c-a));
		double q=(r-1)*(s-1)*(t-1);

		//std::cout << "  " << fa << ' ' << fb << ' ' << fc << '\n';
		//std::cout << "  " << r << ' ' << s << ' ' << t << ' ' << p << ' ' << q << '\n';

		//collapse the interval
		if((fc>0.0) != increasing){
			e=(c-a)/(b-a);
			a=c;
			fa=fc;
		}
		else{
			e=(b-c)/(b-a);
			b=c;
			fb=fc;
		}

		if(verbose){
			std::cout << " interval now [" << a << ',' << b << ']' << std::endl;
			std::cout << " interpolated guess = " << c+p/q << std::endl;
		}

		//accept the interpolation only if it falls within the current boundaries
		if((c+(p/q))>=a && (c+(p/q))<=b){
			//would like to interpolate, but only do so if recent convergence has been suitably rapid
			if(e<0.5){ //it has not; bisect instead
				if(verbose)
					std::cout << " will bisect (too slow)" << std::endl;
				c=(a+b)/2;
			}
			else{ //it has; use the interpolation
				if(verbose)
					std::cout << " will interpolate" << std::endl;
				c+=p/q;
			}
		}
		else{ //otherwise, bisect
			if(verbose)
				std::cout << " will bisect (out-of bounds)" << std::endl;
			c=(a+b)/2;
		}
	}
	return(c);
}

struct bracketResult{
	double xin, fin;
	double xout, fout;
};

//Given a function which is known to have a root in the domain [a,limit) or (limit,a]
//and an initial guess b which is on the same side of a as the root (>a or <a, respectively)
//finds an interval which strictly brackets the root.
//The return value is a pair consisting of the other endpoint of the interval, and the function value there
///\param a a value of the indepedent variable known to be on one side of the root
///\param b a value of the indepedent variable which is in the direction of the
///       root from a, and may be on the other side
///\param fa the value of func(a)
///\param limit the most extreme value which should be considered when searching
///       for the outer end of the bracketing interval, or NaN if +/- infinity
///       should be automatically selected
///\param verbose whether to print debugging information to stdout
template<typename FuncType>
bracketResult bracketRoot(FuncType func, double target,
                                     double a, double b,
                                     double fa,
                                     double limit=std::numeric_limits<double>::quiet_NaN(),
                                     bool verbose=false){
	bool dir=b>a;
	if(std::isnan(limit))
		limit=(dir?1:-1)*std::numeric_limits<double>::infinity();
	if(verbose){
		std::cout << " bracketing domain is [" << a << ',' << limit << ')' << std::endl;
		std::cout << " guess is " << b << std::endl;
	}
	fa-=target;
	double c=a, fc=fa;
	double clast, fclast;
	double fact=1;
	while((fc>0)==(fa>0)){
		clast=c;
		fclast=fc;
		c=a+fact*(b-a);
		if((dir && c>limit) || (!dir && c<limit)){
			if(verbose)
				std::cout << " clipping endpoint to limit" << std::endl;
			c=limit;
		}
		if(verbose)
			std::cout << " testing endpoint " << c << std::endl;
		fc=func(c)-target;
		if(verbose)
			std::cout << "  (" << fc << ')' << std::endl;
		if(c==limit)
			break;
		fact*=2;
	}
	//return(std::make_pair(c,fc+target));
	bracketResult r;
	r.xin=clast;
	r.fin=fclast+target;
	r.xout=c;
	r.fout=fc+target;
	if(verbose){
		std::cout << " bracketing recommends interval: [" << r.xin << ',' << r.xout
		<< "] (function values " << r.fin << ", " << r.fout << ")\n";
	}
	return(r);
}

} //namespace phys_tools

#endif
