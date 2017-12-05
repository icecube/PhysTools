#include <PhysTools/histogram.h>
#include <iostream>

struct valueWithError{
	double val, err;
	valueWithError():val(0),err(0){}
	valueWithError(double v):val(v),err(sqrt(v)){}
	valueWithError(double v, double e):val(v),err(e){}
	valueWithError& operator+=(const valueWithError& other){
		val+=other.val;
		err=sqrt(err*err+other.err*other.err);
		return(*this);
	}
	valueWithError& operator*=(double a){
		val*=a;
		err*=a;
		return(*this);
	}
	valueWithError& operator/=(double a){
		val/=a;
		err/=a;
		return(*this);
	}
};

valueWithError operator*(const valueWithError& v, double a){
	return(valueWithError(v)*=a);
}
valueWithError operator*(double a, const valueWithError& v){
	return(valueWithError(v)*=a);
}
valueWithError operator/(const valueWithError& v, double a){
	return(valueWithError(v)/=a);
}

std::ostream& operator<<(std::ostream& os, const valueWithError& v){
	return(os << v.val << '[' << v.err << ']');
}

namespace phys_tools{
	namespace histograms{
		namespace detail{
			template<>
			class histogramTraits<valueWithError>{
			public:
				using amount=phys_tools::histograms::detail::amount_type<valueWithError>;
				static void defaultData(valueWithError* data, unsigned int count){/* Nothing to do */}
				static double unit(){ return 1.0; }
				constexpr static bool enable_automatic_amount_handling=true;
			};
		}
	}
}

int main(){
	using namespace phys_tools::histograms;
	
	histogram<2> h1(LinearAxis(0,1),LinearAxis(0,1));
	h1.add(0,0,amount(1));
	h1.add(0,1,amount(2));
	h1.add(1,0,amount(3));
	h1.add(1,1,amount(4));
	std::cout << h1 << std::endl;
	
	h1/=10.;
	std::cout << h1 << std::endl;
	
	h1*=5.;
	std::cout << h1 << std::endl;
	
	histogram<2,valueWithError> h2(LinearAxis(0,1),LinearAxis(0,1));
	h2.add(0,0,amount(valueWithError(1)));
	h2.add(0,1,amount(valueWithError(2)));
	h2.add(1,0,amount(valueWithError(3)));
	h2.add(1,1,amount(valueWithError(4)));
	std::cout << h2 << std::endl;
	
	h2/=10.;
	std::cout << h2 << std::endl;
	
	h2*=5.;
	std::cout << h2 << std::endl;
}