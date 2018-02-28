// Copyright (c) 2014-2017, Christopher Weaver
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef AUTODIFF_H_INCLUDED
#define AUTODIFF_H_INCLUDED

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <iterator>

#include <boost/math/constants/constants.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

namespace phys_tools{
///Automatic differentiation
namespace autodiff{
	
	constexpr unsigned int Dynamic=0;
	
	//fwd decl
	template <unsigned int nVars, typename T>
	class FD;
	
	namespace detail{
		
		//Provide a uniform mechanism for assessing the dimensionality of a FD's gradient
		//For fixed-size FDs, just return the constant size
		template <template<unsigned int,class> class F, int D, class T>
		struct dimensionExtractor{
			static unsigned int nVars(const F<D,T>&){
				return(F<D,T>::N);
			}
		};
		
		//For dynamically-sized FDs, actually query the object's current dimensionality
		template <template<unsigned int,class> class F, class T>
		struct dimensionExtractor<F,Dynamic,T>{
			static unsigned int nVars(const F<Dynamic,T>& f){
				return(f.nVars);
			}
		};
		
	}
	
	///\brief A type for forward-mode automatic differentiation
	///
	///If the number of independent variables derivatives (or better, the number
	///of independent variables for which the derivatives of the computed values
	///are interesting) is known at compile time, this should be supplied as the
	///first template parameter for this class, allowing dynamic memory allocation
	///to be avoided. This can make computation more efficient by up to a factor
	///of 20-30. Otherwise, the special value `Dynamic` may be used, which casues
	///storage for gradient vectors to be allocated on the heap.
	///
	///Each of the independent variables for which partial derivatives are to be
	///evaluated must be assigned a unique position within the gradient vector,
	///using assignIndex() (on the independent variable). Finally, the partial
	///derivatives of a dependent variable with respect to the independent
	///variables are read out after the calculation using derivative() on the
	///dependent variable, where the index argument to derivative() is the index
	///assigned to the independent variable with respect to which the partial
	///derivative is desired.
	///
	///Example:
	/// Suppose we want to compute f(x,y) = 2*x*y + 4*x + 6*y, df/dx, and df/dy
	/// at x=7.22, y=9.4, in single precision.
	/// First, we create the independent variables, and assign each an index in
	/// the gradient vector:
	/// \code
	///  autodiff::FD<2,float> x(7.22);
	///  autodiff::FD<2,float> y(9.4);
	///  x.assignIndex(0);
	///  y.assignIndex(1);
	/// \endcode
	/// Or, more concisely:
	/// \code
	///  autodiff::FD<2,float> x(7.22,0);
	///  autodiff::FD<2,float> x(9.4,1);
	/// \endcode
	/// Then we calculate the result (and implicitly calculate the derivatives):
	/// \code
	///  autodiff::FD<2,float> f = 2.0f*x*y + 4.0f*x + 6.0f*z;
	/// \endcode
	/// We can then get the value and the derivatives:
	/// \code
	///  std::cout << "f = " << f.value() << std::endl;
	///  std::cout << "df/dx = " << f.derivative(0) << std::endl;
	///  std::cout << "df/dy = " << f.derivative(1) << std::endl;
	/// \endcode
	/// The same can be done even if the total number of independent variables is
	/// not known ahead of time (perhaps because this is part of some larger
	/// calculation for which the number of included variables can change):
	/// \code
	///  autodiff::FD<autodiff::Dynamic,float> x(7.22,0);
	///  autodiff::FD<autodiff::Dynamic,float> x(9.4,1);
	///  autodiff::FD<autodiff::Dynamic,float> f = 2.0f*x*y + 4.0f*x + 6.0f*z;
	///  std::cout << "f = " << f.value() << std::endl;
	///  std::cout << "df/dx = " << f.derivative(0) << std::endl;
	///  std::cout << "df/dy = " << f.derivative(1) << std::endl;
	/// \endcode
	/// Note that it is still necessary to assign unique indices to the independent
	/// variables, however, this works perfectly well with indices chosen at runtime.
	template <unsigned int nVars, typename T=double>
	class FD{
	private:
		T v; //value
		T g[nVars]; //gradient
		
		struct noInit_tag{};
		static constexpr noInit_tag noInit(){ return(noInit_tag{}); };
		
		explicit FD<nVars,T>(const noInit_tag&){}
		
		FD<nVars,T>(T t, const noInit_tag&):
		v(t){}
		
	public:
		typedef T BaseType;
		enum {N=nVars};
		
		FD<nVars,T>(){
			std::fill(g,g+nVars,0);
		}
		
		FD<nVars,T>(T t):
		v(t){
			std::fill(g,g+nVars,0);
		}
		
		FD<nVars,T>(T t, unsigned int index):
		v(t){
			assert(index<nVars);
			std::fill(g,g+nVars,0);
			g[index]=T(1);
		}
		
        FD<nVars,T>(const FD<nVars,T>& f):
		v(f.v){
			std::copy(f.g,f.g+nVars,g);
		}
		
        template<typename U>
		FD<nVars,T>(const FD<nVars,U>& f):
		v(f.value()){
            f.copyGradient(g);
		}
		
		FD<nVars,T>& operator=(const FD<nVars,T>& f){
			if(&f!=this){
				v=f.v;
				std::copy(f.g,f.g+nVars,g);
			}
			return(*this);
		}
		
        template<typename U>
        FD<nVars,T>& operator=(const FD<nVars,U>& f){
			v=f.value();

            f.copyGradient(g);
			return(*this);
		}
		
		FD<nVars,T>& operator=(const T& t){
			v=t;
			std::fill(g,g+nVars,0);
			return(*this);
		}
		
		void assignIndex(unsigned int index){
			assert(index<nVars);
			std::fill(g,g+nVars,0);
			g[index]=T(1);
		}
		
		const T& value() const{
			return(v);
		}
		
		const T& derivative(unsigned int index) const{
			assert(index<nVars);
			return(g[index]);
		}
		
		///Manually set a component of the derivative to a particular value.
		///This is mostly only useful when interfacing with other code which
		///somehow computes gradient informaton manually and one wishes to then
		///propagate it further using automatic differentiation.
		void setDerivative(unsigned int index, T d){
			assert(index<nVars);
			g[index]=d;
		}
		
        template<typename Iter>
        void copyGradient(Iter grad) const{
            using val_t=typename std::iterator_traits<Iter>::value_type;
            std::transform(g, g+nVars, grad, [](T c) -> val_t { return static_cast<val_t>(c); });
		}
		
		void copyGradient(T grad[nVars]) const{
			std::copy(g,g+nVars,grad);
		}
		
		//unary +
		FD<nVars,T> operator+() const{
			return(*this);
		}
		
		//unary -
		FD<nVars,T> operator-() const{
			FD<nVars,T> r(*this);
			r.v=-r.v;
			std::transform(r.g,r.g+nVars,r.g,std::negate<T>());
			return(r);
		}
		
		//addition
		template <typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T>& operator+=(const U& u){
			v+=T(u);
			return(*this);
		}
		
		FD<nVars,T>& operator+=(const FD<nVars,T>& f){
			v+=f.v;
			std::transform(g,g+nVars,f.g,g,std::plus<T>());
			return(*this);
		}
		
		template <typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T> operator+(const U& u) const{
			return(FD<nVars,T>(*this)+=u);
		}
		
		FD<nVars,T> operator+(const FD<nVars,T>& f) const{
			FD<nVars,T> result(v+f.v,noInit());
			std::transform(g,g+nVars,f.g,result.g,std::plus<T>());
			return(result);
		}
		
		//subtraction
		template <typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T>& operator-=(const U& u){
			v-=T(u);
			return(*this);
		}
		
		FD<nVars,T>& operator-=(const FD<nVars,T>& f){
			v-=f.v;
			std::transform(g,g+nVars,f.g,g,std::minus<T>());
			return(*this);
		}
		
		template <typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T> operator-(const U& u) const{
			return(FD<nVars,T>(*this)-=u);
		}
		
		FD<nVars,T> operator-(const FD<nVars,T>& f) const{
			FD<nVars,T> result(v-f.v,noInit());
			std::transform(g,g+nVars,f.g,result.g,std::minus<T>());
			return(result);
		}
		
		//multiplication
		template <typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T>& operator*=(const U& u){
			T t(u);
			v*=t;
			std::transform(g,g+nVars,g,[&](const T& gi){ return(gi*t); });
			return(*this);
		}
		
		FD<nVars,T>& operator*=(const FD<nVars,T>& f){
			for(unsigned int i=0; i<nVars; i++)
				g[i] = g[i]*f.v + f.g[i]*v;
			v*=f.v;
			return(*this);
		}
		
		template<typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T> operator*(const U& u) const{
			return(FD<nVars,T>(*this)*=u);
		}
		
		FD<nVars,T> operator*(const FD<nVars,T>& f) const{
			FD<nVars,T> result(v*f.v,noInit());
			for(unsigned int i=0; i<nVars; i++)
				result.g[i] = g[i]*f.v + f.g[i]*v;
			return(result);
		}
		
		//division
		template<typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T>& operator/=(const U& u){
			T t(u);
			v/=t;
			std::transform(g,g+nVars,g,[&](const T& gi){ return(gi/t); });
			return(*this);
		}
		
		FD<nVars,T>& operator/=(const FD<nVars,T>& f){
			v/=f.v;
			for(unsigned int i=0; i<nVars; i++)
				g[i] = (g[i] - v*f.g[i])/f.v;
			return(*this);
		}
		
		template<typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<nVars,T> operator/(const U& u) const{
			return(FD<nVars,T>(*this)/=u);
		}
		
		FD<nVars,T> operator/(const FD<nVars,T>& f) const{
			FD<nVars,T> result(v/f.v,noInit());
			for(unsigned int i=0; i<nVars; i++)
				result.g[i] = (g[i] - result.v*f.g[i])/f.v;
			return(result);
		}
	
		template <template<unsigned int,class> class F, int D, class T_>
		friend struct detail::dimensionExtractor;
		template <unsigned int nVars_, typename T_, typename U>
		friend FD<nVars_,T_> operator/(const U& u, const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> cos(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> sin(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> tan(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> acos(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> asin(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> atan(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> cosh(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> sinh(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> tanh(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> acosh(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> asinh(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> atanh(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> exp(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> log(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> log10(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> log2(const FD<nVars_,T_>& f);
		template <typename FD_t, typename U>
		friend FD_t pow(const FD_t& b, const U& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type);
		template <typename FD_t, typename U>
		friend FD_t pow(const U& b, const FD_t& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type);
		template <typename FD_t, typename U>
		friend FD_t pow(const U& b, const FD_t& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> pow(const FD<nVars_,T_>& b, const FD<nVars_,T_>& e);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> sqrt(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> cbrt(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> floor(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> ceil(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> abs(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> fabs(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> tgamma(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> lgamma(const FD<nVars_,T_>& f);
	};
	
	template <unsigned int nVars, typename T, typename U>
	FD<nVars,T> operator+(const U& u, const FD<nVars,T>& f){
		return(f+u);
	}

	template <unsigned int nVars, typename T, typename U>
	FD<nVars,T> operator+(const U& u, const FD<nVars,T>&& f){
		return(std::move(f)+u);
	}

	template <unsigned int nVars, typename T, typename U>
	FD<nVars,T> operator-(const U& u, const FD<nVars,T>& f){
		return(-(FD<nVars,T>(f)-=u));
	}
	
	template <unsigned int nVars, typename T, typename U>
	FD<nVars,T> operator*(const U& u, const FD<nVars,T>& f){
		return(f.operator*(u));
	}

	template <unsigned int nVars, typename T, typename U>
	FD<nVars,T> operator*(const U& u, FD<nVars,T>&& f){
		return(std::move(f)*u);
	}

	template <unsigned int nVars, typename T, typename U>
	FD<nVars,T> operator/(const U& u, const FD<nVars,T>& f){
		FD<nVars,T> result(f);
		T t(u);
		result.v=t/result.v;
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		T m(-result.v/f.v);
		std::transform(result.g,result.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> cos(const FD<nVars,T>& f){
		using std::cos;
		using std::sin;
		FD<nVars,T> result(cos(f.v),FD<nVars,T>::noInit());
		const T m(-sin(f.v));
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}

	template <typename T>
	FD<Dynamic,T> cos(const FD<Dynamic,T>& f){
		using std::cos;
		using std::sin;
		FD<Dynamic,T> result(cos(f.v),f.nVars,FD<Dynamic,T>::noInit());
		const T m(-sin(f.v));
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}

	template <unsigned int nVars, typename T>
	FD<nVars,T> sin(const FD<nVars,T>& f){
		using std::cos;
		using std::sin;
		FD<nVars,T> result(sin(f.v),FD<nVars,T>::noInit());
		const T m(cos(f.v));
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> sin(const FD<Dynamic,T>& f){
		using std::cos;
		using std::sin;
		FD<Dynamic,T> result(sin(f.v),f.nVars,FD<Dynamic,T>::noInit());
		const T m(cos(f.v));
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> tan(const FD<nVars,T>& f){
		using std::tan;
		FD<nVars,T> result(tan(f.v),FD<nVars,T>::noInit());
		T m=T(1)+result.v*result.v;
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> tan(const FD<Dynamic,T>& f){
		using std::tan;
		FD<Dynamic,T> result(tan(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=T(1)+result.v*result.v;
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> acos(const FD<nVars,T>& f){
		using std::acos;
		using std::sqrt;
		FD<nVars,T> result(acos(f.v),FD<nVars,T>::noInit());
		T m=-T(1)/sqrt(T(1)-f.v*f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> acos(const FD<Dynamic,T>& f){
		using std::acos;
		using std::sqrt;
		FD<Dynamic,T> result(acos(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=-T(1)/sqrt(T(1)-f.v*f.v);
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> asin(const FD<nVars,T>& f){
		using std::asin;
		using std::sqrt;
		FD<nVars,T> result(asin(f.v), FD<nVars,T>::noInit());
		T m=T(1)/sqrt(T(1)-f.v*f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> asin(const FD<Dynamic,T>& f){
		using std::asin;
		using std::sqrt;
		FD<Dynamic,T> result(asin(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=T(1)/sqrt(T(1)-f.v*f.v);
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> atan(const FD<nVars,T>& f){
		using std::atan;
		FD<nVars,T> result(atan(f.v),FD<nVars,T>::noInit());
		T m=T(1)/(T(1)+f.v*f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> atan(const FD<Dynamic,T>& f){
		using std::atan;
		FD<Dynamic,T> result(atan(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=T(1)/(T(1)+f.v*f.v);
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> cosh(const FD<nVars,T>& f){
		using std::cosh;
		using std::sinh;
		FD<nVars,T> result(cosh(f.v),FD<nVars,T>::noInit());
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		T m(sinh(f.v));
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> cosh(const FD<Dynamic,T>& f){
		using std::cosh;
		using std::sinh;
		FD<Dynamic,T> result(cosh(f.v),f.nvars,FD<Dynamic,T>::noInit());
		T m(sinh(f.v));
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> sinh(const FD<nVars,T>& f){
		using std::cosh;
		using std::sinh;
		FD<nVars,T> result(sinh(f.v),FD<nVars,T>::noInit());
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		T m(cosh(f.v));
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> sinh(const FD<Dynamic,T>& f){
		using std::cosh;
		using std::sinh;
		FD<Dynamic,T> result(sinh(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=cosh(f.v);
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> tanh(const FD<nVars,T>& f){
		using std::tanh;
		FD<nVars,T> result(tanh(f.v),FD<nVars,T>::noInit());
		T m=1-result.v*result.v;
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> tanh(const FD<Dynamic,T>& f){
		using std::tanh;
		FD<Dynamic,T> result(tanh(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=1-result.v*result.v;
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> acosh(const FD<nVars,T>& f){
		using std::acosh;
		using std::sqrt;
		FD<nVars,T> result(acosh(f.v),FD<nVars,T>::noInit());
		T m=1/(sqrt(f.v+1)*sqrt(f.v-1));
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> acosh(const FD<Dynamic,T>& f){
		using std::acosh;
		using std::sqrt;
		FD<Dynamic,T> result(acosh(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=1/(sqrt(f.v+1)*sqrt(f.v-1));
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> asinh(const FD<nVars,T>& f){
		using std::asinh;
		using std::sqrt;
		FD<nVars,T> result(asinh(f.v),FD<nVars,T>::noInit());
		T m=1/(sqrt(f.v*f.v));
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> asinh(const FD<Dynamic,T>& f){
		using std::asinh;
		using std::sqrt;
		FD<Dynamic,T> result(asinh(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=1/(sqrt(f.v*f.v));
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> atanh(const FD<nVars,T>& f){
		using std::atanh;
		FD<nVars,T> result(atanh(f.v),FD<nVars,T>::noInit());
		T m=1/(1-f.v*f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> atanh(const FD<Dynamic,T>& f){
		using std::atanh;
		FD<Dynamic,T> result(atanh(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=1/(1-f.v*f.v);
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> exp(const FD<nVars,T>& f){
		using std::exp;
		FD<nVars,T> result(exp(f.v),FD<nVars,T>::noInit());
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g,[&](const T& gi){ return(gi*result.v); });
		return(result);
	}

	template <typename T>
	FD<Dynamic,T> exp(const FD<Dynamic,T>& f){
		using std::exp;
		const unsigned int n=detail::dimensionExtractor<FD,Dynamic,T>::nVars(f);
		FD<Dynamic,T> result(f,n,FD<Dynamic,T>::noInit());
		result.v=exp(f.v);
		std::transform(f.g,f.g+f.nVars,result.g,[&](const T& gi){ return(gi*result.v); });
		return(result);
	}

	template <typename T>
	FD<Dynamic,T> exp(FD<Dynamic,T>&& f){
		using std::exp;
		FD<Dynamic,T> result(exp(f.v),std::move(f));
		std::transform(result.g,result.g+result.nVars,result.g,[&](const T& gi){ return(gi*result.v); });
		return(result);
	}

	template <unsigned int nVars, typename T>
	FD<nVars,T> log(const FD<nVars,T>& f){
		using std::log;
		FD<nVars,T> result(log(f.v),FD<nVars,T>::noInit());
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g,[&](const T& gi){ return(gi/f.v); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> log(const FD<Dynamic,T>& f){
		using std::log;
		FD<Dynamic,T> result(log(f.v),f.nVars,FD<Dynamic,T>::noInit());
		std::transform(f.g,f.g+f.nVars,result.g,[&](const T& gi){ return(gi/f.v); });
		return(result);
	}

	//TODO: this may not be advantageous?
	template <typename T>
	FD<Dynamic,T> log(FD<Dynamic,T>&& f){
		using std::log;
		FD<Dynamic,T> result(log(f.v),std::move(f));
		for(unsigned int i=0; i<result.nVars; i++)
			result.g[i] /= f.v;
		return(result);
	}

	template <unsigned int nVars, typename T>
	FD<nVars,T> log10(const FD<nVars,T>& f){
		using std::log;
		using std::log10;
		FD<nVars,T> result(log10(f.v),FD<nVars,T>::noInit());
		const T d=log(T(10))*f.v;
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> log10(const FD<Dynamic,T>& f){
		using std::log;
		using std::log10;
		FD<Dynamic,T> result(log10(f.v),f.nVars,FD<Dynamic,T>::noInit());
		const T d=log(T(10))*f.v;
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
    template <unsigned int nVars, typename T>
	FD<nVars,T> log2(const FD<nVars,T>& f){
		using std::log;
		using std::log2;
		FD<nVars,T> result(log2(f.v),FD<nVars,T>::noInit());
		const T d=log(T(2))*f.v;
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> log2(const FD<Dynamic,T>& f){
		using std::log;
		using std::log2;
		FD<Dynamic,T> result(log2(f.v),f.nVars,FD<Dynamic,T>::noInit());
		const T d=log(T(2))*f.v;
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
	//TODO: update all pow functions with noInit and Dynamic specialization
	template <typename FD_t, typename U>
	FD_t pow(const FD_t& b, const U& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type = 0){
		using std::pow;
		typedef typename FD_t::BaseType T;
		FD_t result(b);
		T te(e);
		result.v=pow(b.v,te);
		const unsigned int n=detail::dimensionExtractor<FD,FD_t::N,typename FD_t::BaseType>::nVars(result);
		std::transform(result.g,result.g+n,result.g,
					   std::bind2nd(std::multiplies<T>(),pow(b.v,te-T(1))*e));
		return(result);
	}
	
	template <typename FD_t, typename U>
	FD_t pow(const U& b, const FD_t& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type = 0){
		using std::pow;
		using std::log;
		typedef typename FD_t::BaseType T;
		FD_t result(e);
		T tb(b);
		result.v=pow(tb,e.v);
		const unsigned int n=detail::dimensionExtractor<FD,FD_t::N,typename FD_t::BaseType>::nVars(result);
		std::transform(result.g,result.g+n,result.g,
					   std::bind2nd(std::multiplies<T>(),result.v*log(tb)));
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> pow(const FD<nVars,T>& b, const FD<nVars,T>& e){
		using std::pow;
		using std::log;
		FD<nVars,T> result(pow(b.v,e.v));
		T c1(e.v*pow(b.v,e.v-T(1)));
		T c2(result.v*log(b.v));
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		for(unsigned int i=0; i<n; i++)
			result.g[i] = c1*b.g[i] + c2*e.g[i];
		return(result);
	}
	template <typename T>
	FD<Dynamic,T> pow(const FD<Dynamic,T>& b, const FD<Dynamic,T>& e){
		using std::pow;
		using std::log;
		FD<Dynamic,T> result(pow(b.v,e.v));
		result.changeGradientSize(std::max(detail::dimensionExtractor<FD,Dynamic,T>::nVars(b),
		                                   detail::dimensionExtractor<FD,Dynamic,T>::nVars(e)));
		T c1(e.v*pow(b.v,e.v-T(1)));
		T c2(result.v*log(b.v));
		//TODO: this doesn't properly respect mismatched input gradient sizes
		for(unsigned int i=0; i<result.nVars; i++)
			result.g[i] = c1*b.g[i] + c2*e.g[i];
		return(result);
	}
	template <typename T>
	FD<Dynamic,T> pow(const FD<Dynamic,T>& b, FD<Dynamic,T>&& e){
		assert(b.nVars==e.nVars); //TODO: deal with mismatched sizes
		using std::pow;
		using std::log;
		//FD<Dynamic,T> result(std::move(e));
		//result.v=pow(b.v,e.v);
		FD<Dynamic,T> result(pow(b.v,e.v),b.nVars,FD<Dynamic,T>::noInit());
		//result.changeGradientSize(std::max(detail::dimensionExtractor<FD,Dynamic,T>::nVars(b),
		//                                   detail::dimensionExtractor<FD,Dynamic,T>::nVars(e)));
		T c1(e.v*pow(b.v,e.v-T(1)));
		T c2(result.v*log(b.v));
		//TODO: this doesn't properly respect mismatched input gradient sizes
		for(unsigned int i=0; i<result.nVars; i++)
			//result.g[i] = c1*b.g[i] + c2*result.g[i];
			result.g[i] = c1*b.g[i] + c2*e.g[i];
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> sqrt(const FD<nVars,T>& f){
		using std::sqrt;
		FD<nVars,T> result(sqrt(f.v),FD<nVars,T>::noInit());
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		T d=2*result.v;
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> sqrt(const FD<Dynamic,T>& f){
		using std::sqrt;
		FD<Dynamic,T> result(sqrt(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T d=2*result.v;
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> cbrt(const FD<nVars,T>& f){
		using std::cbrt;
		FD<nVars,T> result(cbrt(f.v),FD<nVars,T>::noInit());
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		T d=3*result.v*result.v;
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> cbrt(const FD<Dynamic,T>& f){
		using std::cbrt;
		FD<Dynamic,T> result(cbrt(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T d=3*result.v*result.v;
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi/d); });
		return(result);
	}
	
	//TODO: hypot

	//TODO: check all atan2 functions for problems with swapped variable names
	template <typename FD_t, typename U>
	FD_t atan2(const FD_t& x, const U& y, typename boost::enable_if< boost::is_arithmetic< U >, int >::type = 0){
		using std::atan2;
		if(x==0){
			if(y>0)
				return(FD_t(boost::math::constants::half_pi<typename FD_t::BaseType>()));
			if(y<0)
				return(FD_t(-boost::math::constants::half_pi<typename FD_t::BaseType>()));
			if(std::numeric_limits<typename FD_t::BaseType>::has_quiet_NaN)
				return(FD_t(std::numeric_limits<typename FD_t::BaseType>::quiet_NaN()));
			if(std::numeric_limits<typename FD_t::BaseType>::has_signaling_NaN)
				return(FD_t(std::numeric_limits<typename FD_t::BaseType>::signaling_NaN()));
			throw std::domain_error("x==0 and y==0 in call to atan2()");
		}
		FD_t result=atan2(x,y);
		if(x<0){
			if(y>=0)
				result+=boost::math::constants::pi<typename FD_t::BaseType>();
			else
				result-=boost::math::constants::pi<typename FD_t::BaseType>();
		}
		return(result);
	}

	template <typename FD_t, typename U>
	FD_t atan2(const U& x, const FD_t& y, typename boost::enable_if< boost::is_arithmetic< U >, int >::type = 0){
		using std::atan2;
		if(x==0){
			if(y>0)
				return(FD_t(boost::math::constants::half_pi<typename FD_t::BaseType>()));
			if(y<0)
				return(FD_t(-boost::math::constants::half_pi<typename FD_t::BaseType>()));
			if(std::numeric_limits<typename FD_t::BaseType>::has_quiet_NaN)
				return(FD_t(std::numeric_limits<typename FD_t::BaseType>::quiet_NaN()));
			if(std::numeric_limits<typename FD_t::BaseType>::has_signaling_NaN)
				return(FD_t(std::numeric_limits<typename FD_t::BaseType>::signaling_NaN()));
			throw std::domain_error("x==0 and y==0 in call to atan2()");
		}
		FD_t result=atan2(x,y);
		if(x<0){
			if(y>=0)
				result+=boost::math::constants::pi<typename FD_t::BaseType>();
			else
				result-=boost::math::constants::pi<typename FD_t::BaseType>();
		}
		return(result);
	}

	template <unsigned int nVars, typename T>
	FD<nVars,T> atan2(const FD<nVars,T>& x, const FD<nVars,T>& y){
		using std::atan2;
		if(x==0){
			if(y>0)
				return(FD<nVars,T>(boost::math::constants::half_pi<T>()));
			if(y<0)
				return(FD<nVars,T>(-boost::math::constants::half_pi<T>()));
			if(std::numeric_limits<T>::has_quiet_NaN)
				return(FD<nVars,T>(std::numeric_limits<T>::quiet_NaN()));
			if(std::numeric_limits<T>::has_signaling_NaN)
				return(FD<nVars,T>(std::numeric_limits<T>::signaling_NaN()));
			throw std::domain_error("x==0 and y==0 in call to atan2()");
		}
		FD<nVars,T> result=atan2(x,y);
		if(x<0){
			if(y>=0)
				result+=boost::math::constants::pi<T>();
			else
				result-=boost::math::constants::pi<T>();
		}
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> erf(const FD<nVars,T>& f){
		using std::erf;
		using std::exp;
		FD<nVars,T> result(erf(f.v),FD<nVars,T>::noInit());
		T m=2*exp(-f.v*f.v)/boost::math::constants::root_pi<T>();
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> erf(const FD<Dynamic,T>& f){
		using std::erf;
		using std::exp;
		FD<Dynamic,T> result(erf(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=2*exp(-f.v*f.v)/boost::math::constants::root_pi<T>();
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> erfc(const FD<nVars,T>& f){
		using std::erfc;
		using std::exp;
		FD<nVars,T> result(erfc(f.v),FD<nVars,T>::noInit());
		T m=-2*exp(-f.v*f.v)/boost::math::constants::root_pi<T>();
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> erfc(const FD<Dynamic,T>& f){
		using std::erfc;
		using std::exp;
		FD<Dynamic,T> result(erfc(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=-2*exp(-f.v*f.v)/boost::math::constants::root_pi<T>();
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> tgamma(const FD<nVars,T>& f){
		using boost::math::tgamma;
		using boost::math::digamma;
		FD<nVars,T> result(tgamma(f.v),FD<nVars,T>::noInit());
		T m=result.v*digamma(f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> tgamma(const FD<Dynamic,T>& f){
		using boost::math::tgamma;
		using boost::math::digamma;
		FD<Dynamic,T> result(tgamma(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=result.v*digmma(f.v);
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> lgamma(const FD<nVars,T>& f){
		using boost::math::lgamma;
		using boost::math::digamma;
		FD<nVars,T> result(lgamma(f.v),FD<nVars,T>::noInit());
		T m=digamma(f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		std::transform(f.g,f.g+n,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <typename T>
	FD<Dynamic,T> lgamma(const FD<Dynamic,T>& f){
		using boost::math::lgamma;
		using boost::math::digamma;
		FD<Dynamic,T> result(lgamma(f.v),f.nVars,FD<Dynamic,T>::noInit());
		T m=digmma(f.v);
		std::transform(f.g,f.g+f.nVars,result.g, [&](const T& gi){ return(gi*m); });
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> ceil(const FD<nVars,T>& f){
		using std::ceil;
		FD<nVars,T> result(f);
		result.v=ceil(f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		//TODO: is it true that (ignoring discontinuities) the derivative of ceil is insensitive to variations in all variables?
		std::fill(result.g,result.g+n,T(0));
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> abs(const FD<nVars,T>& f){
		using std::abs;
		FD<nVars,T> result(f);
		bool flip=f.v<0;
		if(flip){
			result.v=-result.v;
			const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
			std::transform(result.g,result.g+n,result.g,
						   std::bind2nd(std::multiplies<T>(),T(-1)));
		}
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> fabs(const FD<nVars,T>& f){
		using std::fabs;
		FD<nVars,T> result(f);
		bool flip=f.v<0;
		if(flip){
			result.v=-result.v;
			const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
			std::transform(result.g,result.g+n,result.g,
						   std::bind2nd(std::multiplies<T>(),T(-1)));
		}
		return(result);
	}
	
	template <unsigned int nVars, typename T>
	FD<nVars,T> floor(const FD<nVars,T>& f){
		using std::floor;
		FD<nVars,T> result(f);
		result.v=floor(f.v);
		const unsigned int n=detail::dimensionExtractor<FD,nVars,T>::nVars(result);
		//TODO: is it true that (ignoring discontinuities) the derivative of floor is insensitive to variations in all variables?
		std::fill(result.g,result.g+n,T(0));
		return(result);
	}
	
	//comparison operators:
	template <unsigned int nVars, typename T>
	bool operator!(const FD<nVars,T>& f){
		return(!f.value());
	}
	
	//equality
	template <unsigned int nVars, typename T, typename U>
	bool operator==(const FD<nVars,T>& f, const U& u){
		return(f.value()==u);
	}
	
	template <unsigned int nVars, typename T, typename U>
	bool operator==(const U& u, const FD<nVars,T>& f){
		return(u==f.value());
	}
	
	template <unsigned int nVars, typename T>
	bool operator==(const FD<nVars,T>& f1, const FD<nVars,T>& f2){
		return(f1.value()==f2.value());
	}
	
	//inequality
	template <unsigned int nVars, typename T, typename U>
	bool operator!=(const FD<nVars,T>& f, const U& u){
		return(f.value()!=u);
	}
	
	template <unsigned int nVars, typename T, typename U>
	bool operator!=(const U& u, const FD<nVars,T>& f){
		return(u!=f.value());
	}
	
	template <unsigned int nVars, typename T>
	bool operator!=(const FD<nVars,T>& f1, const FD<nVars,T>& f2){
		return(f1.value()!=f2.value());
	}
	
	//greater-than
	template <unsigned int nVars, typename T, typename U>
	bool operator>(const FD<nVars,T>& f, const U& u){
		return(f.value()>u);
	}
	
	template <unsigned int nVars, typename T, typename U>
	bool operator>(const U& u, const FD<nVars,T>& f){
		return(u>f.value());
	}
	
	template <unsigned int nVars, typename T>
	bool operator>(const FD<nVars,T>& f1, const FD<nVars,T>& f2){
		return(f1.value()>f2.value());
	}
	
	//greater-than-or-equal
	template <unsigned int nVars, typename T, typename U>
	bool operator>=(const FD<nVars,T>& f, const U& u){
		return(f.value()>=u);
	}
	
	template <unsigned int nVars, typename T, typename U>
	bool operator>=(const U& u, const FD<nVars,T>& f){
		return(u>=f.value());
	}
	
	template <unsigned int nVars, typename T>
	bool operator>=(const FD<nVars,T>& f1, const FD<nVars,T>& f2){
		return(f1.value()>=f2.value());
	}
	
	//less-than
	template <unsigned int nVars, typename T, typename U>
	bool operator<(const FD<nVars,T>& f, const U& u){
		return(f.value()<u);
	}
	
	template <unsigned int nVars, typename T, typename U>
	bool operator<(const U& u, const FD<nVars,T>& f){
		return(u<f.value());
	}
	
	template <unsigned int nVars, typename T>
	bool operator<(const FD<nVars,T>& f1, const FD<nVars,T>& f2){
		return(f1.value()<f2.value());
	}
	
	//less-than-or-equal
	template <unsigned int nVars, typename T, typename U>
	bool operator<=(const FD<nVars,T>& f, const U& u){
		return(f.value()<=u);
	}
	
	template <unsigned int nVars, typename T, typename U>
	bool operator<=(const U& u, const FD<nVars,T>& f){
		return(u<=f.value());
	}
	
	template <unsigned int nVars, typename T>
	bool operator<=(const FD<nVars,T>& f1, const FD<nVars,T>& f2){
		return(f1.value()<=f2.value());
	}
	
	//stream operators
	template <unsigned int nVars, typename T>
	std::ostream& operator<<(std::ostream& os, const FD<nVars,T>& f){
		os << f.value() << " [";
		for(unsigned int i=0; i<nVars; i++)
			os << (i?",":"") << f.derivative(i);
		os << ']';
		return(os);
	}
	
	template <unsigned int nVars, typename T>
	std::istream& operator>>(std::istream& is, FD<nVars,T>& f){
		T v;
		is >> v;
		f = v;
		return(is);
	}
	
	//TODO: could this be made allocator aware?
	///\brief A type for forward-mode automatic differentiation with a number of variables known only at runtime
	template <typename T>
	class FD<Dynamic,T>{
	private:
		T v; //value
		T* g; //gradient
		unsigned int nVars;
		
		inline void changeGradientSize(unsigned int newVars){
			if(newVars<=nVars)
				return;
			T* newg = new T[newVars];
			//copy over any data we already had
			std::copy(g,g+nVars,newg);
			//zero fill the new portion of the gradient
			std::fill(newg+nVars,newg+newVars,T(0));
			delete[] g;
			g=newg;
			nVars=newVars;
		}
	
		struct noInit_tag{};
		static constexpr noInit_tag noInit(){ return(noInit_tag{}); };
		
		FD<Dynamic,T>(T t, noInit_tag):
		v(t),g(nullptr),nVars(0){}
		
		FD<Dynamic,T>(T t, unsigned int size, noInit_tag):
		v(t),g(nullptr),nVars(size){
			if(nVars)
				g=new T[nVars];
		}
		
		FD<Dynamic,T>(const FD<Dynamic,T>& f, noInit_tag):
		v(f.v),g(nullptr),nVars(f.nVars){
			if(nVars)
				g=new T[nVars];
		}
		
	public:
		typedef T BaseType;
		enum{N=0};
		
		FD<Dynamic,T>():g(nullptr),nVars(0){}
		
		FD<Dynamic,T>(T t):
		v(t),g(nullptr),nVars(0){}
		
		FD<Dynamic,T>(T t, unsigned int index):
		v(t),g(nullptr),nVars(0){
			changeGradientSize(index+1);
			g[index]=T(1);
		}
		
		FD<Dynamic,T>(T t, unsigned int index, unsigned int size):
		v(t),g(nullptr),nVars(0){
			assert(size>index);
			changeGradientSize(size);
			g[index]=T(1);
		}
		
		FD<Dynamic,T>(const FD<Dynamic,T>& f):
		v(f.v),g(nullptr),nVars(0){
			changeGradientSize(f.nVars);
			std::copy(f.g,f.g+nVars,g);
		}
		
		FD<Dynamic,T>(FD<Dynamic,T>&& f):
		v(f.v),g(f.g),nVars(f.nVars){
			//report that we've stolen the contents of f
			f.g=nullptr;
			f.nVars=0;
		}
		
		FD<Dynamic,T>(T t, FD<Dynamic,T>&& f):
		v(t),g(f.g),nVars(f.nVars){
			//report that we've stolen the contents of f
			f.g=nullptr;
			f.nVars=0;
		}
		
		FD<Dynamic,T>& operator=(const FD<Dynamic,T>& f){
			if(&f!=this){
				v=f.v;
				changeGradientSize(f.nVars);
				std::copy(f.g,f.g+nVars,g);
			}
			return(*this);
		}
		
		FD<Dynamic,T>& operator=(FD<Dynamic,T>&& f){
			if(&f!=this){
				v=f.v;
				nVars=f.nVars;
				g=f.g;
				//report that we've stolen the contents of f
				f.g=nullptr;
				f.nVars=0;
			}
			return(*this);
		}
		
		FD<Dynamic,T>& operator=(const T& t){
			v=t;
			std::fill(g,g+nVars,0);
			return(*this);
		}
		
		~FD<Dynamic,T>(){
			delete[] g;
		}
		
		void assignIndex(unsigned int index){
			changeGradientSize(index+1);
			assert(index<nVars);
			g[index]=T(1);
		}
		
		const T& value() const{
			return(v);
		}
		
		const T& derivative(unsigned int index) const{
			assert(index<nVars);
			return(g[index]);
		}
		
		void copyGradient(T grad[]) const{
			std::copy(g,g+nVars,grad);
		}
		
		///Manually set a component of the derivative to a particular value.
		///This is mostly only useful when interfacing with other code which
		///somehow computes gradient informaton manually and one wishes to then
		///propagate it further using automatic differentiation.
		void setDerivative(unsigned int index, T d){
			assert(index<nVars);
			g[index]=d;
		}
		
		//unary +
		FD<Dynamic,T> operator+() const{
			return(*this);
		}
		
		//unary -
		FD<Dynamic,T> operator-() const{
			FD<Dynamic,T> r(*this);
			r.v=-r.v;
			r.changeGradientSize(nVars);
			std::transform(r.g,r.g+nVars,r.g,std::negate<T>());
			return(r);
		}
		
		//addition
		template <typename U>
		FD<Dynamic,T>& operator+=(const U& u){
			v+=T(u);
			return(*this);
		}
		
		FD<Dynamic,T>& operator+=(const FD<Dynamic,T>& f){
			v+=f.v;
			changeGradientSize(f.nVars);
			std::transform(g,g+f.nVars,f.g,g,std::plus<T>());
			return(*this);
		}
		
		template <typename U>
		FD<Dynamic,T> operator+(const U& u) const &{
			FD<Dynamic,T> result(v+u,nVars,noInit());
			std::copy(g,g+nVars,result.g);
			return(result);
		}
		
		template <typename U>
		FD<Dynamic,T> operator+(const U& u) &&{
			return(FD<Dynamic,T>(std::move(*this))+=u);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator+(const FD<Dynamic,T>& f) const &{
			FD<Dynamic,T> result(v+f.v,std::max(nVars,f.nVars),noInit());
			unsigned int n=std::min(nVars,f.nVars);
			for(unsigned int i=0; i<n; i++)
				result.g[i] = g[i] + f.g[i];
			for(unsigned int i=n; i<nVars; i++)
				result.g[i] = g[i];
			for(unsigned int i=n; i<f.nVars; i++)
				result.g[i] = f.g[i];
			return(result);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator+(FD<Dynamic,T>&& f) const &{
			if(f.nVars>nVars){
				FD<Dynamic,T> result(std::move(f));
				result.v+=v;
				for(unsigned int i=0; i<nVars; i++)
					result.g[i] += g[i];
				return(result);
			}
			FD<Dynamic,T> result(v+f.v,nVars,noInit());
			for(unsigned int i=0; i<f.nVars; i++)
				result.g[i] = g[i] + f.g[i];
			for(unsigned int i=f.nVars; i<nVars; i++)
				result.g[i] = g[i];
			return(result);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator+(const FD<Dynamic,T>& f) &&{
			if(nVars>f.nVars){
				FD<Dynamic,T> result(std::move(*this));
				result.v+=f.v;
				for(unsigned int i=0; i<f.nVars; i++)
					result.g[i] += f.g[i];
				return(result);
			}
			FD<Dynamic,T> result(v+f.v,f.nVars,noInit());
			for(unsigned int i=0; i<nVars; i++)
				result.g[i] = g[i] + f.g[i];
			for(unsigned int i=nVars; i<f.nVars; i++)
				result.g[i] = f.g[i];
			return(result);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator+(FD<Dynamic,T>&& f) &&{
			if(nVars>f.nVars){
				FD<Dynamic,T> result(std::move(*this));
				result.v+=f.v;
				for(unsigned int i=0; i<f.nVars; i++)
					result.g[i] += f.g[i];
				return(result);
			}
			FD<Dynamic,T> result(std::move(f));
			result.v+=v;
			for(unsigned int i=0; i<nVars; i++)
				result.g[i] += g[i];
			return(result);
		}
		
		//subtraction
		template <typename U>
		FD<Dynamic,T>& operator-=(const U& u){
			v-=T(u);
			return(*this);
		}
		
		FD<Dynamic,T>& operator-=(const FD<Dynamic,T>& f){
			v-=f.v;
			changeGradientSize(f.nVars);
			std::transform(g,g+f.nVars,f.g,g,std::minus<T>());
			return(*this);
		}
		
		template <typename U>
		FD<Dynamic,T> operator-(const U& u) const{
			FD<Dynamic,T> result(v-u,nVars,noInit());
			std::copy(g,g+nVars,result.g);
			return(result);
		}
		
		//TODO: replace this with optimized implementations as for addition
		FD<Dynamic,T> operator-(const FD<Dynamic,T>& f) const{
			return(FD<Dynamic,T>(*this)-=f);
		}
		
		//multiplication
		template <typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<Dynamic,T>& operator*=(const U& u){
			T t(u);
			v*=t;
			std::transform(g,g+nVars,g,[&](const T& gi){ return(gi*t); });
			return(*this);
		}
		
		FD<Dynamic,T>& operator*=(const FD<Dynamic,T>& f){
			changeGradientSize(f.nVars);
			for(unsigned int i=0; i<f.nVars; i++)
				g[i] = g[i]*f.v + f.g[i]*v;
			v*=f.v;
			return(*this);
		}
		
		template<typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		FD<Dynamic,T> operator*(const U& u) const &{
			FD<Dynamic,T> result(v*u,nVars,noInit());
			for(unsigned int i=0; i<nVars; i++)
				result.g[i]=u*g[i];
			return(result);
		}
		
		template<typename U, typename=typename std::enable_if<std::is_constructible<T,U>::value>::type>
		__attribute__((always_inline)) FD<Dynamic,T> operator*(const U& u) &&{
			FD<Dynamic,T> result(std::move(*this));
			result.v*=u;
			for(unsigned int i=0; i<result.nVars; i++)
				result.g[i]*=u;
			return(result);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator*(const FD<Dynamic,T>& f) const &{
			FD<Dynamic,T> result(v*f.v,std::max(nVars,f.nVars),noInit());
			unsigned int n=std::min(nVars,f.nVars);
			for(unsigned int i=0; i<n; i++)
				result.g[i] = g[i]*f.v + f.g[i]*v;
			for(unsigned int i=n; i<nVars; i++)
				result.g[i] = g[i]*f.v;
			for(unsigned int i=n; i<f.nVars; i++)
				result.g[i] = f.g[i]*v;
			return(result);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator*(FD<Dynamic,T>&& f) const &{
			if(f.nVars>nVars){
				FD<Dynamic,T> result(std::move(f));
				result.v*=v;
				for(unsigned int i=0; i<nVars; i++)
					result.g[i] = g[i]*f.v + result.g[i]*v;
				for(unsigned int i=nVars; i<result.nVars; i++)
					result.g[i] *= v;
				return(result);
			}
			FD<Dynamic,T> result(v*f.v,nVars,noInit());
			for(unsigned int i=0; i<f.nVars; i++)
				result.g[i] = g[i]*f.v + f.g[i]*v;
			for(unsigned int i=f.nVars; i<nVars; i++)
				result.g[i] = g[i]*f.v;
			return(result);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator*(const FD<Dynamic,T>& f) &&{
			if(nVars>f.nVars){
				FD<Dynamic,T> result(std::move(*this));
				result.v*=f.v;
				for(unsigned int i=0; i<f.nVars; i++)
					result.g[i] = result.g[i]*f.v + f.g[i]*v;
				for(unsigned int i=f.nVars; i<result.nVars; i++)
					result.g[i] *= f.v;
				return(result);
			}
			FD<Dynamic,T> result(v*f.v,f.nVars,noInit());
			for(unsigned int i=0; i<nVars; i++)
				result.g[i] = g[i]*f.v + f.g[i]*v;
			for(unsigned int i=nVars; i<f.nVars; i++)
				result.g[i] = f.g[i]*v;
			return(result);
		}
		
		__attribute__((always_inline)) FD<Dynamic,T> operator*(FD<Dynamic,T>&& f) &&{
			if(nVars>f.nVars){
				FD<Dynamic,T> result(std::move(*this));
				result.v*=f.v;
				for(unsigned int i=0; i<f.nVars; i++)
					result.g[i] = result.g[i]*f.v + f.g[i]*v;
				for(unsigned int i=f.nVars; i<result.nVars; i++)
					result.g[i] *= f.v;
				return(result);
			}
			FD<Dynamic,T> result(std::move(f));
			result.v*=v;
			for(unsigned int i=0; i<nVars; i++)
				result.g[i] = g[i]*f.v + result.g[i]*v;
			for(unsigned int i=nVars; i<result.nVars; i++)
				result.g[i] *= v;
			return(result);
		}
		
		//division
		template<typename U>
		FD<Dynamic,T>& operator/=(const U& u){
			T t(u);
			v/=t;
			std::transform(g,g+nVars,g,[&](const T& gi){ return(gi/t); });
			return(*this);
		}
		
		FD<Dynamic,T>& operator/=(const FD<Dynamic,T>& f){
			v/=f.v;
			changeGradientSize(f.nVars);
			for(unsigned int i=0; i<f.nVars; i++)
				g[i] = (g[i] - v*f.g[i])/f.v;
			return(*this);
		}
		
		template<typename U>
		FD<Dynamic,T> operator/(const U& u) const &{
			FD<Dynamic,T> result(v/u,nVars,noInit());
			for(unsigned int i=0; i<nVars; i++)
				result.g[i]=g[i]/u;
			return(result);
		}
		
		template<typename U>
		__attribute__((always_inline)) FD<Dynamic,T> operator/(const U& u) &&{
			return(FD<Dynamic,T>(std::move(*this))/=u);
			/*FD<Dynamic,T> result(v/u,std::move(*this));
			for(unsigned int i=0; i<result.nVars; i++)
				result.g[i]/=u;
			return(result);*/
		}
		
		//TODO: replace this with optimized implementations as for multiplication
		FD<Dynamic,T> operator/(const FD<Dynamic,T>& f) const{
			return(FD<Dynamic,T>(*this)/=f);
		}
		
		template <template<unsigned int,class> class F, int D, class T_>
		friend struct detail::dimensionExtractor;
		template <unsigned int nVars_, typename T_, typename U>
		friend FD<nVars_,T_> operator/(const U& u, const FD<nVars_,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> cos(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> sin(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> tan(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> acos(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> asin(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> atan(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> cosh(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> sinh(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> tanh(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> acosh(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> asinh(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> atanh(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> exp(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> exp(FD<Dynamic,T_>&& f);
		template <typename T_>
		friend FD<Dynamic,T_> log(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> log(FD<Dynamic,T_>&& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> log10(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> log2(const FD<nVars_,T_>& f);
		template <typename FD_t, typename U>
		friend FD_t pow(const FD_t& b, const U& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type);
		template <typename FD_t, typename U>
		friend FD_t pow(const U& b, const FD_t& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type);
		template <typename FD_t, typename U>
		friend FD_t pow(const U& b, const FD_t& e, typename boost::enable_if< boost::is_arithmetic< U >, int >::type);
		template <typename T_>
		friend FD<Dynamic,T_> pow(const FD<Dynamic,T_>& b, const FD<Dynamic,T_>& e);
		template <typename T_>
		friend FD<Dynamic,T_> pow(const FD<Dynamic,T_>& b, FD<Dynamic,T_>&& e);
		template <typename T_>
		friend FD<Dynamic,T_> sqrt(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> cbrt(const FD<Dynamic,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> floor(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> ceil(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> abs(const FD<nVars_,T_>& f);
		template <unsigned int nVars_, typename T_>
		friend FD<nVars_,T_> fabs(const FD<nVars_,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> tgamma(const FD<Dynamic,T_>& f);
		template <typename T_>
		friend FD<Dynamic,T_> lgamma(const FD<Dynamic,T_>& f);
	};
	
} //namespace autodiff
} //namespace phys_tools
	
namespace std{
	//an autodiff::FD<nVars,T> has all the same properties as a T
	template <unsigned int nVars, typename T>
	class numeric_limits<phys_tools::autodiff::FD<nVars,T> > : public numeric_limits<T>{};
	
	template <unsigned int nVars, typename T>
	bool isnan(const phys_tools::autodiff::FD<nVars,T>& f){
        bool gradnan = false;
		for(unsigned int i=0; i<nVars; i++)
			gradnan |= isnan(f.derivative(i));
		return(isnan(f.value()) || gradnan);
	}
}
	
#endif //AUTODIFF_H_INCLUDED
