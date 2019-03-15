#ifndef LF_LIKELIHOOD_H
#define LF_LIKELIHOOD_H

#define BOOST_RESULT_OF_USE_DECLTYPE

#include <cmath>
#include <functional>
#include <numeric>
#include <random>

#ifdef __APPLE__
    #include <xmmintrin.h>
#elif __linux__
    #include <stdexcept>
#endif

#include <type_traits>
#include <vector>
#include <deque>
#include <cfenv>

#include <boost/math/constants/constants.hpp>
#include <boost/mpl/if.hpp>

#include "../lbfgsb/lbfgsb.h"

#include "../histogram.h"
#include "../brent.h"
#include "../residuals.h"

#include "../autodiff.h"
#include "../ThreadPool.h"

namespace phys_tools{
///Tools for performing binned maximum likelihood fits
namespace likelihood{

	//some useful type traits
	template<typename T>
	struct remove_reference_wrapper{ using type=T; };

	template<typename T>
	struct remove_reference_wrapper<std::reference_wrapper<T>>{ using type=T; };

	template<typename T>
	struct ensure_reference_wrapper{ using type=typename std::reference_wrapper<T>; };

	template<typename T>
	struct ensure_reference_wrapper<std::reference_wrapper<T>>{
		using type=typename std::reference_wrapper<T>;
	};

    template <class InIt>
    typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
        typedef typename std::iterator_traits<InIt>::value_type real;
        real sum = real(0);
        real running_error = real(0);
        real temp;
        real difference;

        for (; begin != end; ++begin) {
            difference = *begin;
            difference -= running_error;
            temp = sum;
            temp += difference;
            running_error = temp;
            running_error -= sum;
            running_error -= difference;
            sum = std::move(temp);
        }
        return sum;
    };

	//Frequently used likelihood functions

	struct poissonLikelihood{
		template <typename T>
		T operator()(double dataCount, T const & lambda, T const & w2_sum) const{
			if(lambda==0)
                return(dataCount==0?0:-std::numeric_limits<T>::max());
			//return(0);
			//would this be more correct?
			T sum(lambda);
			sum+=lgamma(dataCount+1);
			return(dataCount*log(lambda)-sum);
		}
	};

	struct dimaLikelihood{
		template<typename T>
		T operator()(double dataCount, const std::vector<T>& simulationWeights, int n_events) const{
			//std::cout << "   " << dataCount << " observed events compared to " << simulationWeights.size() << " simulated events" << '\n';
			using std::abs; //we want access to std::abs for things like doubles so that we can also find things under the same name via ADL

			if(simulationWeights.empty())
				return(T(0));
			//to save on divides later, transform the weights into times: t_ij = 1/w_ij
			std::vector<T> t_i;
			t_i.reserve(simulationWeights.size());
			for(const auto& w_ij : simulationWeights)
				if(w_ij!=0) //drop events with weight zero; they should be irrelevant
					t_i.emplace_back(1./w_ij);

			T lambda_i(0), lambda_i_last(1);
			T ssum(0);
			//first we need to compute this bin's individual lagrange multiplier, lamdba_i, using the
			//Newton-Raphson method. However, the constraint function has poles, which the derivative
			//extrapolations may cause the N-R method to cross, leading to cycling or divergence. So, we
			//also compute the bounds on the allowed range for lambda_i, and if the extrapolation would take
			//us outside those bounds we take a bisection step toward the indicated bound.
			T lambda_i_max(1);
			T lambda_i_min(-std::numeric_limits<T>::infinity());
//			for(auto w_ij : simulationWeights){
//				if(-1/w_ij>lambda_i_min)
//					lambda_i_min=-1/w_ij;
//				ssum+=1/(1/w_ij+lambda_i); //while here also prepare to compute initial value of R_i
//			}
			for(auto t_ij : t_i){
				assert(t_ij>0.0);
				if(-t_ij>lambda_i_min){
					lambda_i_min=-t_ij;
					//std::cout << "   found t_ij=" << t_ij << ", changed lambda_i_min to " << lambda_i_min << std::endl;
				}
				ssum+=1/(t_ij+lambda_i); //while here also prepare to compute initial value of R_i
			}
			//std::cout << "   lambda_i domain: [" << lambda_i_min << ',' << lambda_i_max << ']' << std::endl;
			T R_i_last=dataCount/(1-lambda_i)-ssum;
			if(dataCount==0){
				lambda_i=1;
				R_i_last=0; //skip the next loop
			}
			//std::cout << "   initial lambda_i=" << lambda_i << std::endl;

			//TODO: should convergence criteria be more rigorous?
			const double /*change_tol=1e-5,*/ root_tol=1e-10;
			unsigned int n=0;
			while(/*abs(lambda_i-lambda_i_last)>change_tol*abs(lambda_i_last) &&*/ abs(R_i_last)>root_tol){
				lambda_i_last=lambda_i;

				ssum=0;
				T ssum2(0);
//				for(auto w_ij : simulationWeights){
//					T t=1/(1/w_ij+lambda_i);
//					ssum+=t;
//					ssum2+=t*t;
//				}
				for(auto t_ij : t_i){
					T t=1/(t_ij+lambda_i);
					ssum+=t;
					ssum2+=t*t;
				}

				T z=1-lambda_i;
				if(z==0) //don't even try the N-R step
					lambda_i=(lambda_i_max+lambda_i_last)/2;
				else{
					//Try the N-R step
					lambda_i-=(dataCount/z-ssum)/(dataCount/(z*z)+ssum2);
					//std::cout << "    N-R stepped to " << std::setprecision(16) << lambda_i << std::setprecision(6) << std::endl;
					if(lambda_i>lambda_i_max)
						lambda_i=(lambda_i_max+lambda_i_last)/2;
					else if(lambda_i<lambda_i_min)
						lambda_i=(lambda_i_min+lambda_i_last)/2;
				}
				R_i_last=dataCount/(1-lambda_i_last)-ssum;

				//ssum=0;
				////for(auto w_ij : simulationWeights)
				////	ssum+=1/(1/w_ij+lambda_i);
				//for(auto t_ij : t_i)
				//	ssum+=1/(t_ij+lambda_i);
				//T R_i=(dataCount/(1-lambda_i)-ssum);

				//std::cout << "    labmda_i now " << std::setprecision(16) << lambda_i << std::setprecision(6) << std::endl;
				//std::cout << "     R_i_last = " << R_i_last << std::endl;
				//assert(!std::isinf((double)R_i_last));
				//std::cout << "     R_i = " << R_i << std::endl;
				if(n++>100){
					assert(abs(R_i_last)<10*root_tol); //root must not be too poor
					break;
				}
			}
			//std::cout << "   lambda_i=" << lambda_i << '\n';

			T llh(0);
//			for(auto w_ij : simulationWeights){
//				//log(t_ij)-log(t_ij+lambda_i) = -log(w_ij)-log(1/w_ij+lambda_i)
//				llh-=log(w_ij)+log(1/w_ij+lambda_i);
//				//std::cout << "    w: " << w_ij << ' ' << llh << std::endl;
//			}
			for(auto t_ij : t_i)
				llh+=log(t_ij)-log(t_ij+lambda_i);
			if(dataCount)
				llh+=dataCount*-log(1-lambda_i);
			//std::cout << "    dataCount=" << dataCount << " log(1-lambda_i)=" << log(1-lambda_i) << '\n';
			//std::cout << "   llh=" << llh << '\n';
			return(llh);
		}
	};

	struct chi2Likelihood{
		template<typename T>
		T operator()(unsigned int dataCount, const std::vector<T>& expectationWeights, int n_events) const{
			T exp=std::accumulate(expectationWeights.begin(),expectationWeights.end(),T(0),std::plus<T>());
			if(exp>0){
				T diff=dataCount-exp;
				return((diff*diff)/exp);
			}
			return(0.0);
		}
	};

	struct saturatedPoissonLikelihood{
		template <typename T>
		T operator()(double dataCount, const std::vector<T>& simulationWeights, int n_events) const{
			if(dataCount==0)
				return(0);
			T sum(dataCount);
			sum+=lgamma(dataCount+1);
			return(dataCount*log(dataCount)-sum);
		}
	};

	///From Bohm and Zech, "Comparison of experimental data to Monte Carlo
	///simulation—Parameter estimation and goodness-of-fit testing with
	///weighted events", 2012
	///http://dx.doi.org/10.1016/j.nima.2012.06.021
	struct bohmZechLikelihood{
		template<typename T>
		T operator()(double dataCount, const std::vector<T>& simulationWeights) const{
			//having an observation where there is no simulation is a problem, but
			//we can't really solve it in his context, so we return 0 to leave
			//the likelihood unaffected
			if(simulationWeights.size()==0)
				return(0);
			double n=dataCount;
			T m(0), wk2(0);
			for(const T& wk : simulationWeights){
				m+=wk;
				wk2+=wk*wk;
			}
			//this is the effective number of simulated events
			T mt(m*m/wk2); //\tilde{m}
			//this is the effective scaling factor for the simulation,
			//the average weight of each effective event
			T cmt(m/mt); //\tilde{c_m}
			T lambdat=((n+mt)/(1+1/cmt)); //\tilde{\lambda}
			T slambdat=lambdat/cmt; //\tilde{\lambda} \over \tilde{c_m}

			T sum(lambdat+slambdat);
			sum+=lgamma(n+1);
			sum+=lgamma(mt+1);
			return(n*log(lambdat)+mt*log(slambdat)-sum);
		}
	};

    template<typename DataType>
    DataType compute_barlow_ti(double di, const std::vector<unsigned int>& ai, const std::vector<DataType>& wi) {
        using T=DataType;
        auto max_wi_it = std::max_element(wi.begin(), wi.end());
        unsigned int max_wi_i = std::distance(wi.begin(), max_wi_it);
        T lower_bound = -T(1.0) / T(*max_wi_it);
        T upper_bound(1.0);
        T tol(1e-10);

        std::function<T(T)> func = [&](T ti)->T{
            T result = T(di) / (T(1) - ti);
            std::vector<T> results(ai.size(), T(0));
            for(unsigned int j=0; j<ai.size(); ++j) {
                results[j] = (wi[j]*ai[j]) / (T(1)+wi[j]*ti);
            }
            return result - std::accumulate(results.begin(), results.end(), T(0), std::plus<T>());
        };

        T ti(brent::zero(lower_bound, upper_bound, tol, func));
        return ti;
    };

    template<typename DataType>
    struct compute_barlow_LLH {
        using T=DataType;
        DataType operator()(double di, const std::vector<unsigned int>& ai, const std::vector<DataType>& wi) {
            std::vector<T> Ai(ai.size());
            T ti(1);
            if(di > 0) {
                std::function<T(T)> func = [&](T ti)->T{
                    T result = T(di) / (T(1) - ti);
                    std::vector<T> results(ai.size(), T(0));
                    for(unsigned int j=0; j<ai.size(); ++j) {
                        results[j] = (wi[j]*ai[j]) / (T(1)+wi[j]*ti);
                    }
                    return result - std::accumulate(results.begin(), results.end(), T(0), std::plus<T>());
                };

                auto max_wi_it = std::max_element(wi.begin(), wi.end());
                unsigned int max_wi_i = std::distance(wi.begin(), max_wi_it);
                T lower_bound = -T(1.0) / T(*max_wi_it);
                T upper_bound(1.0);
                T tol(1e-10);


                ti = compute_barlow_ti(di, ai, wi);

                // Special case: MC source with the largest strength has zero MC events
                if(ai[max_wi_i] == 0) {
                    T Aki = T(di) / (T(1.0) + T(*max_wi_it));
                    for(unsigned int j=0; j<ai.size(); ++j) {
                        if(j == max_wi_i) {
                            continue;
                        }
                        if(*max_wi_it == wi[j]) {
                            Aki = T(0);
                            break;
                        }
                        Aki -= wi[j]*ai[j]/(*max_wi_it - wi[j]);
                    }
                    if(Aki < 0)
                        Aki = T(0);
                    Ai[max_wi_i] = Aki;
                    if(Aki > T(0)) {
                        ti = lower_bound;
                        for(unsigned int j=0; j<ai.size(); ++j) {
                            if(j == max_wi_i)
                                continue;
                            Ai[j] = ai[j] / (T(1)+wi[j]*ti);
                        }
                    }
                }
                else {
                    for(unsigned int j=0; j<ai.size(); ++j) {
                        Ai[j] = ai[j] / (T(1)+wi[j]*ti);
                    }
                }
            }
            else {
                //std::cout << "no data" << std::endl;
                for(unsigned int j=0; j<ai.size(); ++j) {
                    Ai[j] = ai[j] / (T(1)+wi[j]*ti);
                }
            }

            T f(0);
            T sum(0);

            for(unsigned int j=0; j<ai.size(); ++j) {
                if(Ai[j] < T(0))
                    Ai[j] = T(0);
                f += wi[j]*Ai[j];
                sum += ai[j]*log(Ai[j]) - Ai[j] - lgamma(ai[j]+T(1));
            }
            sum += T(di)*log(f) - f - lgamma(T(di)+T(1));

            return sum;
        }
    };

    template<unsigned int Dim, typename T>
    struct compute_barlow_LLH<phys_tools::autodiff::FD<Dim, T>> {
        using result_type=phys_tools::autodiff::FD<Dim, T>;
        result_type operator()(double di, const std::vector<unsigned int>& ai, const std::vector<result_type>& dwi) {
            std::vector<T> wi(dwi.size());
            for(unsigned int j; j<dwi.size(); ++j) {
                wi[j] = dwi[j].value();
            }

            bool special_case = false;
            auto max_wi_it = std::max_element(wi.begin(), wi.end());
            unsigned int max_wi_i = std::distance(wi.begin(), max_wi_it);
            T lower_bound = -T(1.0) / T(*max_wi_it);
            T ti(1);
            result_type dti(1);

            if(di == 0) {
                ti = T(1);
                dti = ti;
                const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,T>::nVars(dti);
                for(unsigned int i; i<n; ++i) {
                    dti.setDerivative(i, T(0));
                }
            }
            else
                ti = compute_barlow_ti<T>(di, ai, wi);

            std::vector<T> Ai(ai.size());

            // Special case: MC source with the largest strength has zero MC events
            if(ai[max_wi_i] == 0) {
                T Aki = T(di) / (T(1.0) + T(*max_wi_it));
                for(unsigned int j=0; j<ai.size(); ++j) {
                    if(j == max_wi_i) {
                        continue;
                    }
                    if(*max_wi_it == wi[j]) {
                        Aki = T(0);
                        break;
                    }
                    Aki -= wi[j]*ai[j]/(*max_wi_it - wi[j]);
                }
                if(Aki < 0)
                    Aki = T(0);
                Ai[max_wi_i] = Aki;
                if(Aki > T(0)) {
                    ti = lower_bound;
                    special_case = true;
                    for(unsigned int j=0; j<ai.size(); ++j) {
                        if(j == max_wi_i)
                            continue;
                        Ai[j] = ai[j] / (T(1)+wi[j]*ti);
                    }
                }
                else {
                    for(unsigned int j=0; j<ai.size(); ++j) {
                        Ai[j] = ai[j] / (T(1)+wi[j]*ti);
                    }
                }
            }
            else {
                for(unsigned int j=0; j<ai.size(); ++j) {
                    Ai[j] = ai[j] / (T(1)+wi[j]*ti);
                }
            }
            dti = ti;

            T fi(0);
            T Li(0);
            for(unsigned int j=0; j<ai.size(); ++j) {
                if(Ai[j] < T(0))
                    Ai[j] = T(0);
                fi += wi[j]*Ai[j];
                Li += ai[j]*log(Ai[j]) - Ai[j] - lgamma(ai[j]+T(1));
            }
            Li += T(di)*log(fi) - fi - lgamma(T(di)+T(1));

            result_type dfi(fi);
            result_type dLi(Li);
            std::vector<result_type> dAi;
            for(unsigned int j=0; j<ai.size(); ++j) {
                dAi.push_back(result_type(Ai[j]));
            }

            const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,T>::nVars(dti);
            for(unsigned int i; i<n; ++i) {
                T dt(0);
                T dt_den(0);
                T dA(0);
                T tmp(0);
                if(di>0) {
                    if(special_case) {
                        dt = T(1)/pow(wi[max_wi_i],T(2)) * dwi[max_wi_i].derivative(i);
                    }
                    else{
                        if(di > T(0)) {
                            for(unsigned int j=0; j<ai.size(); ++j) {
                                tmp = (wi[j]*ai[j])/pow(T(1)+wi[j]*ti, T(2));
                                dt += dwi[j].derivative(i) * (Ai[j] - tmp*ti);
                                dt_den += (pow(fi, T(2))/di + tmp*wi[j]);
                            }
                            dt /= dt_den;
                        }
                    }
                }
                for(unsigned int j=0; j<ai.size(); ++j) {
                    if(!special_case || j!=max_wi_i) {
                        dAi[j].setDerivative(i,-T(ai[j])/pow(T(1)+wi[j]*ti,T(2))*(dwi[j].derivative(i)*ti+dt));
                        dA += dAi[j].derivative(i);
                    }
                }
                if(special_case) {
                    dAi[max_wi_i].setDerivative(i, -di/pow(T(1)+wi[max_wi_i], T(2)) - dA);
                }

                dti.setDerivative(i, dt);
                if(di>T(0)) {
                    dfi.setDerivative(i, dt*pow(fi,T(2))/di);
                }
                else {
                    T df(0);
                    for(unsigned int j=0; j<ai.size(); ++j) {
                        df += dwi[j].derivative(i) * Ai[j] + wi[j] * dAi[j].derivative(i);
                    }
                    dfi.setDerivative(i, df);
                }
                T dL((di/fi-T(1))*dfi.derivative(i));
                for(unsigned int j=0; j<ai.size(); ++j) {
                    assert(Ai[j]>0||ai[j]==0);
                    if(Ai[j] > T(0)) {
                        dL += (ai[j]/Ai[j] - T(1))*dAi[j].derivative(i);
                    }
                }
                dLi.setDerivative(i, dL);
            }
            return dLi;
        }
    };

    struct barlowSimpleLikelihood {
        template<typename T>
        T operator()(double dataCount, const std::vector<T>& simulationWeights, const std::vector<unsigned int>& categories) const {
            auto it = simulationWeights.begin();
            auto end = simulationWeights.begin();
            std::vector<T> wi(1);
            const std::vector<unsigned int> ai(1, simulationWeights.size());
            for(unsigned int j=0; j<wi.size(); ++j) {
                end += simulationWeights.size();
                wi[j] = std::accumulate(it, end, T(0), std::plus<T>()) / ai[j];
                it = end;
            }
            if(ai.size()>0) {
                return compute_barlow_LLH<T>()(dataCount, ai, wi);
            }
            else {
                return T(0);
            }
        }
    };

    struct barlowLikelihood {
        template<typename T>
        T operator()(double dataCount, const std::vector<T>& simulationWeights, const std::vector<unsigned int>& categories) const {
            auto it = simulationWeights.begin();
            auto end = simulationWeights.begin();
            std::vector<T> wi(categories.size());
            const std::vector<unsigned int> &ai = categories;
            for(unsigned int j=0; j<categories.size(); ++j) {
                end += categories[j];
                wi[j] = std::accumulate(it, end, T(0), std::plus<T>()) / ai[j];
                it = end;
            }
            if(ai.size()>0) {
                return compute_barlow_LLH<T>()(dataCount, ai, wi);
            }
            else {
                return T(0);
            }
        }
    };

    struct barlowExtendedLikelihood {
        template<typename T>
        T operator()(double dataCount, const std::vector<T>& simulationWeights, const std::vector<unsigned int>& categories) const {
            auto it = simulationWeights.begin();
            auto end = simulationWeights.begin();
            const std::vector<T> & wi = simulationWeights;
            const std::vector<unsigned int> ai(simulationWeights.size(), 1);
            if(ai.size()>0) {
                return compute_barlow_LLH<T>()(dataCount, ai, wi);
            }
            else {
                return T(0);
            }
        }
    };

    // compute log(1+x) without losing precision for small values of x
    template<typename T>
    T LogOnePlusX(T x)
    {
        if (x <= -1.0)
        {
            std::stringstream os;
            os << "Invalid input argument (" << x
               << "); must be greater than -1.0";
            throw std::invalid_argument( os.str() );
        }

        if (fabs(x) > 1e-4)
        {
            // x is large enough that the obvious evaluation is OK
            return log(1.0 + x);
        }

        // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
        // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
        T x2 = x*x;
        T x3 = x2*x;
        T x4 = x3*x;
        return x-x2/2.0+x3/3.0-x4/4.0;
    };

    template<typename T>
    T LogOneMinusX(T x)
    {
        if (x <= -1.0)
        {
            std::stringstream os;
            os << "Invalid input argument (" << x
               << "); must be greater than -1.0";
            throw std::invalid_argument( os.str() );
        }

        if (fabs(x) > 1e-4)
        {
            // x is large enough that the obvious evaluation is OK
            return log(1.0 - x);
        }

        // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
        // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
        T x2 = x*x;
        T x3 = x2*x;
        T x4 = x3*x;
        return -x-x2/2.0-x3/3.0-x4/4.0;
    };

    template<typename T>
    T SFromW(T w) {
        if(fabs(w) >= 1e-4) {
            return 1.0 / (1.0+1.0/w);
        }

        return w - pow(w, 2.0) + pow(w, 3.0) - pow(w, 4.0) + pow(w, 5.0);
    };

    template<class T, class InputItA, class InputItB>
    T LogSumExp(InputItA a_first, InputItA a_last, InputItB b_first, InputItB b_last) {
        assert(std::distance(a_first, a_last) == std::distance(b_first, b_last));
        T zero(0);
        T minf(-std::numeric_limits<double>::infinity());
        std::vector<T> a_copy(a_first, a_last);
        if(a_copy.size() == 1)
            return T(*a_first + log(*b_first));
        auto a_it = a_copy.begin();
        InputItB b_it = b_first;
        for(; a_it != a_copy.end() && b_it != b_last; ++a_it, ++b_it) {
            if(*b_it == 0)
                *a_it = minf;
        }
        T a_max = *std::max_element(a_copy.begin(), a_copy.end());
        a_it = a_copy.begin();
        b_it = b_first;
        for(; a_it != a_copy.end() && b_it != b_last; ++a_it, ++b_it) {
            *a_it = *b_it * exp(*a_it - a_max);
        }
        return a_max + log(accumulate(a_copy.begin(), a_copy.end()));
    }
    /*
    template<typename T>
    T LogSumExp(const & std::vector<T> a, const & std::vector<T> b) {
        assert(a.size() == b.size());
        T zero(0);
        T minf(-std::numeric_limits<double>::infinity());
        std::vector<T> a_copy(a);
        for(unsigned int i=0; i<a.size(); ++i) {
            if(b[i] == 0)
                a_copy[i] = minf;
        }
        T a_max = std::max_element(a.begin(), a.end());
        for(unsigned int i=0; i<a.size(); ++i) {
            a_copy[i] = b[i]*exp(a_copy[i] - a_max[i]);
        }
        return a_max + log(accumulate(a_copy.begin(), a_copy.end(), zero, std::plus<T>()));
    }*/

    template<class T, class InputIt>
    T LogSumExp(InputIt first, InputIt last) {
        T zero(0);
        T minf(-std::numeric_limits<double>::infinity());
        std::vector<T> a_copy(first, last);
        if(a_copy.size() == 1)
            return a_copy[0];
        T a_max = *std::max_element(a_copy.begin(), a_copy.end());
        for(auto it = a_copy.begin(); it !=a_copy.end(); ++it) {
            *it = exp(*it - a_max);
        }
        T val = accumulate(a_copy.begin(), a_copy.end());
        return a_max + log(val);
    }
    /*
    template<typename T>
    T LogSumExp(const & std::vector<T> a) {
        T zero(0);
        T minf(-std::numeric_limits<double>::infinity());
        std::vector<T> a_copy(a);
        T a_max = std::max_element(a.begin(), a.end());
        for(unsigned int i=0; i<a.size(); ++i) {
            a_copy[i] = exp(a_copy[i] - a_max[i]);
        }
        return a_max + log(accumulate(a_copy.begin(), a_copy.end(), zero, std::plus<T>()));
    }*/

    template<typename U>
    struct PhysToolsCast {
        template<typename T>
        struct temp {
            U operator()(T val) {
                return static_cast<U>(val);
            }
        };

        template<unsigned int Dim, typename T>
        struct temp<phys_tools::autodiff::FD<Dim, T>> {
            using result_type=phys_tools::autodiff::FD<Dim, T>;
            U operator()(result_type val) {
                return static_cast<U>(val.value());
            }
        };

        template<typename T>
        U operator()(T val) {
            return temp<T>()(val);
        }
    };
/*
    template<typename U, typename T>
    struct PhysToolsCast {
        U operator()(T val) {
            return static_cast<U>(val.value());
        }
    };

    template<typename U, unsigned int Dim, typename T>
    struct PhysToolsCast<U, phys_tools::autodiff::FD<Dim, T>> {
        using result_type=phys_tools::autodiff::FD<Dim, T>;
        U operator()(result_type val) {
            return static_cast<U>(val.value());
        }
    };
*/

    struct thorstenLikelihood {
        template<typename T>
        T operator()(double k, const std::vector<T>& raw_w, int n_events) const {
#ifdef __APPLE__
            //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW);
            _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW);
#elif __linux__
            //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
            feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

            double log_percentage = log(0.99);

            double prior_factor = 0;

            unsigned int count = 1;
            std::vector<unsigned int> kmc;
            std::vector<T> w;

            T zero(0);
            T wN(std::numeric_limits<double>::infinity());
            auto starting_it = raw_w.begin();
            for(; starting_it!=raw_w.end(); ++starting_it) {
                if(zero < *starting_it) {
                    break;
                }
            }

            // Get the unique weights and their counts
            if(starting_it != raw_w.end()) {
                T wi(*starting_it);
                wN = wi;
                for(auto it=starting_it+1; it!=raw_w.end(); ++it) {
                    if(*it != wi) {
                        w.push_back(wi);
                        wi = *it;
                        kmc.push_back(count);
                        count = 0;
                        if(wi < wN) {
                            wN = wi;
                        }
                    }
                    ++count;
                }
                if(count >= 1) {
                    w.push_back(wi);
                    kmc.push_back(count);
                    count = 0;
                }
            }

            double kmcs = std::accumulate(kmc.begin(), kmc.end(), (unsigned int)(0), std::plus<unsigned int>()) + prior_factor;

            std::vector<T> lgammak;
            std::vector<T> ldelta = {T(0)};

            auto eta_gen = boost::adaptors::transform(w, [&](T ww){return LogOneMinusX(wN / ww);});

            std::vector<T> eta(eta_gen.begin()+1, eta_gen.end());

            double alpha_factor = prior_factor / double(w.size());

            auto alpha_gen = boost::adaptors::transform(kmc, [&](unsigned int kmci)->double{return alpha_factor + kmci;});

            std::vector<double> alpha(alpha_gen.begin()+1, alpha_gen.end());

            T lwN = log(wN);
            auto lC_gen = boost::adaptors::transform(boost::irange((unsigned int)0, (unsigned int)w.size()), [&](unsigned int i)->T{return (lwN - log(w[i]))*alpha_gen[i];});
            //auto lC_gen = boost::adaptors::transform(w, alpha, [&](T wi, double alphai)->T{return (lwN - log(wi))*alphai;});

            T lC = accumulate(lC_gen.begin(), lC_gen.end());

            std::function<T(unsigned int)> get_log_gamma = [&] (unsigned int i)->T{
                std::function<double(double)> b_func([&](double kmci)->double{return kmci/double(i);});
                std::function<T(T)> a_func([&](T e)->T{return e*i;});
                auto b_gen = boost::adaptors::transform(alpha, b_func);
                auto a_gen = boost::adaptors::transform(eta, a_func);
                return LogSumExp<T>(a_gen.begin(), a_gen.end(), b_gen.begin(), b_gen.end());
            };

            std::function<void()> gen_next_log_delta = [&] ()->void{
                lgammak.push_back(get_log_gamma(lgammak.size()+1));
                unsigned int i = lgammak.size();
                std::function<T(unsigned int)> get_lg_plus_li = [&](unsigned int j)->T{return lgammak[j]+log(j+1);};
                auto lg_plus_li = boost::adaptors::transform(boost::irange((unsigned int)(0), (unsigned int)(i)), get_lg_plus_li);
                unsigned int ldelta_size = ldelta.size();
                std::function<T(unsigned int)> get_ld_plus_lg_plus_li = [&](unsigned int i)->T{return ldelta[ldelta_size-i-1]+lg_plus_li[i];};
                auto ld_plus_lg_plus_li = boost::adaptors::transform(boost::irange((unsigned int)0, ldelta_size), get_ld_plus_lg_plus_li);
                T next_ldelta = LogSumExp<T>(ld_plus_lg_plus_li.begin(), ld_plus_lg_plus_li.end()) - log(i);
                ldelta.push_back(next_ldelta);
            };

            std::function<double()> get_log_percentage = [&] ()->double{
                auto ldelta_gen = boost::adaptors::transform(ldelta, [&](T delta)->double{return PhysToolsCast<double>()(delta);});
                double ldelta_sum = LogSumExp<double>(ldelta_gen.begin(), ldelta_gen.end());
                double lC_double = PhysToolsCast<double>()(lC);
                double index = ldelta_sum + lC_double;
                std::cout << index << std::endl;
                return index;
            };

            while(get_log_percentage() < log_percentage) {
                gen_next_log_delta();
            }

            std::vector<T> terms(5);

            auto log_sum_terms_gen = boost::adaptors::transform(boost::irange((unsigned int)0, (unsigned int)ldelta.size()), [&](unsigned int i)->T{
                    terms[0] = ldelta[i];
                    terms[1] = lgamma(k+i+kmcs);
                    terms[2] = -i*lwN;
                    terms[3] = -lgamma(i + kmcs);
                    terms[4] = -i*LogOnePlusX(1.0/wN);
                    return accumulate(terms.begin(), terms.end());});

            T delta_sum = LogSumExp<T>(log_sum_terms_gen.begin(), log_sum_terms_gen.end());

            terms[0] = lC;
            terms[1] = -kmcs*lwN;
            terms[2] = -lgamma(1+k);
            terms[3] = -(k+kmcs)*LogOnePlusX(1.0/wN);
            terms[4] = delta_sum;

            return accumulate(terms.begin(), terms.end());

#ifdef __APPLE__
            _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~( _MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW));
#elif __linux__
            fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
        }
    };
/*
    struct thorstenLikelihood {
        template<typename T>
        T operator()(double k, const std::vector<T>& raw_w, int n_events) const {
#ifdef __APPLE__
            _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW);
#elif __linux__
            feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif

            unsigned int count = 1;
            std::vector<unsigned int> kmc;
            std::vector<T> w;

            std::cout << "RAW WEIGHTS BEGIN" << std::endl;

            for(unsigned int i=0; i<raw_w.size(); ++i) {
                std::cout << raw_w[i] << std::endl;

            }
            std::cout << "RAW WEIGHTS END" << std::endl;

            T zero(0);
            auto starting_it = raw_w.begin();
            for(; starting_it!=raw_w.end(); ++starting_it) {
                if(zero < *starting_it && SFromW(*starting_it) != zero && LogOnePlusX(*starting_it) != zero) {
                    break;
                }
            }
            if(starting_it != raw_w.end()) {
                //T wi = raw_w[0];
                T wi(*starting_it);
                for(auto it=starting_it+1; it!=raw_w.end(); ++it) {
                    if(*it != wi) {
                        w.push_back(wi);
                        wi = *it;
                        kmc.push_back(count);
                        if(count > 1) {
                            std::cout << "count: " << count << std::endl;
                        }
                        count = 0;
                    }
                    ++count;
                }
                if(count > 0) {
                    w.push_back(wi);
                    if(count > 1) {
                        std::cout << "count: " << count << std::endl;
                    }
                    kmc.push_back(count);
                }
                assert(std::distance(w.begin(), std::unique(w.begin(), w.end(), std::equal_to<T>())) == w.size());

                const unsigned int kmc_tot = std::accumulate(kmc.begin(), kmc.end(), (unsigned int)(0), std::plus<unsigned int>());
                const unsigned int M = w.size(); // Number of distinct weights


                //auto s_gen = boost::adaptors::transform(boost::irange((unsigned int)(0),(unsigned int)(M)),[&](unsigned int i){return SFromW(w[i]);});
                //std::vector<T> s(s_gen.begin(), s_gen.end());
                //T L(static_cast<T>(residuals::thorsten_fast<T>()(w,s,k)));
                //std::cout << "L = " << L << std::endl;
                T L(0);
                //std::cout << "L = " << L << std::endl;
                const std::vector<T> z(1, T(0));
                const std::vector<unsigned int> n(1, (unsigned int)(k+kmc_tot-1));
                std::vector<T> s(w.size());
                for(unsigned int i=0; i<M; ++i) {
                    s[i] = SFromW(w[i]);
                    T f = LogOnePlusX(w[i]);
                    L -= f;
                }
                L *= kmc_tot;
                //std::cout << "L *= " << kmc_tot << std::endl;
                //std::cout << "L = " << L << std::endl;
                const std::vector<unsigned int>& m = kmc;
                T f = static_cast<T>(residuals::contour_integral_bignum<10000, T>()(z, n, s, m));
                //T f = residuals::contour_integral<T>()(z, n, s, m);
                //T f = 1;
                //std::cout << "f = " << f << std::endl;
                assert(f > 0);
                T lf = log(f);
                L += lf;
                //std::cout << "L += " << lf << std::endl;
                //std::cout << "L = " << L << std::endl;

#ifdef __APPLE__
                _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~( _MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW));
#elif __linux__
                fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
                return L;
            }
            else {
#ifdef __APPLE__
                _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~( _MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW));
#elif __linux__
                fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
                return T(0);
            }
        }
    };
*/

    struct gammaPriorPoissonLikelihood {
        template<typename T>
        T operator()(double k, T const & alpha, T const & beta) {
            std::vector<T> items(5);
            items[0] = alpha*log(beta);
            items[1] = lgamma(k+alpha);
            items[2] = -lgamma(k+1);
            items[3] = -(k+alpha)*LogOnePlusX(beta);
            items[4] = -lgamma(alpha);
            return accumulate(items.begin(), items.end());
        }
    };

    struct thorstenEqualWeightsLikelihood {
        template<typename T>
        T operator()(double k, int n_events, T const & w_sum, T const & w2_sum) const {
            double kmc = n_events;
            if(kmc <= 0) {
                return T(0);
            }

            T one(1);
            T zero(0);
            if(w_sum == zero) {
                if(k == 0) {
                    return T(0);
                }
                else {
                    return T(-std::numeric_limits<double>::infinity());
                }
            }

            T s = sqrt(w2_sum);

            T alpha = kmc;
            T beta = kmc / s;
            T L = gammaPriorPoissonLikelihood()(k, alpha, beta);

            return L;
        }
    };

    struct SAYLikelihood {
        template<typename T>
        T operator()(double k, T const & w_sum, T const & w2_sum) const {
            if(w_sum <= 0 || w2_sum < 0) {
                return(k==0?0:-std::numeric_limits<T>::max());
            }

            if(w2_sum == 0) {
                return poissonLikelihood()(k, w_sum, w2_sum);
            }

            T one(1);
            T zero(0);
            if(w_sum == zero) {
                if(k == 0) {
                    return zero;
                }
                else {
                    return T(-std::numeric_limits<double>::infinity());
                }
            }

            T alpha = w_sum*w_sum/w2_sum + one;
            T beta = w_sum/w2_sum;
            T L = gammaPriorPoissonLikelihood()(k, alpha, beta);

            return L;
        }
    };

	//A bin type for keeping events in their histogram bins

	template<typename DataType>
	struct entryStoringBin{
	private:
		std::vector<DataType> data;
	public:
		entryStoringBin()=default;

		//allow constructing 'zero' values
		explicit entryStoringBin(double d){
			assert(d==0);
		}

		//allow assigning 'zero' values
		entryStoringBin& operator=(double d){
			assert(d==0);
            data.clear();
            return(*this);
		}

		entryStoringBin& operator+=(const DataType& e){
			data.push_back(e);
			return(*this);
		}
		entryStoringBin& operator+=(const entryStoringBin& other){
			data.insert(data.end(),other.data.begin(),other.data.end());
			return(*this);
		}
		//scaling does nothing
		entryStoringBin& operator*=(double scale){
			return(*this);
		}
		const entryStoringBin operator+(DataType a) const{
			return(entryStoringBin(*this)+=a);
		}
		const entryStoringBin operator+(const entryStoringBin& other) const{
			return(entryStoringBin(*this)+=other);
		}
		//scaling does nothing
		const entryStoringBin operator*(double scale) const{
			return(entryStoringBin(*this));
		}
		typename std::vector<DataType>::const_iterator begin() const{
			return(data.begin());
		}
		typename std::vector<DataType>::const_iterator end() const{
			return(data.end());
		}
		const std::vector<DataType>& entries() const{
			return(data);
		}
		size_t size() const{
			return(data.size());
		}
		bool empty() const{
			return(data.empty());
		}
        void clear() {
            data.clear();
        }
	};

} //namespace likelihood

namespace histograms{
	namespace detail{
		template<typename T>
		class histogramTraits<likelihood::entryStoringBin<T>>{
		public:
			struct amount{
				T value;
				amount(T v):value(v){}

				template <typename OtherType>
				amount(const OtherType& other):value(other.value){}
			};
			static void defaultData(likelihood::entryStoringBin<T>* data, unsigned int count){/* Nothing to do */}
		};
	} //namespace detail
} //namespace histograms

namespace likelihood{

	//for organizing simulation sets according to the same binning as the observed data

	template<typename Event, typename Container, int HistDim, typename HistStoreType>
	std::vector<phys_tools::histograms::histogram<HistDim,entryStoringBin<Event>>>
	binSimulation(phys_tools::histograms::histogram<HistDim,HistStoreType> observed, const std::vector<Container>& simulationSets){
		using phys_tools::histograms::histogram;
		using ResultHistType=histogram<HistDim,entryStoringBin<Event>>;
		using amount=typename ResultHistType::amount;
		using HistDimExt = phys_tools::histograms::detail::dimensionExtractor<histogram,HistDim,HistStoreType>;
		const unsigned int actualDim=HistDimExt::getDim(&observed);
		std::vector<ResultHistType> binnedSimulation;

		for(auto& simulation : simulationSets){
			binnedSimulation.emplace_back();
			auto& hist=binnedSimulation.back();
			phys_tools::histograms::detail::conditionalDimensionSetter<HistDim,HistStoreType>::setDimensions(hist,actualDim);
			for(unsigned int i=0; i<actualDim; i++)
				hist.setAxis(i, observed.getAxis(i)->copy());

			for(auto& e : simulation)
				hist.add(e,amount(std::cref(e)));
		}
		return(binnedSimulation);
	}

    template<unsigned int ...>
    struct parameters { };

    template<typename func, typename DataType, unsigned int ...S>
    DataType callFunc(func f, const std::vector<DataType>& params, parameters<S...>) {
       return f((params[S]) ...);
    }

    template<typename PriorPairs>
    class ArbitraryPriorSet {
    private:
        static constexpr unsigned int size=std::tuple_size<PriorPairs>::value;
        PriorPairs priorPairs;

        template<typename func, typename DataType>
        struct priorFunctionEvaluator {
            template<unsigned int S0, unsigned int... S1>
            DataType operator()(func f, const std::vector<DataType>& params, parameters<S0, S1...> indices) const {
                return callFunc(f, params, indices);
            }

            DataType operator()(func f, const std::vector<DataType>& params, parameters<> indices) const {
                return f(params);
            }
        };

        template<typename DataType, unsigned int index>
        struct priorEvaluator {
            DataType operator()(const PriorPairs& priors, const std::vector<DataType>& params) const {
                using priorPairType = typename std::tuple_element<size-index, PriorPairs>::type;
                const auto priorPair = std::get<size-index>(priors);
                using priorType = typename std::tuple_element<0, priorPairType>::type;
                const auto prior = std::get<0>(priorPair);
                using indicesType = typename std::tuple_element<1, priorPairType>::type;
                const auto indices = std::get<1>(priorPair);
                return priorFunctionEvaluator<priorType, DataType>()(prior, params, indices) + priorEvaluator<DataType, index-1>()(priors, params);
            }
        };

        template<typename DataType>
        struct priorEvaluator<DataType, 0> {
            DataType operator()(const PriorPairs& priors, const std::vector<DataType>& params) const {
                return(DataType(1));
            }
        };

    public:
        ArbitraryPriorSet(PriorPairs priorPairs):priorPairs(priorPairs){}

        template<typename DataType>
        DataType operator()(const std::vector<DataType>& params) const {
            return(priorEvaluator<DataType, size>()(priorPairs, params));
        }
    };

    template<typename T, typename U, unsigned int n, unsigned int size>
    struct zip_2_impl {
        auto operator()(T t, U u) ->
        decltype(std::tuple_cat(std::make_tuple(std::tuple<typename std::tuple_element<n,T>::type, typename std::tuple_element<n,U>::type>(std::get<n>(t), std::get<n>(u))), zip_2_impl<T,U,n+1,size>()(t, u))) {
            return std::tuple_cat(std::make_tuple(std::tuple<typename std::tuple_element<n,T>::type, typename std::tuple_element<n,U>::type>(std::get<n>(t), std::get<n>(u))), zip_2_impl<T,U,n+1,size>()(t, u));
        }
    };

    template<typename T, typename U, unsigned int size>
    struct zip_2_impl<T, U, size, size> {
        std::tuple<> operator()(T t, U u) {
            return std::make_tuple();
        }
    };

    template<typename T, typename U>
    auto zip_2(T t, U u) -> decltype(zip_2_impl<T,U,0,std::tuple_size<T>::value>()(t, u)) {
        return zip_2_impl<T,U,0,std::tuple_size<T>::value>()(t, u);
    }

    // Empty basecase
    template<typename...>
    std::tuple<> group_2(){return std::make_tuple();};

    // General case
    template<typename T, typename U, typename... Args>
    auto group_2(T t, U u, Args... args) ->
    decltype(std::tuple_cat(std::make_tuple(std::tuple<T,U>(t,u)), group_2(args...))) {
        return std::tuple_cat(std::make_tuple(std::tuple<T,U>(t,u)), group_2(args...));
    }

    template<typename Indices, typename... PriorTypes>
    auto makeArbitraryPriorSet(PriorTypes... priors) ->
    ArbitraryPriorSet<decltype(zip_2(std::tuple<PriorTypes...>(priors...), Indices()))> {
        return ArbitraryPriorSet<decltype(zip_2(std::tuple<PriorTypes...>(priors...), Indices()))>(
                zip_2(std::tuple<PriorTypes...>(priors...), Indices()));
    }

    template<typename... PriorTypes>
    auto makeArbitraryPriorSetFromArgs(PriorTypes... priors) ->
    ArbitraryPriorSet<decltype(group_2<PriorTypes...>(priors...))> {
        return ArbitraryPriorSet<decltype(group_2<PriorTypes...>(priors...))>(group_2(priors...));
    }

    template<typename Indices, typename... PriorTypes>
    struct ArbitraryPriorType {
        typedef decltype(makeArbitraryPriorSet<Indices>(std::declval<PriorTypes>()...)) type;
    };

	template<typename... PriorTypes>
	class FixedSizePriorSet{
	private:
		using Storage=std::tuple<PriorTypes...>;
		static constexpr unsigned int size=std::tuple_size<Storage>::value;
		Storage priors;

		template<typename DataType, int index>
		struct priorEvaluator{
			DataType operator()(const Storage& priors, const std::vector<DataType>& params) const{
				if(std::isnan(params[size-index]))
					return(std::numeric_limits<DataType>::quiet_NaN());
				return(std::get<size-index>(priors)(params[size-index])+priorEvaluator<DataType,index-1>()(priors,params));
			}
		};

		template<typename DataType>
		struct priorEvaluator<DataType,0>{
			DataType operator()(const Storage& priors, const std::vector<DataType>& params) const{
				return(DataType(1));
			}
		};

	public:
		FixedSizePriorSet(PriorTypes... priors):priors(priors...){}

		template<typename DataType>
		DataType operator()(const std::vector<DataType>& params) const{
			//if(params.size()!=size)
			//	std::cerr << "Got " << params.size() << " parameters but " << size << " were expected" << std::endl;
			assert(params.size()==size);
			return(priorEvaluator<DataType,size>()(priors,params));
		}
	};

	template<typename... PriorTypes>
	FixedSizePriorSet<PriorTypes...> makePriorSet(PriorTypes... priors){
		return(FixedSizePriorSet<PriorTypes...>(priors...));
	}

	//computes weights for observed data in the most obvious way: each event is given weight 1
	struct simpleDataWeighter{
		using result_type=double;
		template<typename Event>
		result_type operator()(const Event& e) const{ return(1.0); }
	};

    struct simpleDataWeightConstructor{
		using result_type=simpleDataWeighter;
		result_type operator()(const std::vector<double>&) const{ return(simpleDataWeighter()); }
	};

	// fundamental wrapper object for setting up likelihood fits for models by comparing observed data
	// events to simulated events
	template<typename Event, typename HistogramsType, typename DataWeighterConstructor, typename SimulationWeighterConstructor, typename CPrior, typename LFunc, int MaxDerivativeDimension=-1>
	struct LikelihoodProblem{
		using RawEvent = typename remove_reference_wrapper<Event>::type;
		static constexpr int DerivativeDimension=MaxDerivativeDimension;

		//static_assert(std::is_convertible<typename SimulationContainer::value_type,const RawEvent&>::value, "The simulation must be containers of objects convertible to Event&");

		//'observed' data
		HistogramsType observation;

		//expected data for every discrete nuisance parameter combination
		std::vector<HistogramsType> simulations;

		//prior distributions for continuous parameters
		CPrior continuousPrior;

		//prior distribution for discrete nuisance parameters
		std::vector<double> discretePrior;

		//object used to compute the weights of the observed data events
		//this may be a non-trivial calculation if the 'observed data' is a MC sample
		DataWeighterConstructor dataWeightC;

		//object used to construct the weighter for expected events from a set of model parameters
		//must be callable with an std::vector<double> containing the model parameters
		SimulationWeighterConstructor simWeightC;

		//the function used to compute the likelihood of observing a particular set of events in a bin
		//given a population of simulated, expected events in that bin, over some livetime
		LFunc likelihoodFunction;

		//the initial values for all of the model parameters,
		//including the index of the discrete nuisance parameter combination
		std::vector<double> parameterSeeds;

		//warning, evil
		mutable size_t lastBestDiscreteIndex;

		size_t evaluationThreadCount;

		LikelihoodProblem(HistogramsType observation, const std::vector<HistogramsType>& simulations,
						  CPrior continuousPrior, std::vector<double> discretePrior,
						  const DataWeighterConstructor& dataWeightC,
						  const SimulationWeighterConstructor& simWeightC,
						  const LFunc& likelihoodFunction,
						  std::vector<double> parameterSeeds,
						  unsigned int evaluationThreadCount=1):
		observation(observation),
		simulations(simulations),
		continuousPrior(continuousPrior),
		discretePrior(discretePrior),
		dataWeightC(dataWeightC),
		simWeightC(simWeightC),
		likelihoodFunction(likelihoodFunction),
		parameterSeeds(parameterSeeds),
		lastBestDiscreteIndex(0),
		evaluationThreadCount(evaluationThreadCount)
		{}

		template<typename AltLFunc>
		LikelihoodProblem<Event, HistogramsType, DataWeighterConstructor, SimulationWeighterConstructor, CPrior, AltLFunc, MaxDerivativeDimension>
		makeAlternateLikelihood(const AltLFunc& altlikelihoodFunction) const{
			using result_type=LikelihoodProblem<Event, HistogramsType, DataWeighterConstructor, SimulationWeighterConstructor, CPrior, AltLFunc, MaxDerivativeDimension>;
			return(result_type(observation,simulations,continuousPrior,discretePrior,dataWeightC,simWeightC,altlikelihoodFunction,parameterSeeds,evaluationThreadCount));
		}

		std::vector<double> getSeed() const{
			return(parameterSeeds);
		}

        void setSeed(const std::vector<double>& newSeed) {
			assert(newSeed.size() == parameterSeeds.size());
            parameterSeeds = newSeed;
		}

		void setObservation(const HistogramsType& newObs){
			observation=newObs;
		}

		void setEvaluationThreadCount(size_t count){
			if(count==0)
				throw std::runtime_error("At least one evaluation thread is required");
			evaluationThreadCount=count;
		}

		SimulationWeighterConstructor& getSimulationWeighterConstructor(){ return(simWeightC); }
		const SimulationWeighterConstructor& getSimulationWeighterConstructor() const{ return(simWeightC); }

		//evaluate the (non-prior) contribution to the likelihood from one observation,expectation histogram pair
		template<typename DataType, typename SimulationWeighter, typename TDataWeighter, typename HistogramType>
		void evaluateLikelihoodCore(const HistogramType& observation, SimulationWeighter& weighter, TDataWeighter& dweighter, const HistogramType& simulation,
								    std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions) const{
			const auto& likelihoodFunction=this->likelihoodFunction;
			auto dataWeightAccumulator=[&dweighter](double t, const Event& e){ return(t+dweighter(e)); };

			auto likelihoodContribution=[dataWeightAccumulator,&weighter,&simulation,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
				entryStoringBin<Event> obs=*it;

				double observationAmount=std::accumulate(obs.begin(),obs.end(),0.0,dataWeightAccumulator);

				std::vector<DataType> expectationWeights;
                std::vector<DataType> expectationSqWeights;
                DataType w;
                DataType w2;
                int n_events=0;
				typename HistogramType::const_iterator expIt=simulation.findBinIterator(it);
				if(expIt!=simulation.end()){
					/*{
						std::lock_guard<std::mutex> lck(printMtx);
						std::cout << "   exp coords:";
						for(unsigned int i=0; i<2; i++)
							std::cout << ' ' << expIt.getBinEdge(i);
						std::cout << std::endl;
					}*/
					const std::vector<Event>& exp=((entryStoringBin<Event>)*expIt).entries();
					expectationWeights.reserve(((entryStoringBin<Event>)*expIt).size());
					expectationSqWeights.reserve(((entryStoringBin<Event>)*expIt).size());
					for(const RawEvent& e : ((entryStoringBin<Event>)*expIt)){
                        w = weighter(e);
                        assert(w >= 0);
                        w2 = w * w;
                        assert(w2 >= 0);
                        assert(e.num_events > 0);
                        n_events += e.num_events;
						expectationWeights.push_back(w);
						expectationSqWeights.push_back(w2/e.num_events);
						if(std::isnan(expectationWeights.back()) || expectationWeights.back()<0.0){
							std::lock_guard<std::mutex> lck(printMtx);
							std::cout << "Bad weight: " << expectationWeights.back() << "\nEvent:\n" << e << std::endl;
							//std::cout << e.energy << ' ' << e.year << ' ' << expectationWeights.back() << std::endl;
						}
						if(std::isnan(expectationSqWeights.back()) || expectationSqWeights.back()<0.0){
							std::lock_guard<std::mutex> lck(printMtx);
							std::cout << "Bad weightSq: " << expectationSqWeights.back() << "\nEvent:\n" << e << std::endl;
							std::cout << "Bad weight: " << expectationWeights.back() << "\nEvent:\n" << e << std::endl;
							//std::cout << e.energy << ' ' << e.year << ' ' << expectationWeights.back() << std::endl;
						}
						if(std::isnan(e.num_events) || e.num_events<0.0){
							std::lock_guard<std::mutex> lck(printMtx);
							std::cout << "Bad num_events: " << e.num_events << "\nEvent:\n" << e << std::endl;
							//std::cout << e.energy << ' ' << e.year << ' ' << expectationWeights.back() << std::endl;
						}
						//std::cout << "    " << expectationWeights.back() << std::endl;
					}
				}

                //std::sort(expectationWeights.begin(), expectationWeights.end(), std::less<DataType>());

                auto w_sum = accumulate(expectationWeights.begin(), expectationWeights.end());
                auto w2_sum = accumulate(expectationSqWeights.begin(), expectationSqWeights.end());

                if(observationAmount > 0 && w_sum <= 0) {
                    std::cout << "BAD BIN" << std::endl;
                    std::cout << "Printing weights" << std::endl;
                    for(auto w : expectationWeights) {
                        std::cout << w << std::endl;
                    }
                    std::cout << "Printing events" << std::endl;
					for(const RawEvent& e : ((entryStoringBin<Event>)*expIt)){
                        std::cout << e << std::endl;
                    }
                    std::cout << "Printing data" << std::endl;
					for(auto e : obs) {
                        std::cout << e << std::endl;
                    }
                }

                auto contribution=likelihoodFunction(observationAmount,w_sum,w2_sum);

				/*{
					std::lock_guard<std::mutex> lck(printMtx);
					DataType expectationAmount=std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0));
					std::cout << "   obs coords:";
					for(unsigned int i=0; i<3; i++)
						std::cout << ' ' << it.getBinEdge(i);
					std::cout << ' ' << contribution
					<< ' ' << observationAmount << ' ' << expectationAmount;
					std::cout << " [";
					for(const auto& w : expectationWeights)
						std::cout << ' ' << w;
					std::cout << ']' << std::endl;
				}*/
                //std::cout << "contribution: " << contribution << std::endl;
				return(contribution);
			};
			auto likelihoodContributionNoObs=[&weighter,&observation,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
				auto obsIt=observation.findBinIterator(it);

				//only proceed if this bin does not exist in the observation
				if(obsIt==observation.end()){
					std::vector<DataType> expectationWeights;
					std::vector<DataType> expectationSqWeights;
                    DataType w;
                    DataType w2;
                    int n_events=0;
					const std::vector<Event>& exp=((entryStoringBin<Event>)*it).entries();
					expectationWeights.reserve(((entryStoringBin<Event>)*it).size());
					expectationSqWeights.reserve(((entryStoringBin<Event>)*it).size());
					for(const RawEvent& e : ((entryStoringBin<Event>)*it)) {
                        w = weighter(e);
                        assert(w >= 0);
                        w2 = w*w;
                        assert(w2 >= 0);
                        assert(e.num_events > 0);
                        n_events += e.num_events;
						expectationWeights.push_back(w);
					    expectationSqWeights.push_back(w2/e.num_events);
                    }

                    //std::sort(expectationWeights.begin(), expectationWeights.end(), std::less<DataType>());

					//auto contribution=likelihoodFunction(0,expectationWeights,n_events);
				    auto contribution=likelihoodFunction(0,accumulate(expectationWeights.begin(), expectationWeights.end()),accumulate(expectationSqWeights.begin(), expectationSqWeights.end()));
					/*{
						std::lock_guard<std::mutex> lck(printMtx);
						DataType expectationAmount=std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0));
						std::cout << "   exp coords:";
						for(unsigned int i=0; i<DataDimension; i++)
							std::cout << ' ' << it.getBinEdge(i);
						std::cout << ' ' << contribution
						<< ' ' << 0.0 << ' ' << expectationAmount;
						std::cout << " [";
						for(const auto& w : expectationWeights)
							std::cout << ' ' << w;
						std::cout << ']' << std::endl;
					}*/
                    //std::cout << "noobs contribution: " << contribution << std::endl;
					return(contribution);
				}
				/*{
					std::lock_guard<std::mutex> lck(printMtx);
					std::cout << "   exp coords:";
					for(unsigned int i=0; i<DataDimension; i++)
						std::cout << ' ' << it.getBinEdge(i);
					std::cout << ' ' << DataType(0)
					<< ' ' << 0.0 << ' ' << DataType(0);
					std::cout << " []" << std::endl;
				}*/
				return(0);
			};

			//handle all bins where there is observed data
			//std::cout << " bins with observations" << std::endl;
			for(auto it=observation.begin(), end=observation.end(); it!=end; it++)
				contributions.push_back(pool.enqueue(likelihoodContribution,it));

			//handle all bins where there is expected data, but which are not in the observation histogram
			//std::cout << " bins without observations" << std::endl;
			for(auto it=simulation.begin(), end=simulation.end(); it!=end; it++)
				contributions.push_back(pool.enqueue(likelihoodContributionNoObs,it));
		}

//        //evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
//        //[recursive/iterative case]
//		template<unsigned int Counter, typename DataType, typename SimulationWeighter>
//		void evaluateLikelihoodIterate(const HistogramsType& observation, SimulationWeighter weighter, const HistogramsType& simulation,
//								  std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
//			//evaluate for this pair
//			constexpr unsigned int idx=std::tuple_size<HistogramsType>::value - Counter;
//			evaluateLikelihoodCore(std::get<idx>(observation), weighter, std::get<idx>(simulation),
//								   printMtx, pool, contributions);
//			//evaluate for the next pair
//			evaluateLikelihoodIterate<Counter-1>(observation, weighter, simulation,
//												 printMtx, pool, contributions);
//		}
//        //evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
//        //[base case]
//		template<typename DataType, typename SimulationWeighter>
//		void evaluateLikelihoodIterate<0,DataType,SimulationWeighter>(const HistogramsType& observation, SimulationWeighter weighter, const HistogramsType& simulation,
//																 std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
//			//do nothing
//		}

		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
	    //[recursive/iterative case]
		template<unsigned int Counter, typename Likelihood, typename DataType, typename SimulationWeighter, typename TDataWeighter>
		struct evaluateLikelihoodIterator{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, TDataWeighter& dweighter, const HistogramsType& simulation,
							std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
				//evaluate for this pair
				constexpr unsigned int idx=std::tuple_size<HistogramsType>::value - Counter;
				like.evaluateLikelihoodCore(std::get<idx>(observation), weighter, dweighter, std::get<idx>(simulation),
									   printMtx, pool, contributions);
				//evaluate for the next pair
				evaluateLikelihoodIterator<Counter-1, Likelihood, DataType, SimulationWeighter, TDataWeighter>{}
					(like, observation, weighter, dweighter, simulation, printMtx, pool, contributions);
			}
		};
		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
	    //[base case]
		template<typename Likelihood, typename DataType, typename SimulationWeighter, typename TDataWeighter>
		struct evaluateLikelihoodIterator<0,Likelihood,DataType,SimulationWeighter,TDataWeighter>{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, TDataWeighter& dweighter, const HistogramsType& simulation,
							std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
				//do nothing
			}
		};

		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
		//[top level]
		template<typename DataType, typename SimulationWeighter, typename TDataWeighter>
		DataType evaluateLikelihood(SimulationWeighter weighter, TDataWeighter dweighter, const HistogramsType& simulation) const{
			std::mutex printMtx;
	        ThreadPool pool(evaluationThreadCount);
			std::vector<std::future<DataType>> contributions;
			const HistogramsType& observation=this->observation;

			//Iterate over the observeation and expectation tuples,
			//firing off tasks for computing the per-bin likelihoods.
			//Futures for these results will accumulate in contributions.
			evaluateLikelihoodIterator<std::tuple_size<HistogramsType>::value,decltype(*this),DataType,SimulationWeighter,TDataWeighter>{}
			(*this, observation, weighter, dweighter, simulation, printMtx, pool, contributions);

			//Iterate over the futures collecting the results and combining them when they are availiable.
			DataType llh(0.0);
			for(auto& future : contributions){
				auto contribution=future.get();
				llh+=contribution;
			}
            //std::cout << "LLH: " << llh << std::endl;
			return(llh);
		}

//		//evaluates the data (non-prior) portion of the likelihood
//		template<typename DataType, typename SimulationWeighter>
//		DataType evaluateLikelihood(SimulationWeighter weighter, const HistogramTypes& simulation) const{
//			DataType llh(0.0);
//
//			auto dataWeightAccumulator=[this](double t, const Event& e){ return(t+this->dataWeighter(e)); };
//
//			std::mutex printMtx;
//			const auto& likelihoodFunction=this->likelihoodFunction;
//			const auto& observation=this->observation;
//			auto likelihoodContribution=[&weighter,&simulation,&dataWeightAccumulator,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
//				entryStoringBin<Event> obs=*it;
//				double observationAmount=std::accumulate(obs.begin(),obs.end(),0.0,dataWeightAccumulator);
//
//				std::vector<DataType> expectationWeights;
//				typename HistogramType::const_iterator expIt=simulation.findBinIterator(it);
//				if(expIt!=simulation.end()){
//					/*std::cout << "   exp coords:";
//					for(unsigned int i=0; i<DataDimension; i++)
//						std::cout << ' ' << expIt.getBinEdge(i);
//					std::cout << std::endl;*/
//					const std::vector<Event>& exp=((entryStoringBin<Event>)*expIt).entries();
//					expectationWeights.reserve(((entryStoringBin<Event>)*expIt).size());
//					for(const RawEvent& e : ((entryStoringBin<Event>)*expIt)){
//						expectationWeights.push_back(weighter(e));
//						if(std::isnan(expectationWeights.back()) || expectationWeights.back()<0.0){
//							std::lock_guard<std::mutex> lck(printMtx);
//                            std::cout << "Bad weight: " << expectationWeights.back() << "\nEvent:\n" << e << std::endl;
//							//std::cout << e.energy << ' ' << e.year << ' ' << expectationWeights.back() << std::endl;
//						}
//						//std::cout << "    " << expectationWeights.back() << std::endl;
//					}
//				}
//
//				auto contribution=likelihoodFunction(observationAmount,expectationWeights);
//				/*{
//					std::lock_guard<std::mutex> lck(printMtx);
//					DataType expectationAmount=std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0));
//					std::cout << "   obs coords:";
//					for(unsigned int i=0; i<3; i++)
//						std::cout << ' ' << it.getBinEdge(i);
//					std::cout << ' ' << contribution
//					<< ' ' << observationAmount << ' ' << expectationAmount;
//					std::cout << " [";
//					for(const auto& w : expectationWeights)
//						std::cout << ' ' << w;
//					std::cout << ']' << std::endl;
//				}*/
//				return(contribution);
//			};
//			auto likelihoodContributionNoObs=[&weighter,&observation,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
//				auto obsIt=observation.findBinIterator(it);
//
//				//only proceed if this bin does not exist in the observation
//				if(obsIt==observation.end()){
//					std::vector<DataType> expectationWeights;
//					const std::vector<Event>& exp=((entryStoringBin<Event>)*it).entries();
//					expectationWeights.reserve(((entryStoringBin<Event>)*it).size());
//					for(const RawEvent& e : ((entryStoringBin<Event>)*it))
//						expectationWeights.push_back(weighter(e));
//
//					auto contribution=likelihoodFunction(0,expectationWeights);
//					/*{
//						std::lock_guard<std::mutex> lck(printMtx);
//						DataType expectationAmount=std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0));
//						std::cout << "   exp coords:";
//						for(unsigned int i=0; i<DataDimension; i++)
//							std::cout << ' ' << it.getBinEdge(i);
//						std::cout << ' ' << contribution
//						<< ' ' << 0.0 << ' ' << expectationAmount;
//						std::cout << " [";
//						for(const auto& w : expectationWeights)
//							std::cout << ' ' << w;
//						std::cout << ']' << std::endl;
//					}*/
//					return(contribution);
//				}
//				/*{
//					std::lock_guard<std::mutex> lck(printMtx);
//					std::cout << "   exp coords:";
//					for(unsigned int i=0; i<DataDimension; i++)
//						std::cout << ' ' << it.getBinEdge(i);
//					std::cout << ' ' << DataType(0)
//					<< ' ' << 0.0 << ' ' << DataType(0);
//					std::cout << " []" << std::endl;
//				}*/
//				return(0);
//			};
//
//			ThreadPool pool(evaluationThreadCount);
//			std::vector<std::future<DataType>> contributions;
//
//			std::vector<DataType> expectationWeights;
//
//			//handle all bins where there is observed data
//			//std::cout << " bins with observations" << std::endl;
//			for(auto it=observation.begin(), end=observation.end(); it!=end; it++)
//				contributions.push_back(pool.enqueue(likelihoodContribution,it));
//
//			//handle all bins where there is expected data, but which are not in the observation histogram
//			//std::cout << " bins without observations" << std::endl;
//			for(auto it=simulation.begin(), end=simulation.end(); it!=end; it++)
//				contributions.push_back(pool.enqueue(likelihoodContributionNoObs,it));
//
//			for(auto& future : contributions){
//				auto contribution=future.get();
//				llh+=contribution;
//			}
//
//			/*//finally, check the overflow bins
//			{
//				//std::cout << " overflow bins" << std::endl;
//				const auto& obs=observation.getOverflow().entries();
//				const auto& exp=simulation.getOverflow().entries();
//				if(!obs.empty() || !exp.empty()){
//					expectationWeights.clear();
//					double observationAmount=std::accumulate(obs.begin(),obs.end(),0.0,dataWeightAccumulator);
//					expectationWeights.reserve(exp.size());
//					for(const auto& e : exp)
//						expectationWeights.push_back(weighter(e));
//					//std::cout << "  " << observationAmount << " obs " << std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0),std::plus<DataType>()) << " exp";
//					llh+=likelihoodFunction(observationAmount,expectationWeights);
//					//std::cout << "  " << llh << std::endl;
//				}
//			}*/
//
//			return(llh);
//		}

		template<typename DataType, typename SimulationWeighter, typename TDataWeighter, typename DataHistogramType, typename ResultHistogramType>
		void evaluateLikelihoodContributionsCore(const DataHistogramType& observation, SimulationWeighter& weighter, TDataWeighter& dweighter, const DataHistogramType& simulation, ResultHistogramType& result) const{
			const auto& likelihoodFunction=this->likelihoodFunction;
			auto dataWeightAccumulator=[&dweighter](double t, const Event& e){ return(t+dweighter(e)); };

			using HistDimExt = phys_tools::histograms::detail::dimensionExtractor<phys_tools::histograms::histogram,DataHistogramType::dimensions,typename DataHistogramType::dataType>;
			const unsigned int histDimensions=HistDimExt::getDim(&observation);
			result.setUseContentScaling(false);

			auto likelihoodContribution=[dataWeightAccumulator,&weighter,&simulation,&likelihoodFunction,&result,histDimensions](typename DataHistogramType::const_iterator it)->void{
				entryStoringBin<Event> obs=*it;

				double observationAmount=std::accumulate(obs.begin(),obs.end(),0.0,dataWeightAccumulator);

				std::vector<DataType> expectationWeights;
				typename DataHistogramType::const_iterator expIt=simulation.findBinIterator(it);
				if(expIt!=simulation.end()){
					const std::vector<Event>& exp=((entryStoringBin<Event>)*expIt).entries();
					expectationWeights.reserve(((entryStoringBin<Event>)*expIt).size());
					for(const RawEvent& e : ((entryStoringBin<Event>)*expIt))
						expectationWeights.push_back(weighter(e));
				}
				DataType contribution=likelihoodFunction(observationAmount,expectationWeights);

				std::vector<double> coordinates(histDimensions,0);
				for(unsigned int i=0; i<histDimensions; i++)
					coordinates[i]=it.getBinCenter(i);
				result.add(coordinates.data(),contribution);
			};
			auto likelihoodContributionNoObs=[&weighter,&observation,&likelihoodFunction,&result,histDimensions](typename DataHistogramType::const_iterator it)->void{
				auto obsIt=observation.findBinIterator(it);

				DataType contribution=0;
				//only proceed if this bin does not exist in the observation
				if(obsIt==observation.end()){
					std::vector<DataType> expectationWeights;
					const std::vector<Event>& exp=((entryStoringBin<Event>)*it).entries();
					expectationWeights.reserve(((entryStoringBin<Event>)*it).size());
					for(const RawEvent& e : ((entryStoringBin<Event>)*it))
						expectationWeights.push_back(weighter(e));
					contribution=likelihoodFunction(0,expectationWeights);
				}

				std::vector<double> coordinates(histDimensions,0);
				for(unsigned int i=0; i<histDimensions; i++)
					coordinates[i]=it.getBinCenter(i);
				result.add(coordinates.data(),contribution);
			};

			//handle all bins where there is observed data
			//std::cout << " bins with observations" << std::endl;
			for(auto it=observation.begin(), end=observation.end(); it!=end; it++)
				likelihoodContribution(it);

			//handle all bins where there is expected data, but which are not in the observation histogram
			//std::cout << " bins without observations" << std::endl;
			for(auto it=simulation.begin(), end=simulation.end(); it!=end; it++)
				likelihoodContributionNoObs(it);
		}

		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
		//[recursive/iterative case]
		template<unsigned int Counter, typename Likelihood, typename DataType, typename SimulationWeighter, typename ResultHistogramType>
		struct evaluateLikelihoodContributionsIterator{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, const HistogramsType& simulation,
							ResultHistogramType& result){
				//evaluate for this pair
				constexpr unsigned int idx=std::tuple_size<HistogramsType>::value - Counter;
				like.template evaluateLikelihoodContributionsCore<DataType>(std::get<idx>(observation), weighter, std::get<idx>(simulation),std::get<idx>(result));
				//evaluate for the next pair
				evaluateLikelihoodContributionsIterator<Counter-1, Likelihood, DataType, SimulationWeighter, ResultHistogramType>{}
				(like, observation, weighter, simulation, result);
			}
		};
		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
		//[base case]
		template<typename Likelihood, typename DataType, typename SimulationWeighter, typename ResultHistogramType>
		struct evaluateLikelihoodContributionsIterator<0,Likelihood,DataType,SimulationWeighter,ResultHistogramType>{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, const HistogramsType& simulation,
							ResultHistogramType& result){
				//do nothing
			}
		};

		///Fill a collection of histograms with the per-bin likelihood contributions
		///\param weighter the weighter which transforms simulation into the expectation
		///\param simulation the tuple of histograms of simulated data which forms the basis of the expectation
		///\param results the tuple of histograms into which the likelihood contributions will be stored
		template<typename DataType, typename SimulationWeighter, typename ResultsType>
		void evaluateLikelihoodContributions(SimulationWeighter weighter, const HistogramsType& simulation, ResultsType& results) const{
			evaluateLikelihoodContributionsIterator<std::tuple_size<HistogramsType>::value,decltype(*this),DataType,SimulationWeighter,ResultsType>{}
			(*this, observation, weighter, simulation, results);
		}
		//return a pair consting of:
		//	-a histogram whose bins contain the likelihood contributions from the corresponding bins in the data
		//	-a value for the likelihood contribution of the under- and overflow bins

	private:

		template<typename DataType>
		std::pair<std::vector<DataType>,bool> fixUpParams(const std::vector<DataType>& rawParams) const{
			std::vector<DataType> params;
			bool scanDNuisance=false;
			//for a fully specified set of paramters just copy wholesale
			if(rawParams.size()==parameterSeeds.size())
				params.insert(params.end(),rawParams.begin(),rawParams.end());
			//a minimizer can't know about the discrete nuisance parameter space
			//so it will call with one fewer parameter
			else if(rawParams.size()==parameterSeeds.size()-1){
				scanDNuisance=true;
				params.insert(params.end(),rawParams.begin(),rawParams.end());
				params.push_back(0);
			}
			//any other number of parameters is an error
			else
				throw std::domain_error("Incorrect number of model parameters for likelihood evaluation");
			return(std::make_pair(params,scanDNuisance));
		}

	public:

//		///Computes the number of terms contributing to the likelihood
//		unsigned int likelihoodTerms() const{
//			unsigned int nTerms=0;
//			bool addedZero=false;
//
//			const auto& simulation=simulations.front();
//
//			for(auto it=observation.begin(), end=observation.end(); it!=end; it++){
//				const entryStoringBin<Event>& obs=*it;
//				if(obs.size()!=0)
//					nTerms++;
//				else if(!addedZero){
//					nTerms++;
//					addedZero=true;
//				}
//				/*auto expIt=simulation.findBinIterator(it);
//				if(expIt!=simulation.end())
//					nTerms++;*/
//			}
//
//			for(auto it=simulation.begin(), end=simulation.end(); it!=end; it++){
//				auto obsIt=observation.findBinIterator(it);
//				//only proceed if this bin does not exist in the observation
//				if(obsIt==observation.end()){
//					nTerms++;
//					addedZero=true;
//					break;
//				}
//			}
//
//			return(nTerms);
//		}

		template<typename DataType>
		DataType evaluateLikelihood(const std::vector<DataType>& rawParams, bool includePriors=true) const{
            //for(auto i: rawParams) {
            //    std::cout << i << ", ";
            //}
            //std::cout << std::endl;
			std::vector<DataType> params;
			bool scanDNuisance=false;
			//TODO: fix this once there are discrete variations
			//std::tie(params,scanDNuisance)=fixUpParams(rawParams);
			params=rawParams;
			params.push_back(DataType(0));
			std::vector<DataType> continuousParams(params.begin(),params.begin()+(params.size()-1));

			/*std::cout << " Evaluating at (";
			for(size_t i=0; i<continuousParams.size(); i++){
				if(i)
					std::cout << ',';
				std::cout << continuousParams[i];
			}
			std::cout << "): ";
			std::cout.flush();*/


			//size_t minDN=(scanDNuisance?0:(size_t)params.back());
			//size_t maxDN=(scanDNuisance?simulations.size():(size_t)params.back()+1);
			size_t minDN=0, maxDN=1; //TODO: fix this once there are discrete variations

			DataType cPrior=continuousPrior.template operator()<DataType>(continuousParams);

			//bail out early if we're in a disallowed part of the parameter space
			//if(cPrior<=0 || std::isnan(cPrior)){
			if(std::isnan(cPrior)){
				//std::cout << "Out-of-bounds prior prob" << std::endl;
				return(-std::numeric_limits<DataType>::max());
			}
			//std::cout << "(prior=" << cPrior << ") ";
			//cPrior=log(cPrior); //convert to log(P) for rest of calculation

			//compute the actual llh for every necessary discrete nuisance index
			//(of which there will only be one unless we've been called by a minimizer)
			//and keep the best value
			DataType bestLLH=-std::numeric_limits<DataType>::max();
			for(size_t dn=minDN; dn<maxDN; dn++){
				DataType llh=(includePriors?cPrior+discretePrior[dn]:0);
				params.back()=dn;
				llh+=evaluateLikelihood<DataType>(simWeightC(continuousParams),dataWeightC(continuousParams),simulations[dn]);

				if(llh>bestLLH){
					bestLLH=llh;
					lastBestDiscreteIndex=dn;
				}
			}
			//std::cout << " llh: " << bestLLH << std::endl;
			return(bestLLH);
		}

//        ///\return a pair consitsting of the LLH per bin in the form of a histogram (with the same dimensions and axes as the data) and the log probability due to the priors
//		template<typename DataType>
//		std::pair<phys_tools::histograms::histogram<DataDimension,DataType>,DataType> evaluateLikelihoodHistogram(const std::vector<DataType>& rawParams, bool includePrior=false) const{
//			std::vector<DataType> params=rawParams;
//			bool scanDNuisance=false;
//			/*std::tie(params,scanDNuisance)=fixUpParams(rawParams);
//
//			if(scanDNuisance){
//				if(simulations.size()==1)
//					params.back()=0;
//				else{
//					//dumb, but this shouldn't be called frequently
//					//find out which discrete index to use
//					evaluateLikelihood(rawParams);
//					params.back()=lastBestDiscreteIndex;
//				}
//			}*/
//
//			std::vector<DataType> continuousParams(params.begin(),params.begin()+(params.size()-1));
//			DataType cPrior=continuousPrior.template operator()<DataType>(continuousParams);
//
//			//bail out early if we're in a disallowed part of the parameter space
//			if(cPrior<=0 || std::isnan(cPrior)){
//				//std::cout << "Out-of-bounds prior prob" << std::endl;
//				return(std::make_pair(phys_tools::histograms::histogram<DataDimension,DataType>(),-std::numeric_limits<DataType>::max()));
//			}
//			cPrior=log(cPrior); //convert to log(P) for rest of calculation
//			std::pair<phys_tools::histograms::histogram<DataDimension,DataType>,DataType> results=
//				evaluateLikelihoodContributions<DataType>(simWeightC(continuousParams),simulations[size_t(params.back())]);
//
//			if(includePrior){
//				for(auto bin : results.first)
//					bin+=cPrior;
//				results.second+=cPrior;
//			}
//			return(results);
//		}

		double operator()(const std::vector<double>& params) const {
			return(-evaluateLikelihood<double>(params));
		}

		//computes the gradient of the negative log likelihood with respect to the model parameters
		std::vector<double> gradient(const std::vector<double>& p) const {
			const size_t size=p.size();
            using GradType=autodiff::FD<MaxDerivativeDimension>;
			std::vector<GradType> params(size);
			for(size_t i=0; i<size; i++)
				params[i]=GradType(p[i],i);
			GradType result=evaluateLikelihood<GradType>(params);
			std::vector<double> grad(size);
			for(unsigned int i=0; i<size; i++)
				grad[i]=-result.derivative(i); //note negation!
			return(grad);
		}

		std::vector<std::vector<double>> hessian(const std::vector<double>& p) const{
			using dtype=autodiff::FD<MaxDerivativeDimension>;
			using htype=autodiff::FD<MaxDerivativeDimension,dtype>;
			const size_t size=p.size();
			std::vector<htype> params(size);
			for(size_t i=0; i<size; i++)
				params[i]=htype(dtype(p[i],i),i);
			htype result=evaluateLikelihood<htype>(params);
			std::vector<std::vector<double>> h(size);
			for(unsigned int i=0; i<size; i++){
				h[i].resize(size);
				for(unsigned int j=0; j<size; j++)
					h[i][j]=-result.derivative(i).derivative(j); //note negation!
			}
			return(h);
		}
	};

	//helper function for constructing a LikelihoodProblem without having to write out all of the template parameters
	//Note that the order of the template parameters has been shuffled, to put those which cannot be deduced from the
	//arguments (Event, DataDimension, and MaxDerivativeDimension) first so that they can be specified while the rest
	//are left to be deduced, while the order of the function arguments is the same as for the LikelihoodProblem constructor
	template<typename Event, int MaxDerivativeDimension=-1, typename HistogramsType, typename DataWeightConstructor, typename SimulationWeighterConstructor, typename CPrior, typename LFunc,
	         typename LikelihoodType=LikelihoodProblem<Event,HistogramsType,DataWeightConstructor,SimulationWeighterConstructor,CPrior,LFunc,MaxDerivativeDimension>>
	LikelihoodType makeLikelihoodProblem(HistogramsType observation, const std::vector<HistogramsType>& simulations,
                                         CPrior continuousPrior, std::vector<double> discretePrior,
                                         const DataWeightConstructor& dataWeightC,
                                         const SimulationWeighterConstructor& simWeightC,
                                         const LFunc& likelihoodFunction,
                                         std::vector<double> parameterSeeds,
                                         unsigned int evaluationThreadCount=1){
		return(LikelihoodType(observation,simulations,
							  continuousPrior,discretePrior,
							  dataWeightC,simWeightC,likelihoodFunction,
							  parameterSeeds,evaluationThreadCount));
	}

	template<typename FuncType>
	class BFGS_Function : public lbfgsb::BFGS_FunctionBase{
	private:
		FuncType func;
	public:
		BFGS_Function(FuncType f):func(f){}
		virtual double evalF(std::vector<double> x) const{
            //std::cout << "EvaluateF: (";
            //for(auto xi : x) std::cout << xi << ", ";
            //std::cout << ")\n";
			return(-func.template evaluateLikelihood<double>(x));
		}
		virtual std::pair<double,std::vector<double>> evalFG(std::vector<double> x) const{
            //std::cout << "EvaluateFG: (";
            //for(auto xi : x) std::cout << xi << ", ";
            //std::cout << ")\n";
			const size_t size=x.size();
            using GradType=autodiff::FD<FuncType::DerivativeDimension>;
			std::vector<GradType> params(size);
			for(size_t i=0; i<size; i++)
				params[i]=GradType(x[i],i);
			GradType result=func.template evaluateLikelihood<GradType>(params);
			std::vector<double> grad(size);
			for(unsigned int i=0; i<size; i++)
				grad[i]=-result.derivative(i); //note negation!
			return(std::make_pair(-result.value(),grad));
		}
	};

	template<typename Container, typename RNG>
	std::vector<typename Container::value_type>
	generateSample(const std::vector<double>& weights, const Container& data, double quantity, RNG& eng){
		assert(!data.empty());
		assert(weights.size()==data.size());
		std::discrete_distribution<size_t> dist(weights.begin(),weights.end());
		//decide how many events to sample
		size_t s_quantity = std::poisson_distribution<size_t>(quantity)(eng);
		//generate the sample
		std::vector<typename Container::value_type> sample;
		sample.reserve(s_quantity);
		for(unsigned int i=0; i<s_quantity; i++){
			size_t idx=dist(eng);
			//std::cout << " got index " << idx << std::endl;
			sample.push_back(std::cref(data[idx]));
		}
		return(sample);
	}

	//
	template<typename Container, typename Weighter, typename RNG>
	std::vector<typename Container::value_type>
	generateSample(Weighter w, const Container& data, double quantity, RNG& eng){
		assert(!data.empty());
		std::vector<double> weights;
		weights.reserve(data.size());
		double wsum=0;
		for(const auto& e : data){
			double weight=w(e);
			if(std::isnan(weight) || std::isinf(weight) || weight<0)
				weight=0;
			wsum+=weight;
			weights.push_back(weight);
		}
		for(auto& weight : weights)
			weight/=wsum;
		return(generateSample(weights,data,quantity,eng));
	}

    ///A prior which gives uniform probability for the parameter to be anywhere
    ///in a fixed range, but forbids it being outside that range.
	struct UniformPrior{
	private:
		double min, max;
	public:
		UniformPrior(double min=-std::numeric_limits<double>::infinity(),
					 double max=std::numeric_limits<double>::infinity()):
		min(min),max(max){}

		template<typename DataType>
		DataType operator()(DataType x) const{
			if(x<min || x>max)
				return(DataType(-std::numeric_limits<double>::infinity()));
			return(DataType(0.0));
		}
	};

    struct PowerPrior{
	private:
		double power, min, max, lnorm;
	public:
		PowerPrior(double power=0.0, double min=-std::numeric_limits<double>::infinity(),
				double max=std::numeric_limits<double>::infinity()):
		        power(power), min(min), max(max) {
            assert(min >= 0.0);
            if(min == 0.0 and power <= 0.0) {
                lnorm = 0.0;
            }
            else if(std::isfinite(max) && std::isfinite(min)) {
                lnorm = log((std::pow(min, power+1.0) - std::pow(max, power+1.0))/(power+1.0));
            }
            else if(std::isfinite(min) && power < 1.0) {
                lnorm = log(std::pow(min, power+1.0)/(power+1.0));
            }
            else {
                lnorm = 0.0;
            }
        }

		template<typename DataType>
		DataType operator()(DataType x) const{
			if(x<min || x>max)
				return(DataType(-std::numeric_limits<double>::infinity()));
            else if(power == 0.0)
                return 0.0;
			return power*log(x) - lnorm;
		}
	};

	struct GaussianPrior{
	private:
		double mean;
		double stddev;
		double norm;
	public:
		GaussianPrior(double mean, double stddev):
		mean(mean),stddev(stddev),
		norm(boost::math::constants::one_div_root_two_pi<double>()/stddev) {
            if(std::isinf(stddev) || std::isnan(stddev)) {
                norm = 0.0;
            }
        }

		template<typename DataType>
		DataType operator()(DataType x) const{
            if(norm == 0.0)
                return 0.0;
			DataType z=(x-mean)/stddev;
			return(log(norm)-z*z/2);
		}
	};

    struct LimitedGaussianPrior{
	private:
		UniformPrior limits;
		GaussianPrior prior;
	public:
		LimitedGaussianPrior(double mean, double stddev, double min, double max):
		limits(min,max),prior(mean,stddev){}

		template<typename DataType>
		DataType operator()(DataType x) const{
			return(limits(x)+prior(x));
		}
	};

    struct Gaussian2DPrior{
	private:
		double mean0;
		double mean1;
		double stddev0;
		double stddev1;
		double correlation;
		double lnorm;
		double prefactor;
	public:
		Gaussian2DPrior(double mean0, double mean1, double stddev0, double stddev1, double correlation):
		mean0(mean0),mean1(mean1),stddev0(stddev0),stddev1(stddev1),correlation(correlation),
		lnorm(log(boost::math::constants::one_div_two_pi<double>()/(stddev0*stddev1*sqrt(1.0-correlation*correlation)))),
        prefactor(-1.0/(2.0*(1.0-correlation*correlation))){
            if(std::isinf(stddev0) || std::isinf(stddev0) || std::isnan(stddev0) || std::isnan(stddev1) || std::isnan(correlation)) {
                lnorm = 0.0;
                prefactor = 0.0;
            }
        }

		template<typename DataType>
		DataType operator()(DataType x0, DataType x1) const{
            if(prefactor == 0.0)
                return lnorm;
			DataType z0=(x0-mean0)/stddev0;
			DataType z1=(x1-mean1)/stddev1;
			return lnorm + prefactor*(z0*z0 + z1*z1 - 2.0*correlation*z0*z1);
		}
	};

	struct LimitedGaussian2DPrior{
	private:
		UniformPrior limits0;
		UniformPrior limits1;
		Gaussian2DPrior prior;
	public:
		LimitedGaussian2DPrior(double mean0, double mean1, double stddev0, double stddev1, double correlation, double min0, double max0, double min1, double max1):
		limits0(min0,max0),limits1(min1,max1),prior(mean0,mean1,stddev0,stddev1,correlation){}

		template<typename DataType>
		DataType operator()(DataType x) const{
			return(limits(x)+prior(x));
		}
	};

    struct ZigZagPrior{
	private:
		double point;
		double scale;
        double m;
	public:
		ZigZagPrior(double point, double scale, bool small):
        point(point), scale(scale), m(int(small)*2.0 - 1.0) {}

		template<typename DataType>
		DataType operator()(DataType x) const{
            return log((tanh(scale*(point-x)*m)+1.)/2.+1e-18)-exp(-scale*(point-x)*m);
		}
	};

} //namespace likelihood
} //namespace phys_tools

#endif
