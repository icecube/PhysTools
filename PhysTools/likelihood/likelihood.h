#ifndef LF_LIKELIHOOD_H
#define LF_LIKELIHOOD_H

#include <cmath>
#include <deque>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include <boost/mpl/if.hpp>

#include "../optimization/lbfgsb/lbfgsb.h"

#include "../histogram.h"

#include "../autodiff.h"
#include "../root_finding.h"
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

// Implementation details for likelihoods
namespace detail {

	// Check for a bad weight
	// Returns true if the weight is bad
	template <typename T>
	struct CheckBadWeight {
		bool operator()(const T& w) const{
			return std::isnan(w) || w<0.0;
		}
	};

	// Special case for autodiff
	// Returns true if the weight or its derivative is bad
	template<unsigned int Dim, typename T>
	struct CheckBadWeight<phys_tools::autodiff::FD<Dim, T>> {
		using result_type=phys_tools::autodiff::FD<Dim, T>;
		bool operator()(const result_type& val) const{
			bool res = CheckBadWeight<T>()(val.value());
			for(unsigned int i=0; i<Dim && !res; i++)
				res |= std::isnan(val.derivative(i));
			return res;
		}
	};

	// Zero weights (please don't use this)
	struct NullWeighter {
		template <typename Event, typename T>
		T operator()(const Event& e) const{
			return T(0);
		}
	};

	// Estimate the weight uncertainty as the square of the weight
	struct WeightsSqUncertaintyWeighter {
		template <typename Event, typename T>
		T operator()(const Event& e, const T& w) const{
			return w*w;
		}
	};

	// Properties of the weighted expectation
	// Weights and their uncertainties
	template<typename DataType>
	struct ExpectationProperties {
		std::vector<DataType> weights;
		std::vector<DataType> weightUncertainties;
		ExpectationProperties():weights(),weightUncertainties(){};
		ExpectationProperties(std::vector<DataType> w):weights(w),weightUncertainties(){};
		ExpectationProperties(std::vector<DataType> w, std::vector<DataType> w2):weights(w),weightUncertainties(w2){};
	};

	// Computes weights and weight uncertainties for a set of events
	template<typename Weighter, typename UncertaintyWeighter>
	struct WeighterUncertaintyWeighterPair {
		const Weighter weighter;
		const UncertaintyWeighter uWeighter;

		WeighterUncertaintyWeighterPair(Weighter weighter, UncertaintyWeighter uWeighter):weighter(weighter),uWeighter(uWeighter){};

		template<typename DataType, typename Event>
		ExpectationProperties<DataType> operator()(std::vector<Event> const & events) const{
			using RawEvent = typename remove_reference_wrapper<Event>::type;
			std::vector<DataType> expectationWeights;
			std::vector<DataType> expectationWeightUncertainties;
			DataType w;
			DataType uW;
			expectationWeights.reserve(events.size());
			expectationWeightUncertainties.reserve(events.size());
			for(const RawEvent& e : events){
				w = weighter(e);
				uW = uWeighter(e,w);
				expectationWeights.push_back(w);
				expectationWeightUncertainties.push_back(uW);
			}
			return ExpectationProperties<DataType>(expectationWeights, expectationWeightUncertainties);
		}
	};

	// Group param-pack of into pairs of weighterMakers and uncertaintyWeighters
	// Default uncertaintyWeighter is WeightsSqUncertaintyWeighter

	// Empty basecase
	template<typename...>
	std::tuple<> group_weighters(){return std::make_tuple();};

	// Single element basecase
	template<typename T>
	std::tuple<std::tuple<T,WeightsSqUncertaintyWeighter>> group_weighters(T t) {
		return std::make_tuple(std::tuple<T,WeightsSqUncertaintyWeighter>(t, WeightsSqUncertaintyWeighter()));
	};

	// General case
	template<typename T, typename U, typename... Args>
	decltype(std::tuple_cat(std::make_tuple(std::tuple<T,U>(std::declval<T>(),std::declval<U>())), group_weighters(std::declval<Args>()...))) group_weighters(T t, U u, Args... args) {
		return std::tuple_cat(std::make_tuple(std::tuple<T,U>(t,u)), group_weighters(args...));
	};

	// Make a useable weighter from  a set of continuousParams, a weighterMaker, and an uncertaintyWeighter.
	// Returned weighter acts on a vector of events so the cost of using std::function is small
	template<typename Event, typename DataType, typename WeighterConstructor, typename UncertaintyWeighter>
	std::function<ExpectationProperties<DataType>(const std::vector<Event>&)> make_property_weighter(std::tuple<WeighterConstructor, UncertaintyWeighter> uwc, std::vector<DataType> continuousParams) {
		WeighterConstructor wc = std::get<0>(uwc);
		UncertaintyWeighter uw = std::get<1>(uwc);
        auto weighter = wc(continuousParams);
		WeighterUncertaintyWeighterPair<decltype(wc(std::declval<std::vector<DataType>>())),UncertaintyWeighter> wp(std::move(weighter), uw);
		return [=](const std::vector<Event>& events)->ExpectationProperties<DataType>{
			return wp.template operator()<DataType>(events);
		};
	}

	template<typename... Weighters>
	struct SwitchableWeighter{
		static_assert(sizeof...(Weighters)>0,"At least one weighter implementation is required");

		std::tuple<Weighters...> implementations;
		unsigned int selection;

		SwitchableWeighter(Weighters&&... ws):implementations(ws...),selection(0){}

		template<unsigned int index, typename Event>
		struct evaluator{
			double operator()(const SwitchableWeighter& sw, const Event& e) const{
				if(index==sw.selection)
					return(std::get<index>(sw.implementations)(e));
				return(evaluator<index+1,Event>()(sw,e));
			}
		};

		template<typename Event> //invalid index case to stop template recursion
		struct evaluator<sizeof...(Weighters),Event>{
			double operator()(const SwitchableWeighter& sw, const Event& e) const{
				__builtin_unreachable();
			}
		};

        template<typename Event>
		double operator()(const Event& e) const{
			return(evaluator<0,Event>()(*this,e));
		}

		void setWeighter(unsigned int index){
			if(index>=sizeof...(Weighters))
				throw std::logic_error("Out of range weighter index");
			selection=index;
		}
	};

	// A Collection of weighterMakers and uncertaintyWeighters stored as pairs(using std::tuple)
	// operator() produces a weighter usable on a set of events
	// setWeighter switches between pairs used to make the weighter
	template<typename... Weighters>
	struct WeighterCollection {
		static_assert(sizeof...(Weighters)>0,"At least one weighter implementation is required");

		typedef decltype(group_weighters(std::declval<Weighters>()...)) tuple_type;
		tuple_type implementations;

		unsigned int selection;

		WeighterCollection(Weighters... ws):implementations(group_weighters(ws...)),selection(0){}

		template<unsigned int index, typename Event, typename DataType>
		struct evaluator{
			std::function<ExpectationProperties<DataType>(const std::vector<Event>&)> operator()(const WeighterCollection& sw, const std::vector<DataType>& continuousParams) const{
				if(index==sw.selection) {
					return make_property_weighter<Event>(std::get<index>(sw.implementations), continuousParams);
				}
				return(evaluator<index+1,Event,DataType>()(sw,continuousParams));
			}
		};

		template<typename Event, typename DataType>
		struct evaluator<std::tuple_size<tuple_type>::value, Event, DataType>{
			std::function<ExpectationProperties<DataType>(const std::vector<Event>&)> operator()(const WeighterCollection& sw, const std::vector<DataType> & continuousParams) const{
				throw std::logic_error("Out of range weighter index");
			 }
		};

		template<typename Event, typename DataType>
		std::function<ExpectationProperties<DataType>(const std::vector<Event>&)> operator()(const std::vector<DataType> continuousParams) const{
			return(evaluator<0,Event,DataType>()(*this,continuousParams));
		}

		void setWeighter(unsigned int index){
			if(index>=sizeof...(Weighters))
				throw std::logic_error("Out of range weighter index");
			selection=index;
		}
	};

	// Compute the sum using the Kahan Summation Algorithm
	// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
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

		using std::abs;
		using std::log;
		using std::pow;
		if (abs(x) > 1e-4)
		{
			// x is large enough that the obvious evaluation is OK
			return log(1.0 + x);
		}

		// Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
		// Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
		return x-pow(x,2.0)/2.0+pow(x,3.0)/3.0-pow(x,4.0)/4.0;
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

		using std::abs;
		using std::log;
		using std::pow;
		if (abs(x) > 1e-4)
		{
			// x is large enough that the obvious evaluation is OK
			return log(1.0 - x);
		}

		// Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
		// Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
		return -x-pow(x,2.0)/2.0-pow(x,3.0)/3.0-pow(x,4.0)/4.0;
	};

	// A poisson likelihood with gamma distribution prior
	struct gammaPriorPoissonLikelihood {
		template<typename T>
		T operator()(double k, T const & alpha, T const & beta) {
			std::array<T,5> items;
			items[0] = alpha*log(beta);
			items[1] = lgamma(k+alpha);
			items[2] = -lgamma(k+1);
			items[3] = -(k+alpha)*LogOnePlusX(beta);
			items[4] = -lgamma(alpha);
			return accumulate(items.begin(), items.end());
		}
	};

} // namespace detail

	// Convenience function for creating a WeighterCollection
	template<typename... Weighters>
	phys_tools::likelihood::detail::WeighterCollection<Weighters...> make_weighter_collection(Weighters... ws) {
		return phys_tools::likelihood::detail::WeighterCollection<Weighters...>(ws...);
	}

	//Frequently used likelihood functions

	struct poissonLikelihood{
		template <typename T>
		T operator()(double dataCount, likelihood::detail::ExpectationProperties<T> const & expect) const {
			std::vector<T> const & simulationWeights = expect.weights;
			T lambda=std::accumulate(simulationWeights.begin(),simulationWeights.end(),T(0),std::plus<T>());
			if(lambda==0)
				return(0);
			//would this be more correct?
			//return(dataCount==0?0:-std::numeric_limits<T>::max());
			T sum(lambda);
			sum+=lgamma(dataCount+1);
			return(dataCount*log(lambda)-sum);
		}
	};

	struct dimaLikelihood{
		template<typename T>
		T operator()(double dataCount, likelihood::detail::ExpectationProperties<T> const & expect) const {
		    std::vector<T> const & simulationWeights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
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
					//std::cout << "	N-R stepped to " << std::setprecision(16) << lambda_i << std::setprecision(6) << std::endl;
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

				//std::cout << "	labmda_i now " << std::setprecision(16) << lambda_i << std::setprecision(6) << std::endl;
				//std::cout << "	 R_i_last = " << R_i_last << std::endl;
				//assert(!std::isinf((double)R_i_last));
				//std::cout << "	 R_i = " << R_i << std::endl;
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
//				//std::cout << "	w: " << w_ij << ' ' << llh << std::endl;
//			}
			for(auto t_ij : t_i)
				llh+=log(t_ij)-log(t_ij+lambda_i);
			if(dataCount)
				llh+=dataCount*-log(1-lambda_i);
			//std::cout << "	dataCount=" << dataCount << " log(1-lambda_i)=" << log(1-lambda_i) << '\n';
			//std::cout << "   llh=" << llh << '\n';
			return(llh);
		}
	};

	struct chi2Likelihood{
		template<typename T>
		T operator()(double dataCount, likelihood::detail::ExpectationProperties<T> const & expect) const {
		    std::vector<T> const & expectationWeights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
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
		T operator()(double dataCount, likelihood::detail::ExpectationProperties<T> const & expect) const {
		    std::vector<T> const & simulationWeights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
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
		T operator()(double dataCount, likelihood::detail::ExpectationProperties<T> const & expect) const {
		    std::vector<T> const & simulationWeights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
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

    struct glusenkampEqualWeightsLikelihood {
        template<typename T>
		T operator()(double k, likelihood::detail::ExpectationProperties<T> const & expect) const {
            std::vector<T> const & weights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
			T w_sum = std::accumulate(weights.begin(), weights.end(), T(0), std::plus<T>());
			T w2_sum = std::accumulate(weightUncertainties.begin(), weightUncertainties.end(), T(0), std::plus<T>());
            double kmc = weights.size();
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
            T L = likelihood::detail::gammaPriorPoissonLikelihood()(k, alpha, beta);

            return L;
        }
    };

    struct SAYMMSELikelihood {
        template<typename T>
		T operator()(double k, likelihood::detail::ExpectationProperties<T> const & expect) const {
            std::vector<T> const & weights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
			T w_sum = std::accumulate(weights.begin(), weights.end(), T(0), std::plus<T>());
			T w2_sum = std::accumulate(weightUncertainties.begin(), weightUncertainties.end(), T(0), std::plus<T>());
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

            T alpha = pow(w_sum, T(2))/w2_sum;
            T beta = w_sum/w2_sum;
            T L = likelihood::detail::gammaPriorPoissonLikelihood()(k, alpha, beta);

            return L;
        }
    };

    struct SAYMAPLikelihood {
        template<typename T>
		T operator()(double k, likelihood::detail::ExpectationProperties<T> const & expect) const {
            std::vector<T> const & weights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
			T w_sum = std::accumulate(weights.begin(), weights.end(), T(0), std::plus<T>());
			T w2_sum = std::accumulate(weightUncertainties.begin(), weightUncertainties.end(), T(0), std::plus<T>());
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

            const T & mu = w_sum;
            T mu2 = pow(mu, T(2));
            const T & sigma2 = w2_sum;

            T beta = (mu + sqrt(mu2+sigma2*4.0))/(sigma2*2);
            T alpha = (mu*sqrt(mu2+sigma2*4.0)/sigma2 + mu2/sigma2 + 2.0) / 2.0;
            T L = likelihood::detail::gammaPriorPoissonLikelihood()(k, alpha, beta);

            return L;
        }
    };

    struct SAYLikelihood {
        const bool is_map;
        SAYLikelihood():is_map(false) {}
        SAYLikelihood(bool is_map):is_map(is_map) {}
        SAYMAPLikelihood MAP;
        SAYMMSELikelihood MMSE;
        template<typename T>
		T operator()(double k, likelihood::detail::ExpectationProperties<T> const & expect) const {
            if(is_map)
                return MAP(k, expect);
            else
                return MMSE(k, expect);
        }
    };

	///Schneider Arguelles Yuan (empirical Bayes) likelihood
	///Developed by Austin Schneider, Carlos Arguelles, and Tianlu Yuan
	///*citation needed
	struct SAYLikelihood {
		template<typename T>
		T operator()(double k, likelihood::detail::ExpectationProperties<T> const & expect) const {
		    std::vector<T> const & weights = expect.weights;
            std::vector<T> const & weightUncertainties = expect.weightUncertainties;
			T w_sum = std::accumulate(weights.begin(), weights.end(), T(0), std::plus<T>());
			T w2_sum = std::accumulate(weightUncertainties.begin(), weightUncertainties.end(), T(0), std::plus<T>());

			if(w_sum <= 0 || w2_sum < 0) {
				//return(k==0?0:-std::numeric_limits<T>::max());
				return T(0);
			}

			if(w2_sum == 0) {
				return poissonLikelihood()(k, expect);
			}

			T alpha = (w_sum*w_sum)/w2_sum;
			T beta = w_sum/w2_sum;
			T L = likelihood::detail::gammaPriorPoissonLikelihood()(k, alpha, beta);

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
			static likelihood::entryStoringBin<T> neutral(){ return(likelihood::entryStoringBin<T>()); }
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

	// fundamental wrapper object for setting up likelihood fits for models by comparing observed data
	// events to simulated events
	template<typename Event, typename HistogramsType, typename DataWeighterConstructor, typename WCollection, typename CPrior, typename LFunc, int MaxDerivativeDimension=-1>
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
		DataWeighterConstructor dataWeighterC;

		//object used to construct the weighter for expected events from a set of model parameters
		//must be callable with an std::vector<double> containing the model parameters
		WCollection simWeightC;

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
						  const DataWeighterConstructor& dataWeighterC,
						  const WCollection& simWeightC,
						  const LFunc& likelihoodFunction,
						  std::vector<double> parameterSeeds,
						  unsigned int evaluationThreadCount=1):
		observation(observation),
		simulations(simulations),
		continuousPrior(continuousPrior),
		discretePrior(discretePrior),
		dataWeighterC(dataWeighterC),
		simWeightC(simWeightC),
		likelihoodFunction(likelihoodFunction),
		parameterSeeds(parameterSeeds),
		lastBestDiscreteIndex(0),
		evaluationThreadCount(evaluationThreadCount)
		{}

		template<typename AltLFunc>
		LikelihoodProblem<Event, HistogramsType, DataWeighterConstructor, WCollection, CPrior, AltLFunc, MaxDerivativeDimension>
		makeAlternateLikelihood(const AltLFunc& altlikelihoodFunction) const{
			using result_type=LikelihoodProblem<Event, HistogramsType, DataWeighterConstructor, WCollection, CPrior, AltLFunc, MaxDerivativeDimension>;
			return(result_type(observation,simulations,continuousPrior,discretePrior,dataWeightC,simWeightC,altlikelihoodFunction,parameterSeeds,evaluationThreadCount));
		}

		std::vector<double> getSeed() const{
			return(parameterSeeds);
		}

		void setObservation(const HistogramsType& newObs){
			observation=newObs;
		}

		void setEvaluationThreadCount(size_t count){
			if(count==0)
				throw std::runtime_error("At least one evaluation thread is required");
			evaluationThreadCount=count;
		}

		WCollection& getSimulationWeighterCollection(){ return(simWeightC); }
		const WCollection& getSimulationWeighterCollection() const{ return(simWeightC); }

		//evaluate the (non-prior) contribution to the likelihood from one observation,expectation histogram pair
		template<typename DataType, typename DataWeighter typename HistogramType>
		void evaluateLikelihoodCore(const HistogramType& observation, std::function<likelihood::detail::ExpectationProperties<DataType>(const std::vector<Event>&)>& weighter, DataWeighter& dweighter,
									const HistogramType& simulation,
									std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions) const{
			const auto& likelihoodFunction=this->likelihoodFunction;
			auto dataWeightAccumulator=[this](double t, const Event& e){ return(t+dweighter(e)); };

			auto likelihoodContribution=[dataWeightAccumulator,&weighter,&simulation,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
				entryStoringBin<Event> obs=*it;

				double observationAmount=std::accumulate(obs.begin(),obs.end(),0.0,dataWeightAccumulator);

				likelihood::detail::ExpectationProperties<DataType> expect;
				typename HistogramType::const_iterator expIt=simulation.findBinIterator(it);
				if(expIt!=simulation.end())
					expect = weighter(((entryStoringBin<Event>)*expIt).entries());

				auto contribution=likelihoodFunction(observationAmount,expect);
				/*{
					std::lock_guard<std::mutex> lck(printMtx);
					DataType expectationAmount=std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0));
					std::cout << "	 obs coords:";
					for(unsigned int i=0; i<3; i++)
						std::cout << ' ' << it.getBinEdge(i);
					std::cout << ' ' << contribution
					<< ' ' << observationAmount << ' ' << expectationAmount;
					std::cout << " [";
					for(const auto& w : expectationWeights)
						std::cout << ' ' << w;
					std::cout << ']' << std::endl;
				}*/
				return(contribution);
			};
			auto likelihoodContributionNoObs=[&weighter,&observation,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
				auto obsIt=observation.findBinIterator(it);

				//only proceed if this bin does not exist in the observation
				if(obsIt==observation.end()){
					likelihood::detail::ExpectationProperties<DataType> expect = weighter(((entryStoringBin<Event>)*it).entries());

					auto contribution=likelihoodFunction(0,expect);
					/*{
						std::lock_guard<std::mutex> lck(printMtx);
						DataType expectationAmount=std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0));
						std::cout << "	 exp coords:";
						for(unsigned int i=0; i<DataDimension; i++)
							std::cout << ' ' << it.getBinEdge(i);
						std::cout << ' ' << contribution
						<< ' ' << 0.0 << ' ' << expectationAmount;
						std::cout << " [";
						for(const auto& w : expectationWeights)
							std::cout << ' ' << w;
						std::cout << ']' << std::endl;
					}*/
					return(contribution);
				}
				/*{
					std::lock_guard<std::mutex> lck(printMtx);
					std::cout << "	 exp coords:";
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

//		  //evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
//		  //[recursive/iterative case]
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
//		  //evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
//		  //[base case]
//		template<typename DataType, typename SimulationWeighter>
//		void evaluateLikelihoodIterate<0,DataType,SimulationWeighter>(const HistogramsType& observation, SimulationWeighter weighter, const HistogramsType& simulation,
//																 std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
//			//do nothing
//		}

		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
		//[recursive/iterative case]
		template<unsigned int Counter, typename Likelihood, typename DataType, typename SimulationWeighter, typename DataWeighter>
		struct evaluateLikelihoodIterator{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, DataWeighter& dweighter,
							const HistogramsType& simulation,
							std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
				//evaluate for this pair
				constexpr unsigned int idx=std::tuple_size<HistogramsType>::value - Counter;
				like.evaluateLikelihoodCore(std::get<idx>(observation), weighter, dweighter, std::get<idx>(simulation),
									   printMtx, pool, contributions);
				//evaluate for the next pair
				evaluateLikelihoodIterator<Counter-1, Likelihood, DataType, SimulationWeighter>{}
					(like, observation, weighter, dweighter, simulation, printMtx, pool, contributions);
			}
		};
		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
		//[base case]
		template<typename Likelihood, typename DataType, typename SimulationWeighter, typename DataWeighter>
		struct evaluateLikelihoodIterator<0,Likelihood,DataType,SimulationWeighter,DataWeighter>{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, DataWeighter& dweighter, const HistogramsType& simulation,
							std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
				//do nothing
			}
		};

		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
		//[top level]
		template<typename DataType, typename SimulationWeighter, typename DataWeighter>
		DataType evaluateLikelihood(SimulationWeighter weighter, DataWeighter dweighter, const HistogramsType& simulation) const{
			std::mutex printMtx;
			ThreadPool pool(evaluationThreadCount);
			std::vector<std::future<DataType>> contributions;
			const HistogramsType& observation=this->observation;

			//Iterate over the observeation and expectation tuples,
			//firing off tasks for computing the per-bin likelihoods.
			//Futures for these results will accumulate in contributions.
			evaluateLikelihoodIterator<std::tuple_size<HistogramsType>::value,decltype(*this),DataType,SimulationWeighter,DataWeighter>{}
			(*this, observation, weighter, dweighter, simulation, printMtx, pool, contributions);

			//Iterate over the futures collecting the results and combining them when they are availiable.
			DataType llh(0.0);
			for(auto& future : contributions){
				auto contribution=future.get();
				llh+=contribution;
			}
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
//							  std::cout << "Bad weight: " << expectationWeights.back() << "\nEvent:\n" << e << std::endl;
//							//std::cout << e.energy << ' ' << e.year << ' ' << expectationWeights.back() << std::endl;
//						}
//						//std::cout << "	" << expectationWeights.back() << std::endl;
//					}
//				}
//
//				auto contribution=likelihoodFunction(observationAmount,expectationWeights);
//				/*{
//					std::lock_guard<std::mutex> lck(printMtx);
//					DataType expectationAmount=std::accumulate(expectationWeights.begin(),expectationWeights.end(),DataType(0));
//					std::cout << "	 obs coords:";
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
//						std::cout << "	 exp coords:";
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
//					std::cout << "	 exp coords:";
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

		template<typename DataType, typename DataWeighter, typename DataHistogramType, typename ResultHistogramType>
		void evaluateLikelihoodContributionsCore(const DataHistogramType& observation, std::function<likelihood::detail::ExpectationProperties<DataType>(const std::vector<Event>&)>& weighter, DataWeighter dweighter,
				const DataHistogramType& simulation, ResultHistogramType& result) const{
			const auto& likelihoodFunction=this->likelihoodFunction;
			auto dataWeightAccumulator=[&dweighter](double t, const Event& e){ return(t+dweighter(e)); };

			using HistDimExt = phys_tools::histograms::detail::dimensionExtractor<phys_tools::histograms::histogram,DataHistogramType::dimensions,typename DataHistogramType::dataType>;
			const unsigned int histDimensions=HistDimExt::getDim(&observation);
			result.setUseContentScaling(false);

			auto likelihoodContribution=[dataWeightAccumulator,&weighter,&simulation,&likelihoodFunction,&result,histDimensions](typename DataHistogramType::const_iterator it)->void{
				entryStoringBin<Event> obs=*it;

				double observationAmount=std::accumulate(obs.begin(),obs.end(),0.0,dataWeightAccumulator);

				likelihood::detail::ExpectationProperties<DataType> expect;

				typename DataHistogramType::const_iterator expIt=simulation.findBinIterator(it);
				if(expIt!=simulation.end())
					expect = weighter(((entryStoringBin<Event>)*expIt).entries());
				DataType contribution=likelihoodFunction(observationAmount,expect);

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
					likelihood::detail::ExpectationProperties<DataType> expect = weighter(((entryStoringBin<Event>)*it).entries());
					contribution=likelihoodFunction(0,expect);
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
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter,
							const HistogramsType& simulation,
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
			if(std::isnan(cPrior)){
				//std::cout << "Out-of-bounds prior prob" << std::endl;
				return(-std::numeric_limits<DataType>::max());
			}
			//std::cout << "(prior=" << cPrior << ") ";

			//compute the actual llh for every necessary discrete nuisance index
			//(of which there will only be one unless we've been called by a minimizer)
			//and keep the best value
			DataType bestLLH=-std::numeric_limits<DataType>::max();
			for(size_t dn=minDN; dn<maxDN; dn++){
				DataType llh=(includePriors?cPrior+discretePrior[dn]:0);
				params.back()=dn;
				llh+=evaluateLikelihood<DataType>(simWeightC.template operator()<Event, DataType>(continuousParams),dataWeightC(continuousParams),simulations[dn]);

				if(llh>bestLLH){
					bestLLH=llh;
					lastBestDiscreteIndex=dn;
				}
			}
			//std::cout << " llh: " << bestLLH << std::endl;
			return(bestLLH);
		}

//		  ///\return a pair consitsting of the LLH per bin in the form of a histogram (with the same dimensions and axes as the data) and the log probability due to the priors
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

	///Represents the result of evaluating a likelihood at one point in the
	///space over which it is defined
	struct likelihoodPoint{
		///the combination of parameters for which the likelihood was evaluated
		std::vector<double> params;
		///the (log) likelihood value
		double likelihood;
	};

	//helper function for constructing a LikelihoodProblem without having to write out all of the template parameters
	//Note that the order of the template parameters has been shuffled, to put those which cannot be deduced from the
	//arguments (Event, DataDimension, and MaxDerivativeDimension) first so that they can be specified while the rest
	//are left to be deduced, while the order of the function arguments is the same as for the LikelihoodProblem constructor
	template<typename Event, int MaxDerivativeDimension=-1, typename HistogramsType, typename DataWeighter,
			 typename... Weighters,
			 typename CPrior, typename LFunc,
			 typename LikelihoodType=LikelihoodProblem<Event,HistogramsType,DataWeighterConstructor,likelihood::detail::WeighterCollection<Weighters...>,CPrior,LFunc,MaxDerivativeDimension>>
	LikelihoodType makeLikelihoodProblem(HistogramsType observation, const std::vector<HistogramsType>& simulations,
										 CPrior continuousPrior, std::vector<double> discretePrior,
										 const DataWeighterConstructor& dataWeighterConstructor,
										 const likelihood::detail::WeighterCollection<Weighters...>& simWeightC,
										 const LFunc& likelihoodFunction,
										 std::vector<double> parameterSeeds,
										 unsigned int evaluationThreadCount=1) {
		return(LikelihoodType(observation,simulations,
							  continuousPrior,discretePrior,
							  dataWeighterConstructor,simWeightC,likelihoodFunction,
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
	Container
	generateSample(const std::vector<double>& weights, const Container& data, double quantity, RNG& eng){
		assert(!data.empty());
		assert(weights.size()==data.size());
		std::discrete_distribution<size_t> dist(weights.begin(),weights.end());
		//decide how many events to sample
		size_t s_quantity = std::poisson_distribution<size_t>(quantity)(eng);
		//generate the sample
		Container sample;
		//sample.reserve(s_quantity);
		for(unsigned int i=0; i<s_quantity; i++){
			size_t idx=dist(eng);
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

	struct GaussianPrior{
	private:
		double mean;
		double stddev;
		double norm;
	public:
		GaussianPrior(double mean, double stddev):
		mean(mean),stddev(stddev),
		norm(boost::math::constants::one_div_root_two_pi<double>()/stddev){}

		template<typename DataType>
		DataType operator()(DataType x) const{
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

	///A structure containing hints used by findParameterUncertainties to guide
	///its bracketing and root finding of confidence interval endpoints
	struct UncertaintySearchPlan{
		///functions which given the best fit value of the paramter guess a value
		///which defines the scale for searching for the confidence interval's left
		///and right endpoints, respectively
		std::function<double(double)> guessLow, guessHigh;
		///The precision with which the boundaries of the confidence interval
		///should be computed
		double precision;
	};

	///Quickly find confidence intervals on fit parameters by assuming Wilks Theorem.
	///Parameters for which this function is to be used are expected to have
	///UncertaintySearchPlan associated with them in the ParameterSet, under the
	///name "UncertaintySearchPlan".
	///\param prob the likelihood problem for which to compute uncertainties
	///\param ps the parameter information for the likelihood problem
	///\param fit a function which fits for parameters which maximize the
	///           likelihood, subject to parameters being held fixed.
	///\param bestFit the maximum likelihood point (with any constraints)
	///\param indicesToCompute the indices of fit paramters for which confidence
	///       intervals should be computed, or an empty list if all free
	///       parameters should be treated
	///\param confidence the confidence level for which to compute intervals
	///\return a vector of pairs, where each pair contains the lower and upper
	///        bound of the confidence interval for one parameter. For parameters
	///        held fixed or not scanned over, both entries in the pair will be
	///        the best fit value for that parameter.
	template<typename ProblemType>
	std::vector<std::pair<double,double>>
	findParameterUncertainties(ProblemType& prob, const ParameterSet& ps,
							   std::function<likelihoodPoint(const ProblemType&,const ParameterSet&)> fit,
							   const likelihoodPoint& bestFit,
	                           const std::vector<unsigned int>& indicesToCompute={},
	                           double confidence=0.682689492137){
		//somewhat brute force approach due to gcc 4.8 bugs
		auto chi2_cdf=[](double x){
			const double norm=boost::math::constants::root_pi<double>()/tgamma(0.5);
			return(norm*erf(sqrt(x/2)));
		};
		auto br=bracketRoot(chi2_cdf,confidence,0,1,0);
		const double targetLLHDifference=findRoot(chi2_cdf,confidence,br.xin,br.xout,br.fin,br.fout,1e-14)/2;
		//std::cout << "target LLH difference for a confidence of " << confidence << " is " << targetLLHDifference << std::endl;

		//if(!quiet)
		//	std::cout << "Finding one dimensional parameter uncertainites at confidence " << confidence << std::endl;

		auto makeScanFunc=[&](unsigned int idx){
			auto seed=bestFit.params;
			ParameterSet cps(ps);
			cps.fixParameter(idx);
			return([&,seed,idx,cps](double paramValue) mutable{
				cps.setParameterValue(idx,paramValue);
				likelihoodPoint fr = fit(prob,cps);
				return(fr.likelihood);
			});
		};

		auto findParameterError=[&](unsigned int idx, double guessLow, double guessHigh, double precision, double limitLow, double limitHigh){
			const bool verbose=false;
			auto scanFunc=makeScanFunc(idx);
			double middle=bestFit.params[idx];
			if(guessLow<limitLow)
				guessLow=limitLow;
			if(guessHigh>limitHigh)
				guessHigh=limitHigh;

			//left side
			auto bracket=bracketRoot(scanFunc, bestFit.likelihood+targetLLHDifference, middle, guessLow, bestFit.likelihood, limitLow, verbose);
			double lower;
			if(bracket.fout<bestFit.likelihood+targetLLHDifference)
				lower=bracket.xout;
			else
				lower=findRoot(scanFunc, bestFit.likelihood+targetLLHDifference, bracket.xout, bracket.xin, bracket.fout, bracket.fin, precision, verbose);

			//right side
			bracket=bracketRoot(scanFunc, bestFit.likelihood+targetLLHDifference, middle, guessHigh, bestFit.likelihood, limitHigh, verbose);
			double upper;
			if(bracket.fout<bestFit.likelihood+targetLLHDifference)
				upper=bracket.xout;
			else
				upper=findRoot(scanFunc, bestFit.likelihood+targetLLHDifference, bracket.xin, bracket.xout, bracket.fin, bracket.fout, precision, verbose);

			//std::cout << "Error interval for parameter " << idx << " is [" << lower << ',' << upper << ']' << std::endl;
			return(std::make_pair(lower,upper));
		};

		//check whether a given parameter should be scanned
		auto shouldScan=[&indicesToCompute](unsigned int idx){
			return(indicesToCompute.empty() || std::find(indicesToCompute.begin(),indicesToCompute.end(),idx)!=indicesToCompute.end());
		};

		std::vector<std::pair<double,double>> results;
		for(unsigned int i=0; i<ps.numberOfParameters(); i++){
			//skip parameters we're not interested in right now
			//skip paramters which are being held fixed
			//skip parameters which have no search plan attached to them
			if(!shouldScan(i) || ps.isFixed(i) || !ps.parameterHasProperty(i,"UncertaintySearchPlan")){
				results.push_back(std::make_pair(bestFit.params[i],bestFit.params[i]));
				continue;
			}

			const auto& plan=ps.getParameterProperty<UncertaintySearchPlan>(i,"UncertaintySearchPlan");
			double guessLow=plan.guessLow(bestFit.params[i]);
			double guessHigh=plan.guessHigh(bestFit.params[i]);
			auto result =
			findParameterError(i, guessLow, guessHigh, plan.precision,
							   ps.getParameterLowerLimit(i),
							   ps.getParameterUpperLimit(i));
			results.push_back(result);
		}
		return(results);
	}

} //namespace likelihood
} //namespace phys_tools

#endif
