#ifndef LF_LIKELIHOOD_H
#define LF_LIKELIHOOD_H

#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <deque>

#include <boost/math/constants/constants.hpp>

#include <boost/mpl/if.hpp>

#include "../lbfgsb/lbfgsb.h"

#include "../histogram.h"
#include "../brent.h"

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
	
	//Frequently used likelihood functions
	
	struct poissonLikelihood{
		template <typename T>
		T operator()(double dataCount, const std::vector<T>& simulationWeights, std::vector<unsigned int>& categories) const{
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
		T operator()(double dataCount, const std::vector<T>& simulationWeights, std::vector<unsigned int>& categories) const{
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
		T operator()(unsigned int dataCount, const std::vector<T>& expectationWeights, std::vector<unsigned int>& categories) const{
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
		T operator()(double dataCount, const std::vector<T>& simulationWeights, std::vector<unsigned int>& categories) const{
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
		T operator()(double dataCount, const std::vector<T>& simulationWeights, std::vector<unsigned int>& categories) const{
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
    DataType compute_barlow_ti(double di, const std::vector<unsigned int>& ai, std::vector<DataType>& wi) {
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
        DataType operator()(double di, const std::vector<unsigned int>& ai, std::vector<DataType>& wi) {
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
        result_type operator()(double di, const std::vector<unsigned int>& ai, std::vector<result_type>& dwi) {
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

            return compute_barlow_LLH<T>()(dataCount, ai, wi);
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
				return(std::get<size-index>(priors)(params[size-index])*priorEvaluator<DataType,index-1>()(priors,params));
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
	template<typename Event, typename HistogramsType, typename DataWeighter, typename SimulationWeighterConstructor, typename CPrior, typename LFunc, int MaxDerivativeDimension=-1>
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
		DataWeighter dataWeighter;
		
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
						  const DataWeighter& dataWeighter,
						  const SimulationWeighterConstructor& simWeightC,
						  const LFunc& likelihoodFunction,
						  std::vector<double> parameterSeeds,
						  unsigned int evaluationThreadCount=1):
		observation(observation),
		simulations(simulations),
		continuousPrior(continuousPrior),
		discretePrior(discretePrior),
		dataWeighter(dataWeighter),
		simWeightC(simWeightC),
		likelihoodFunction(likelihoodFunction),
		parameterSeeds(parameterSeeds),
		lastBestDiscreteIndex(0),
		evaluationThreadCount(evaluationThreadCount)
		{}
		
		template<typename AltLFunc>
		LikelihoodProblem<Event, HistogramsType, DataWeighter, SimulationWeighterConstructor, CPrior, AltLFunc, MaxDerivativeDimension>
		makeAlternateLikelihood(const AltLFunc& altlikelihoodFunction) const{
			using result_type=LikelihoodProblem<Event, HistogramsType, DataWeighter, SimulationWeighterConstructor, CPrior, AltLFunc, MaxDerivativeDimension>;
			return(result_type(observation,simulations,continuousPrior,discretePrior,simWeightC,altlikelihoodFunction,parameterSeeds,evaluationThreadCount));
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
		
		SimulationWeighterConstructor& getSimulationWeighterConstructor(){ return(simWeightC); }
		const SimulationWeighterConstructor& getSimulationWeighterConstructor() const{ return(simWeightC); }
		
		//evaluate the (non-prior) contribution to the likelihood from one observation,expectation histogram pair
		template<typename DataType, typename SimulationWeighter, typename HistogramType>
		void evaluateLikelihoodCore(const HistogramType& observation, SimulationWeighter& weighter, const HistogramType& simulation,
								    std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions) const{
			const auto& likelihoodFunction=this->likelihoodFunction;
			auto dataWeightAccumulator=[this](double t, const Event& e){ return(t+this->dataWeighter(e)); };
			
			auto likelihoodContribution=[dataWeightAccumulator,&weighter,&simulation,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
				entryStoringBin<Event> obs=*it;
				
				double observationAmount=std::accumulate(obs.begin(),obs.end(),0.0,dataWeightAccumulator);
				
				std::vector<DataType> expectationWeights;
				std::vector<unsigned int> categoryWeights;
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
					categoryWeights.reserve(((entryStoringBin<Event>)*expIt).size());
					for(const RawEvent& e : ((entryStoringBin<Event>)*expIt)){
						expectationWeights.push_back(weighter(e));
						categoryWeights.push_back(e.category);
						if(std::isnan(expectationWeights.back()) || expectationWeights.back()<0.0){
							std::lock_guard<std::mutex> lck(printMtx);
							std::cout << "Bad weight: " << expectationWeights.back() << "\nEvent:\n" << e << std::endl;
							//std::cout << e.energy << ' ' << e.year << ' ' << expectationWeights.back() << std::endl;
						}
						//std::cout << "    " << expectationWeights.back() << std::endl;
					}
				}

                std::vector<unsigned int> categories;
                if(categoryWeights.size()>0) {
                    unsigned int category = categoryWeights[0];
                    std::cout << "categories = [" << category;
                    unsigned int count = 1;

                    for(auto it=categoryWeights.begin()+1; it!=categoryWeights.end(); ++it) {
                        if(*it != category) {
                            category = *it;
                            std::cout << ", " << category;
                            categories.push_back(count);
                            count = 0;
                        }
                        ++count;
                    }
                    categories.push_back(count);
                    std::cout << "]" << std::endl;
                }
				
				auto contribution=likelihoodFunction(observationAmount,expectationWeights,categories);
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
				return(contribution);
			};
			auto likelihoodContributionNoObs=[&weighter,&observation,&likelihoodFunction,&printMtx](typename HistogramType::const_iterator it)->DataType{
				auto obsIt=observation.findBinIterator(it);
				
				//only proceed if this bin does not exist in the observation
				if(obsIt==observation.end()){
					std::vector<DataType> expectationWeights;
					std::vector<unsigned int> categoryWeights;
					const std::vector<Event>& exp=((entryStoringBin<Event>)*it).entries();
					expectationWeights.reserve(((entryStoringBin<Event>)*it).size());
					categoryWeights.reserve(((entryStoringBin<Event>)*it).size());
					for(const RawEvent& e : ((entryStoringBin<Event>)*it)) {
						expectationWeights.push_back(weighter(e));
						categoryWeights.push_back(e.category);
                    }
					
                    std::vector<unsigned int> categories;
                    if(categoryWeights.size()>0) {
                        unsigned int category = categoryWeights[0];
                        std::cout << "categories = [" << category;
                        unsigned int count = 1;

                        for(auto it=categoryWeights.begin()+1; it!=categoryWeights.end(); ++it) {
                            if(*it != category) {
                                category = *it;
                                std::cout << ", " << category;
                                categories.push_back(count);
                                count = 0;
                            }
                            ++count;
                        }
                        categories.push_back(count);
                        std::cout << "]" << std::endl;
                    }
	
					auto contribution=likelihoodFunction(0,expectationWeights,categories);
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
		template<unsigned int Counter, typename Likelihood, typename DataType, typename SimulationWeighter>
		struct evaluateLikelihoodIterator{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, const HistogramsType& simulation,
							std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
				//evaluate for this pair
				constexpr unsigned int idx=std::tuple_size<HistogramsType>::value - Counter;
				like.evaluateLikelihoodCore(std::get<idx>(observation), weighter, std::get<idx>(simulation),
									   printMtx, pool, contributions);
				//evaluate for the next pair
				evaluateLikelihoodIterator<Counter-1, Likelihood, DataType, SimulationWeighter>{}
					(like, observation, weighter, simulation, printMtx, pool, contributions);
			}
		};
		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
	    //[base case]
		template<typename Likelihood, typename DataType, typename SimulationWeighter>
		struct evaluateLikelihoodIterator<0,Likelihood,DataType,SimulationWeighter>{
			void operator()(const Likelihood& like, const HistogramsType& observation, SimulationWeighter& weighter, const HistogramsType& simulation,
							std::mutex& printMtx, ThreadPool& pool, std::vector<std::future<DataType>>& contributions){
				//do nothing
			}
		};
		
		//evaluate the (non-prior) contribution to the likelihood from one pair of observation,expectation histogram tuples
		//[top level]
		template<typename DataType, typename SimulationWeighter>
		DataType evaluateLikelihood(SimulationWeighter weighter, const HistogramsType& simulation) const{
			std::mutex printMtx;
	        ThreadPool pool(evaluationThreadCount);
			std::vector<std::future<DataType>> contributions;
			const HistogramsType& observation=this->observation;
			
			//Iterate over the observeation and expectation tuples,
			//firing off tasks for computing the per-bin likelihoods.
			//Futures for these results will accumulate in contributions.
			evaluateLikelihoodIterator<std::tuple_size<HistogramsType>::value,decltype(*this),DataType,SimulationWeighter>{}
			(*this, observation, weighter, simulation, printMtx, pool, contributions);
			
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
        
		template<typename DataType, typename SimulationWeighter, typename DataHistogramType, typename ResultHistogramType>
		void evaluateLikelihoodContributionsCore(const DataHistogramType& observation, SimulationWeighter& weighter, const DataHistogramType& simulation, ResultHistogramType& result) const{
			const auto& likelihoodFunction=this->likelihoodFunction;
			auto dataWeightAccumulator=[this](double t, const Event& e){ return(t+this->dataWeighter(e)); };
			
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
			if(cPrior<=0 || std::isnan(cPrior)){
				//std::cout << "Out-of-bounds prior prob" << std::endl;
				return(-std::numeric_limits<DataType>::max());
			}
			//std::cout << "(prior=" << cPrior << ") ";
			cPrior=log(cPrior); //convert to log(P) for rest of calculation
			
			//compute the actual llh for every necessary discrete nuisance index
			//(of which there will only be one unless we've been called by a minimizer)
			//and keep the best value
			DataType bestLLH=-std::numeric_limits<DataType>::max();
			for(size_t dn=minDN; dn<maxDN; dn++){
				DataType llh=(includePriors?cPrior+log(discretePrior[dn]):0);
				params.back()=dn;
				llh+=evaluateLikelihood<DataType>(simWeightC(continuousParams),simulations[dn]);
				
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
	template<typename Event, int MaxDerivativeDimension=-1, typename HistogramsType, typename DataWeighter, typename SimulationWeighterConstructor, typename CPrior, typename LFunc,
	         typename LikelihoodType=LikelihoodProblem<Event,HistogramsType,DataWeighter,SimulationWeighterConstructor,CPrior,LFunc,MaxDerivativeDimension>>
	LikelihoodType makeLikelihoodProblem(HistogramsType observation, const std::vector<HistogramsType>& simulations,
                                         CPrior continuousPrior, std::vector<double> discretePrior,
                                         const DataWeighter& dataWeighter,
                                         const SimulationWeighterConstructor& simWeightC,
                                         const LFunc& likelihoodFunction,
                                         std::vector<double> parameterSeeds,
                                         unsigned int evaluationThreadCount=1){
		return(LikelihoodType(observation,simulations,
							  continuousPrior,discretePrior,
							  dataWeighter,simWeightC,likelihoodFunction,
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
				return(DataType(0.0));
			return(DataType(1.0));
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
			return(norm*exp(-z*z/2));
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
			return(limits(x)*prior(x));
		}
	};
	
} //namespace likelihood
} //namespace phys_tools

#endif
