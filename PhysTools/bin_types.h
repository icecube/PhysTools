///\file bin_types.h
///This file defines various useful types for histograms which are alternatives to storing simple numbers

#ifndef PHYSTOOLS_BIN_TYPES_H
#define PHYSTOOLS_BIN_TYPES_H

#include "PhysTools/hdf5_serialization.h"
#include "PhysTools/histogram.h"

namespace phys_tools{
	namespace histograms{
		using namespace phys_tools::detail;
		
		///A simple bin type which accumulated not only the total of the entries added,
		///but also the number of entries which have been made
		template<typename DataType=double>
		struct entryCounter{
		private:
			DataType amount;
			unsigned int count;
		public:
			entryCounter():amount(0),count(0){}
			entryCounter(DataType v):amount(v),count(1){}
			
			entryCounter& operator+=(DataType x){
				amount+=x;
				count++;
				return(*this);
			}
			
			entryCounter& operator+=(const entryCounter& other){
				amount+=other.amount;
				count+=other.count;
				return(*this);
			}
			
			entryCounter& operator*=(DataType a){
				amount*=a;
				return(*this);
			}
			
			entryCounter& operator/=(DataType a){
				amount/=a;
				return(*this);
			}
			
			const entryCounter operator+(const entryCounter& other) const{
				return(entryCounter(*this)+=other);
			}
			
			const entryCounter operator+(DataType a) const{
				return(entryCounter(*this)+=a);
			}
			
			const entryCounter operator*(DataType a) const{
				return(entryCounter(*this)*=a);
			}
			
			const entryCounter operator/(DataType a) const{
				return(entryCounter(*this)/=a);
			}
			
			operator DataType() const{ return(amount); }
			
			const DataType& value() const{ return(amount); }
			unsigned int entries() const{ return(count); }
			
			template<class Archive>
			void serialize(Archive & ar, const unsigned int version){
				ar & boost::serialization::make_nvp("amount",amount);
				ar & boost::serialization::make_nvp("count",count);
			}
		};
		
		namespace detail{
			template<typename T>
			class histogramTraits<entryCounter<T>>{
			public:
				struct amount{
					T value;
					amount(T v):value(v){}
					
					template <typename OtherType>
					amount(const OtherType& other):value(other.value){}
				};
				constexpr static bool enable_automatic_amount_handling=true;
				static void defaultData(entryCounter<T>* data, unsigned int count){}
				static entryCounter<T> unit(){ return(entryCounter<T>(1)); }
			};
		} //namespace detail
		
		//----------------------------------------
		
		///A bin type which keeps a running calculation of the mean and variance of the entries
		///which have been added.
		template<typename DataType=double>
		class meanVarTracker{
		public:
			meanVarTracker():count(0),m(0),t(0),w(0){}
			meanVarTracker(DataType v, DataType w=1):count(1),m(v),t(0),w(w){}
			
			void insert(DataType x, DataType weight=1){
				if(!count){
					m=x;
					w=weight;
				}
				else{
					DataType q=x-m;
					DataType temp=w+weight;
					DataType r=q*weight/temp;
					t+=weight*(x-m)*(x-m-r);
					m+=r;
					w=temp;
				}
				count++;
			}
			
			meanVarTracker& operator+=(const meanVarTracker& other){
				count += other.count;
				m = (m*w + other.w*other.m)/(w+other.w);
				t = (t*w + other.w*other.t)/(w+other.w); //TODO: verify that this expression is correct
				w += other.w;
				return(*this);
			}
			
			meanVarTracker& operator+=(DataType x){
				insert(x);
				return(*this);
			}
			
			meanVarTracker& operator*=(DataType a){
				m *= a;
				t *= a*a;
				w *= a;
				return(*this);
			}
			
			meanVarTracker& operator/=(DataType a){
				m /= a;
				t /= a*a;
				w /= a;
				return(*this);
			}
			
			const meanVarTracker operator+(const meanVarTracker& other) const{
				return(meanVarTracker(*this)+=other);
			}
			
			const meanVarTracker operator+(DataType x) const{
				return(meanVarTracker(*this)+=x);
			}
			
			const meanVarTracker operator*(DataType a) const{
				return(meanVarTracker(*this)*=a);
			}
			
			const meanVarTracker operator/(DataType a) const{
				return(meanVarTracker(*this)/=a);
			}
			
			unsigned int entries() const{ return(count); }
			DataType mean() const{ return(m); }
			///Get the _biased_ sample variance defined as
			// = (1/n) sum w_i*(x_i^2 - mu^2)
			// note that the overall factor is not 1/(n-1)
			DataType variance() const{
				return(count>1?t/w:0.0);
			}
			DataType standardDeviation() const{ return(sqrt(variance())); }
			double errorMin() const{ return(mean()-standardDeviation()); }
			double errorMax() const{ return(mean()+standardDeviation()); }
			
			template<class Archive>
			void serialize(Archive & ar, const unsigned int version){
				ar & boost::serialization::make_nvp("count",count);
				ar & boost::serialization::make_nvp("m",m);
				ar & boost::serialization::make_nvp("t",t);
				ar & boost::serialization::make_nvp("w",w);
			}
			
		private:
			unsigned int count;
			DataType m, t, w;
		};
		
		template<typename DataType>
		meanVarTracker<DataType> operator*(DataType a, const meanVarTracker<DataType>& mvt){
			return(mvt*a);
		}
		
		template<typename DataType>
		std::ostream& operator<<(std::ostream& os, const meanVarTracker<DataType>& t){
			return(os << t.mean() << '[' << t.errorMin() << ',' << t.errorMax() << ']');
		}
		
		namespace detail{
			template<typename T>
			class histogramTraits<meanVarTracker<T>>{
			public:
				struct amount{
					T value, weight;
					amount(T v):value(v),weight(1){}
					amount(T v, T w):value(v),weight(w){}
					
					template <typename OtherType>
					amount(const OtherType& other):value(other.value),weight(1){}
				};
				constexpr static bool enable_automatic_amount_handling=true;
				static void defaultData(meanVarTracker<T>* data, unsigned int count){}
				static meanVarTracker<T> unit(){ return(meanVarTracker<T>(1)); }
			};
			
			template<typename DataType>
			struct amountAdder<meanVarTracker<DataType>>{
				amountAdder(meanVarTracker<DataType>& v, const typename detail::histogramTraits<meanVarTracker<DataType>>::amount& a){
					v.insert(a.value,a.weight);
				}
				static void add(meanVarTracker<DataType>& v, const typename detail::histogramTraits<meanVarTracker<DataType>>::amount& a){
					v.insert(a.value,a.weight);
				}
			};
		} //namespace detail
		
		//----------------------------------------
		
		///A bin type which stores both a sum and its error for entries whose individual errors
		///are assumed to be the square roots of their values. The error on the sum is treated as
		///the square root of the sum of squares of the entry errors (that is, its variance is
		///the sum of the entry variances).
		template<typename DataType=double>
		struct sqErrorValue{
		public:
			DataType amount;
			DataType error;
		public:
			sqErrorValue():amount(0),error(0){}
			sqErrorValue(DataType v):amount(v),error(v*v){}
			///Construct directly with a value and a variance
			///\param v the value
			///\param e the variance
			sqErrorValue(DataType v, DataType e):amount(v),error(e){}
			sqErrorValue(const sqErrorValue<DataType>& other):amount(other.amount),error(other.error){}
			
			sqErrorValue& operator+=(DataType x){
				amount+=x;
				error+=x*x;
				return(*this);
			}
			
			sqErrorValue& operator+=(const sqErrorValue& other){
				amount+=other.amount;
				error+=other.error;
				return(*this);
			}
			
			sqErrorValue& operator*=(DataType a){
				amount*=a;
				error*=a*a;
				return(*this);
			}
			
			sqErrorValue& operator/=(DataType a){
				amount/=a;
				error/=a*a;
				return(*this);
			}
			
			const sqErrorValue operator+(const sqErrorValue& other) const{
				return(sqErrorValue(*this)+=other);
			}
			
			const sqErrorValue operator+(DataType a) const{
				return(sqErrorValue(*this)+=a);
			}
			
			const sqErrorValue operator*(DataType a) const{
				return(sqErrorValue(*this)*=a);
			}
			
			const sqErrorValue operator/(DataType a) const{
				return(sqErrorValue(*this)/=a);
			}
			
			operator DataType() const{ return(amount); }
			
			const DataType& value() const{ return(amount); }
			DataType variance() const{ return(error); }
			DataType standardDeviation() const{ return(sqrt(variance())); }
			double errorMin() const{ return(amount-standardDeviation()); }
			double errorMax() const{ return(amount+standardDeviation()); }
			
			template<class Archive>
			void serialize(Archive & ar, const unsigned int version){
				ar & boost::serialization::make_nvp("amount",amount);
				ar & boost::serialization::make_nvp("error",error);
			}
		};
		
		template<typename DataType>
		sqErrorValue<DataType> operator*(DataType a, const sqErrorValue<DataType>& wet){
			return(wet*a);
		}
		
		template<typename T>
		std::ostream& operator<<(std::ostream& os, const sqErrorValue<T>& t){
			return(os << t.amount << '[' << t.errorMin() << ',' << t.errorMax() << ']');
		}
		
		namespace detail{
			template<typename T>
			class histogramTraits<sqErrorValue<T>>{
			public:
				using amount=amount_type<sqErrorValue<T>>;
				constexpr static bool enable_automatic_amount_handling=true;
				static void defaultData(sqErrorValue<T>* data, unsigned int count){}
				static sqErrorValue<T> unit(){ return(sqErrorValue<T>(1)); }
			};
		}
		
		namespace detail{
			template<typename T>
			struct interval{
				T min,max;
				
				bool contains(const T& x) const{
					return(x>=min && x<=max);
				}
			};
			
			///Compute a Feldman-Cousins acceptance interval for a Poisson
			///process with optional background
			///\param mean the hypothetical true mean
			///\param confidence the confidence level at which to compute the interval
			///\param offset Poisson background mean
			interval<unsigned int> fcPoissonAcceptanceInterval(double mean, double confidence, double offset=0.0);
			
			///Compute a Feldman-Cousins confidence interval for a Poisson
			///process with optional background
			///\param obs the observed number of counts
			///\param confidence the confidence level at which to compute the interval
			///\param offset Poisson background mean
			///\param tol the fractional tolerance with which the endpoints of the interval should be determined
			interval<double> fcPoissonConfidenceInterval(unsigned int obs, double confidence, double offset=0.0, double tol=1e-4);
		}
		
		struct fcErrorValue{
		private:
			double value;
			mutable double errMin, errMax;
			mutable bool errorsUpToDate;
			
			void recomputeErrors() const{
				if(errorsUpToDate)
					return;
				if(value<100){
					detail::interval<double> errRange=detail::fcPoissonConfidenceInterval(value, .6827, 0.0, 1e-3);
					errMin=errRange.min;
					errMax=errRange.max;
				}
				else{
					double r=sqrt(value);
					errMin=value-r;
					errMax=value+r;
				}
				errorsUpToDate=true;
			}
		public:
			fcErrorValue():value(0),errorsUpToDate(false){}
			
			fcErrorValue(double d):value(d),errorsUpToDate(false){}
			
			fcErrorValue& operator+=(const double& d){
				value+=d;
				errorsUpToDate=false;
				return(*this);
			}
			fcErrorValue& operator+=(const fcErrorValue& other){
				value+=other.value;
				errorsUpToDate=false;
				return(*this);
			}
			fcErrorValue& operator*=(double scale){
				value*=scale;
				errorsUpToDate=false;
				return(*this);
			}
			fcErrorValue operator+(double a) const{
				return(fcErrorValue(value+a));
			}
			fcErrorValue operator+(const fcErrorValue& other) const{
				return(fcErrorValue(value+other.value));
			}
			fcErrorValue operator*(double scale) const{
				return(fcErrorValue(value*scale));
			}
			
			operator double() const{ return(value); }
			double errorMin() const{ recomputeErrors(); return(errMin); }
			double errorMax() const{ recomputeErrors(); return(errMax); }
			
			template<class Archive>
			void save(Archive & ar, const unsigned int version) const{
				ar << boost::serialization::make_nvp("value",value);
			}
			
			template<class Archive>
			void load(Archive & ar, const unsigned int version){
				ar >> boost::serialization::make_nvp("value",value);
				errorsUpToDate=false;
			}
			
			BOOST_SERIALIZATION_SPLIT_MEMBER();
		};
		
		fcErrorValue operator*(double scale, const fcErrorValue& v);
		std::ostream& operator<<(std::ostream& os, const fcErrorValue& v);
		
		namespace detail{
			template<>
			class histogramTraits<fcErrorValue>{
			public:
				using amount=amount_type<fcErrorValue>;
				constexpr static bool enable_automatic_amount_handling=true;
				static void defaultData(fcErrorValue* data, unsigned int count){/* Nothing to do */}
				static fcErrorValue unit(){ return 1.0; }
			};
		}
		
	} //namespace histograms
} //namespace phys_tools

#endif
