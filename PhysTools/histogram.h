#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <iosfwd>

#include "PhysTools/axis.h"
#include "PhysTools/detail/param_pack_traits.h"

namespace phys_tools{

	///Namespace for histogramming classes and functions
	namespace histograms{
		using namespace phys_tools::detail;
		
		namespace detail{
			//wrapper type used to tag amounts being added to histograms
			//as different from coordinates
			template<typename DataType>
			struct amount_type{
				DataType value;
				amount_type(DataType v):value(v){}
				
				template<typename OtherDataType>
				amount_type(const amount_type<OtherDataType>& a):value(a.value){}
			};
			
			//traits class for recognizing histograms::amount types
			template<typename T>
			struct is_amount : public std::false_type{};
			template<typename D>
			struct is_amount<amount_type<D>> : public std::true_type{};
			
			///Trait class that identifies whether the final entry in Args is convertible to ExpectedType
			template<typename FirstArg, typename... Args>
			class last_arg_is_amount : public last_arg_is_amount<Args...>{};
			///Base case of last_arg_convertible (just std::is_convertible)
			template<typename LastArg>
			class last_arg_is_amount<LastArg> : public is_amount<LastArg>{};
			
			template<typename DataType, typename=typename std::enable_if<std::is_arithmetic<DataType>::value> >
			class histogramTraits{
			public:
				using amount=phys_tools::histograms::detail::amount_type<DataType>;
				static void defaultData(double* data, unsigned int count){
					memset(data,0,count*sizeof(DataType));
				}
				static DataType unit(){ return DataType(1); }
				constexpr static bool enable_automatic_amount_handling=true;
			};
			
			//this is completely ridiculous
			//Also, this needs to be outside histogramTraits to avoid every specialization of histogramTraits 
			//having to redefine it. Maybe this could be fixed with a base class?
			//	histogramTraits : public histogramTraitsBase<DataType> ?
			template<typename DataType>
			struct amountAdder{
				amountAdder(DataType& v, const typename detail::histogramTraits<DataType>::amount& a){
					v+=a.value;
				}
				static void add(DataType& v, const typename detail::histogramTraits<DataType>::amount& a){
					v+=a.value;
				}
			};
		} //namespace detail
		
		template<typename DataType, typename=typename std::enable_if<detail::histogramTraits<DataType>::enable_automatic_amount_handling> >
		DataType& operator+=(DataType& v, const typename detail::histogramTraits<DataType>::amount& a){
			detail::amountAdder<DataType> adder(v,a);
			return(v);
		}
		
		///A convenience function for adding non-unit amounts to histograms.
		///This function allows one to write
		///\code
		///histogram<2> h(/*some axes*/);
		///h.add(1,2,amount(5));
		///\endcode
		///instead of the lengthier
		///\code
		///h.add(1,2,histogram<2>::amount(5));
		///\endcode
		///Note, however, that this cannot necessarily be used for complicated bin types which
		///communicate extra information in the amount, like meanVarTracker.
		template<typename DataType, typename=typename std::enable_if<detail::histogramTraits<DataType>::enable_automatic_amount_handling> >
		typename detail::histogramTraits<DataType>::amount amount(const DataType& v){
			return(typename detail::histogramTraits<DataType>::amount(v));
		}
		
		///Used to indicate that the dimensions of a histogram are determined at runtime.
		const int Dynamic = -1;
		
		namespace detail{
			///A traits-like class which enables uniformly extracting the dimensionality of
			///compile-time fixed or runtime dynamic histograms.
			template <template<int,class> class H, int D, class S>
			struct dimensionExtractor{
				static unsigned int getDim(const H<D,S>*){
					return(H<D,S>::dimensions);
				}
			};
			
			template <template<int,class> class H, class S>
			struct dimensionExtractor<H,Dynamic,S>{
				static unsigned int getDim(const H<Dynamic,S>* h){
					return(h->dimensions);
				}
			};
		} //namesapce detail
		
		//fwd decl
		///\cond
		template<int N, typename StoreType>
		class histogram;
		///\endcond
		
		///The base type which defines common operations internal to all histograms. Users should not need to use this class directly.
		template<int StaticDim, class StoreType, typename IC, typename EC, typename A>
		class histogramBase{
		protected:
			bool initialized;
			bool useContentScaling;
			typedef histogram<StaticDim,StoreType> Derived;
			typedef detail::dimensionExtractor<histogram,StaticDim,StoreType> DimExt;
			typedef StoreType dataType;
			typedef IC internalCoordinate;
			typedef EC externalCoordinate;
			typedef A amount;
		
		public:
			///Set the axis in dimension d of this distogram to be a.
			///\param the dimension for which to set the axis
			///\param a an axis object allocated with new, ownership of which is transferred to this histogram
			///\warning Currently it is assumed internally that once this function has been called for the
			///         histogram's final dimension it must have been previously called for all other dimensions
			///         and all axes are safely initialized.
			void setAxis(unsigned int d, axis* a){
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				Derived* derivedThis=static_cast<Derived*>(this);
				assert(d<dimensions);
				derivedThis->axes[d]=a;
				useContentScaling|=a->adviseContentScaling();
				derivedThis->count[d]=derivedThis->axes[d]->numBins();
				if(d==dimensions-1){
					initialized=true;
					size_t t=std::accumulate(derivedThis->count,derivedThis->count+dimensions,1U,std::multiplies<size_t>());
					if(t!=0){
						//std::cout << "all axes initialized, allocating " << t << " bins" << std::endl;
						delete[] derivedThis->data;
						derivedThis->data = new dataType[t];
						detail::histogramTraits<dataType>::defaultData(derivedThis->data,t);
					}
				}
			}
			
		protected:
			///Add k bins in dimension d before the existing bins
			void prependBins(unsigned int d, unsigned int k){
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				Derived* derivedThis=static_cast<Derived*>(this);
				
				//std::cout << "prepending " << k << " bin(s) in dimension " << d << std::endl;
				if(derivedThis->count[d]==0){
					//std::cout << "axis is empty, shifting base" << std::endl;
					detail::histogram_access::shiftBase(*derivedThis->axes[d],-k);
					detail::histogram_access::appendBins(*derivedThis->axes[d],1);
					derivedThis->count[d]=derivedThis->axes[d]->numBins();
					internalCoordinate t=std::accumulate(derivedThis->count,derivedThis->count+dimensions,1U,std::multiplies<internalCoordinate>());
					if(t!=0){
						//std::cout << "all axes initialized, allocating" << std::endl;
						delete[] derivedThis->data;
						derivedThis->data = new dataType[t];
						detail::histogramTraits<dataType>::defaultData(derivedThis->data,t);
					}
					return;
				}
				unsigned int c=1;
				if(d>0)
					c=std::accumulate(derivedThis->count,derivedThis->count+d,1,std::multiplies<unsigned int>());
				unsigned int s=c*derivedThis->count[d];
				unsigned int r=1;
				r=std::accumulate(derivedThis->count+(d+1),derivedThis->count+dimensions,1U,std::multiplies<unsigned int>());
				dataType* new_data=new dataType[(k*c+s)*r];
				detail::histogramTraits<dataType>::defaultData(new_data,(k*c+s)*r);
				dataType* new_ptr=new_data, * old_ptr=derivedThis->data;
				while(r!=0){
					new_ptr+=k*c;
					std::copy(old_ptr,old_ptr+s,new_ptr);
					new_ptr+=s;
					old_ptr+=s;
					r--;
				}
				delete[] derivedThis->data;
				derivedThis->data=new_data;
				detail::histogram_access::prependBins(*derivedThis->axes[d],k);
				derivedThis->count[d]+=k;
			}
			
			///Add k bins in dimension d after the existing bins
			void appendBins(unsigned int d, unsigned int k){
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				Derived* derivedThis=static_cast<Derived*>(this);
				
				//std::cout << "appending " << k << " bin(s) in dimension " << d << std::endl;
				if(derivedThis->count[d]==0){
					//std::cout << "axis is empty, shifting base" << std::endl;
					detail::histogram_access::shiftBase(*derivedThis->axes[d],k-1);
					detail::histogram_access::appendBins(*derivedThis->axes[d],1);
					derivedThis->count[d]=derivedThis->axes[d]->numBins();
					unsigned int t=std::accumulate(derivedThis->count,derivedThis->count+dimensions,1U,std::multiplies<unsigned int>());
					if(t!=0){
						//std::cout << "all axes initialized, allocating " << t << " bins" << std::endl;
						delete[] derivedThis->data;
						derivedThis->data = new dataType[t];
						detail::histogramTraits<dataType>::defaultData(derivedThis->data,t);
					}
					return;
				}
				unsigned int c=1;
				if(d>0)
					c=std::accumulate(derivedThis->count,derivedThis->count+d,1,std::multiplies<unsigned int>());
				unsigned int s=c*derivedThis->count[d];
				unsigned int r=1;
				r=std::accumulate(derivedThis->count+(d+1),derivedThis->count+dimensions,1U,std::multiplies<unsigned int>());
				//std::cout << "allocating " << ((k*c+s)*r) << " bins" << std::endl;
				dataType* new_data=new dataType[(k*c+s)*r];
				detail::histogramTraits<dataType>::defaultData(new_data,(k*c+s)*r);
				dataType* new_ptr=new_data, * old_ptr=derivedThis->data;
				while(r!=0){
					std::copy(old_ptr,old_ptr+s,new_ptr);
					new_ptr+=s;
					new_ptr+=k*c;
					old_ptr+=s;
					r--;
				}
				delete[] derivedThis->data;
				derivedThis->data=new_data;
				detail::histogram_access::appendBins(*derivedThis->axes[d],k);
				derivedThis->count[d]+=k;
			}
			
			///Fetch the contents of the bins at the given internal coordinates
			const StoreType& getBin(const internalCoordinate coordinates[]) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				
				internalCoordinate idx=0;
				internalCoordinate sp=1;
				for(unsigned int i=0; i<dimensions; i++){
					idx+=coordinates[i]*sp;
					sp*=derivedThis->count[i];
				}
				return(derivedThis->data[idx]);
			}
			
			//non-const raw access
			///Fetch the contents of the bins at the given internal coordinates
			StoreType& getBin(const internalCoordinate coordinates[]){
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				Derived* derivedThis=static_cast<Derived*>(this);
				
				internalCoordinate idx=0;
				internalCoordinate sp=1;
				for(unsigned int i=0; i<dimensions; i++){
					idx+=coordinates[i]*sp;
					sp*=derivedThis->count[i];
				}
				return(derivedThis->data[idx]);
			}
			
			///convert a set of external coordinates to internal coordinates
			///\return false if the conversion fails due to underflow or overflow
			//TODO: proplerly handle (hide) prepend and append directives
			bool getCoordinates(externalCoordinate coordinates[], internalCoordinate internal[]) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				
				for(unsigned int i=0; i<dimensions; i++){
					std::pair<axis::lookupResult,internalCoordinate> result=derivedThis->axes[i]->findBin(coordinates[i]);
					if(result.first==axis::UNDERFLOW_ || result.first==axis::OVERFLOW_)
						return(false);
					internal[i]=result.second;
				}
				return(true);
			}
			
			///get the data in a bin using external, rather than internal, coordinates
			const StoreType& getBin(externalCoordinate coordinates[]){
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				internalCoordinate internal[dimensions];
				getCoordinates(coordinates,internal);
				return(getBin(internal));
			}
			
			///core implementation of add
			void addCore(const externalCoordinate coordinate[], amount value){
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				Derived* derivedThis=static_cast<Derived*>(this);
				
				assert(initialized);
				
				//std::cout << " inserting value " /*<< value <<*/ " at (";
				//for(unsigned int i=0; i<dimensions; i++){
				//	if(i)
				//		std::cout << ',';
				//	std::cout << coordinate[i];
				//}
				//std::cout << ')' << std::endl;
				
				internalCoordinate internal[dimensions];
				for(unsigned int i=0; i<dimensions; i++){
					std::pair<axis::lookupResult,unsigned int> r=derivedThis->axes[i]->findBin(coordinate[i]);
					switch(r.first){
						case axis::FOUND_BIN:
							//std::cout << "  dimension " << i << ": FOUND_BIN (" << r.second << ")\n";
							internal[i]=r.second;
							break;
						case axis::PREPEND_BINS:
							//std::cout << "  dimension " << i << ": PREPEND_BINS (" << r.second << ")\n";
							prependBins(i,r.second);
							internal[i]=0;
							break;
						case axis::APPEND_BINS:
							//std::cout << "  dimension " << i << ": APPEND_BINS (" << r.second << ")\n";
							appendBins(i,r.second);
							internal[i]=derivedThis->count[i]-1;
							break;
							//the first underflow or overflow encountered wins
							//this shouldn't matter much since people probably are only interested in
							//distinguishing the two in the 1D case, anyway
						case axis::UNDERFLOW_:
							//std::cout << "  dimension " << i << ": UNDERFLOW\n";
							//derivedThis->underflow+=value;
							detail::amountAdder<StoreType>::add(derivedThis->underflow,value);
							return;
						case axis::OVERFLOW_:
							//std::cout << "  dimension " << i << ": OVERFLOW\n";
							//derivedThis->overflow+=value;
							detail::amountAdder<StoreType>::add(derivedThis->overflow,value);
							return;
					}
				}
				
				//std::cout << " internal coordinates are (";
				//for(unsigned int i=0; i<dimensions; i++){
				//	if(i)
				//		std::cout << ',';
				//	std::cout << internal[i];
				//}
				//std::cout << ')' << std::endl;
				
				unsigned int idx=0;
				unsigned int sp=1;
				for(unsigned int i=0; i<dimensions; i++){
					idx+=internal[i]*sp;
					sp*=derivedThis->count[i];
				}
				//std::cout << " bin is " << idx << std::endl;
				//derivedThis->data[idx]+=value;
				detail::amountAdder<StoreType>::add(derivedThis->data[idx],value);
			}
			
		public:
			histogramBase():initialized(false),useContentScaling(false){}
			histogramBase(bool i, bool u):initialized(i),useContentScaling(u){}
			
			///Get the number of dimensions this histogram has
			unsigned int getDimensions() const{
				return(DimExt::getDim(static_cast<const Derived*>(this)));
			}
			
			///Get whether histogram bin contents are scaled by bin volume
			bool getUseContentScaling() const{
				return(useContentScaling);
			}
			///Set whether histogram bin contents are scaled by bin volume.
			///By default the choice to use scaling is made by taking the disjunction of all axes' guidance from axis::adviseContentScaling(). 
			void setUseContentScaling(bool u){
				useContentScaling=u;
			}
			
			bool isInitialized() const{
				return(initialized);
			}
			
			///Get the range of external coordinates currently within the extent of one of the histogram's axes
			///\param dim the dimension for which to fetch the range
			std::pair<externalCoordinate,externalCoordinate> range(unsigned int dim) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				assert(dim<dimensions);
				return(std::make_pair(derivedThis->axes[dim]->binEdge(0),
									  derivedThis->axes[dim]->binEdge(derivedThis->count[dim])));
			}
			
			///Get the number of bins which currently exist in dimension dim
			unsigned int getBinCount(internalCoordinate dim) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				assert(dim<dimensions);
				return(derivedThis->count[dim]);
			}
			
			///Get the external coordinate of the 'left' edge of the given bin in the given dimension
			///\param dim the dimension
			///\param bin the index of the bin in the given dimension
			externalCoordinate getBinEdge(unsigned int dim, internalCoordinate bin) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				assert(dim<dimensions);
				return(derivedThis->axes[dim]->binEdge(bin));
			}
			
			///Get the width of the given bin in the given dimension in external coordinates
			///\param dim the dimension
			///\param bin the index of the bin in the given dimension
			externalCoordinate getBinWidth(unsigned int dim, internalCoordinate bin) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				assert(dim<dimensions);
				return(derivedThis->axes[dim]->binWidth(bin));
			}
			
			///Get the external coordinate of the center of the given bin in the given dimension
			///\param dim the dimension
			///\param bin the index of the bin in the given dimension
			externalCoordinate getBinCenter(unsigned int dim, internalCoordinate bin) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				assert(dim<dimensions);
				return(derivedThis->axes[dim]->binCenter(bin));
			}
			
			///Get the axis corresponding to the given dimension
			const axis* getAxis(unsigned int dim) const{
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				assert(dim<dimensions);
				return(derivedThis->axes[dim]);
			}
			
			///Get the axis corresponding to the given dimension
			axis* getAxis(unsigned int dim){
				const unsigned int dimensions=DimExt::getDim(static_cast<const Derived*>(this));
				const Derived* derivedThis=static_cast<const Derived*>(this);
				assert(dim<dimensions);
				return(derivedThis->axes[dim]);
			}
			
			const dataType* getData() const{
				const Derived* derivedThis=static_cast<const Derived*>(this);
				return(derivedThis->data);
			}
			
			///A proxy object to handle direct stores to bins when content scaling is in effect
			struct modificationProxy{
			private:
				StoreType* value;
				double scale;
				modificationProxy(StoreType& v, double s):value(&v),scale(s){}
				friend Derived;
			public:
				modificationProxy():value(nullptr),scale(0.0){}
				modificationProxy(const modificationProxy& other):value(other.value),scale(other.scale){}
				operator StoreType() const{ return(*value*scale); }
				const modificationProxy& operator=(StoreType v){ *value=v; return(*this); }
				const modificationProxy& operator+=(StoreType v){ *value+=v; return(*this); }
				const modificationProxy& operator-=(StoreType v){ *value-=v; return(*this); }
				template<typename U>
				const modificationProxy& operator*=(const U& v){ *value*=v; return(*this); }
				template<typename U>
				const modificationProxy& operator/=(const U& v){ *value/=v; return(*this); }
			};
			
			///\brief wite this histogram to an open HDF5 file
			///\param container the location within the file to which to write
			///\param name the name used for the histogram
			///\param compressLevel level of gzip compression to apply
			///\param dashiCompat convert all axes to collections of bin edges to increase compatibility with Dashi
			void write(hid_t container, const std::string& name, unsigned int compressLevel=0, bool dashiCompat=false){
				using namespace hdf_interface;
				const Derived* derivedThis=static_cast<const Derived*>(this);
				
				const size_t dim=getDimensions();
				std::unique_ptr<hsize_t[]> dims(new hsize_t[dim]);
				for(unsigned int i=0; i<dim; i++)
					dims[i]=getBinCount(i);
				
				const HDFDatatype& dType=getHDFDatatype<StoreType>();
				hid_t group_id=H5Gcreate(container, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				
				addAttribute(group_id,"rank",dim);
				
				hid_t plist_id=H5P_DEFAULT;
				if(compressLevel){
					assert(compressLevel<=9);
					plist_id = H5Pcreate(H5P_DATASET_CREATE);
					std::unique_ptr<hsize_t[]> chunkDims(new hsize_t[dim]);
					unsigned long size=sizeof(StoreType);
					const unsigned long targetChunkSize=1UL<<16; //64KB
					for(unsigned int i=dim; i>0; i--){
						if(size<targetChunkSize){
							unsigned int dimi=targetChunkSize/size;
							if(dimi>dims[i-1] || dimi==0)
								dimi=dims[i-1];
							chunkDims[i-1]=dimi;
							size*=chunkDims[i-1];
						}
						else
							chunkDims[i-1]=1;
					}
					H5Pset_chunk(plist_id, dim, chunkDims.get());
					H5Pset_deflate(plist_id, compressLevel);
				}
				hid_t dataspace_id = H5Screate_simple(dim, dims.get(), NULL);
				hid_t dataset_id = H5Dcreate(group_id, "_h_bincontent", dType.datatype, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
				
				H5Dwrite(dataset_id, dType.datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, getData());
				
				addAttribute(dataset_id,"overflow",derivedThis->getOverflow());
				addAttribute(dataset_id,"underflow",derivedThis->getUnderflow());
				
				if(compressLevel)
					H5Pclose(plist_id);
				H5Dclose(dataset_id);
				H5Sclose(dataspace_id);
				
				for(unsigned int i=0; i<dim; i++){
					if(dashiCompat)
						FixedUserDefinedAxis(*getAxis(i)).writeHDF(group_id,i);
					else
						getAxis(i)->writeHDF(group_id,i);
				}
				
				H5Gclose(group_id);
			}
			
			///\brief read this histogram from an open HDF5 file
			///\param container the location within the file from which to read
			///\param name the name used for the histogram
			void read(hid_t container, const std::string& name){
				using namespace hdf_interface;
				Derived* derivedThis=static_cast<Derived*>(this);
				
				if(initialized)
					throw std::logic_error("Reading into an already initialized histogram is not supported");
				hid_t group_id=H5Gopen(container, name.c_str(), H5P_DEFAULT);
				
				if(group_id<0)
					throw std::runtime_error("Failed to open group for histogram '"+name+"'");
				
				try{
					//read and check dimension
					size_t dim;
					readAttribute(group_id,"rank",dim);
					if(dim!=getDimensions())
						derivedThis->mismatchedDimensionsOnRead(name,dim);
					
					std::set<std::string> items=groupContents(group_id);
					
					//read axes
					for(unsigned int i=0; i<dim; i++){
						std::string name="_h_axis_"+std::to_string(i);
						if(items.count(name)){ //normal axis
							setAxis(i,deserializeAxis(group_id, name));
							continue;
						}
						name="_h_binedges_"+std::to_string(i);
						if(items.count(name)){ //Dashi compatible axis
							setAxis(i,deserializeAxis(group_id, name));
							continue;
						}
						throw std::runtime_error("Unable to find axis "+std::to_string(i)+" for histogram '"+name+"'");
					}
					
					//read data
					hid_t dataset_id = H5Dopen(group_id, "_h_bincontent", H5P_DEFAULT);
					hid_t dataspace_id = H5Dget_space(dataset_id);
					
					int dataDim = H5Sget_simple_extent_ndims(dataspace_id);
					if(dataDim!=getDimensions()){
						H5Dclose(dataset_id);
						H5Sclose(dataspace_id);
						throw std::runtime_error("Data dimension mismatch for histogram '"+name+
												 "' (expected "+std::to_string(getDimensions())+
												 ", got "+std::to_string(dataDim)+")");
					}
					
					std::unique_ptr<hsize_t[]> extents(new hsize_t[dataDim]);
					H5Sget_simple_extent_dims(dataspace_id,extents.get(),NULL);
					H5Sclose(dataspace_id);
					for(unsigned int i=0; i<getDimensions(); i++){
						if(extents[i]!=getBinCount(i)){
							H5Dclose(dataset_id);
							throw std::runtime_error("Size mismatch for histogram '"+name+
													 "' in dimension "+std::to_string(i)+
													 " (expected "+std::to_string(getBinCount(i))+
													 ", got "+std::to_string(extents[i])+")");
						}
					}
					
					const HDFDatatype& dType=getHDFDatatype<StoreType>();
					herr_t status = H5Dread(dataset_id, dType.datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, derivedThis->data);
					if(status<0){
						H5Dclose(dataset_id);
						throw std::runtime_error("Failed to contents of histogram '"+name+"'");
					}
					H5Dclose(dataset_id);
				}catch(...){
					H5Gclose(group_id);
					throw;
				}
				
				H5Gclose(group_id);
			}
		};
		
		namespace detail{
			//fwd decl
			///\cond
			template<int N, typename StoreType>
			struct histogramPrinter;
			///\endcond
		}
		
		//StoreType must be default constructable, assignable, have a valid StoreType+=double operation
		//It must also be possible to multiply and divide a StoreType by a double for content scaling
		///A multidimensional histogram
		///\param N the number of dimensions covered by the histogram.
		///			The special value Dynamic is used to indicate that the dimensionality will not be determined until runtime. 
		///\param StoreType the type to use for each bin in the histogram
		template<int N, typename StoreType=double>
		class histogram : public histogramBase<N,StoreType,unsigned int,double,typename detail::histogramTraits<StoreType>::amount>{
		public:
			static_assert(N>0,"Histograms must have at least one dimension");
			///The dimensionalitly of the histogram
			enum:unsigned int{ dimensions=N };
			///The type used to represent the histogram bins
			typedef StoreType dataType;
			///The type used to represent bin indices
			typedef unsigned int internalCoordinate;
			///The type used to represent external coordinates
			typedef double externalCoordinate;
			///The type used to specify arbitrary amounts to place in bins
			typedef typename detail::histogramTraits<StoreType>::amount amount;
			typedef histogramBase<N,StoreType,internalCoordinate,externalCoordinate,amount> BaseType;
			
			using BaseType::setAxis;
			using BaseType::prependBins;
			using BaseType::appendBins;
			using BaseType::getBin;
			using BaseType::getCoordinates;
			using BaseType::addCore;
			using BaseType::isInitialized;
			using BaseType::range;
			using BaseType::getBinCount;
			using BaseType::getBinEdge;
			using BaseType::getBinWidth;
			using BaseType::getBinCenter;
			using BaseType::getAxis;
			using modificationProxy=typename BaseType::modificationProxy;
			
			//the base needs access to private members
			friend BaseType;
			//all histograms are friends
			template<int M, typename OtherStoreType>
			friend class histogram;
			
		private:
			axis* axes[N];
			internalCoordinate count[N];
			StoreType* data;
			StoreType underflow, overflow;
		public:
			///Construct a histogram without axes
			histogram():
			data(new StoreType[1]),
			underflow(0),overflow(0)
			{
				std::fill_n(&axes[0],N,nullptr);
				std::fill_n(&count[0],N,0);
			}
			
			///Construct a histogram from a set of axes
			template<typename Axis, typename... Axes, typename = typename
        std::enable_if<!(
            sizeof...(Axes) == 0 &&
            std::is_same<typename
                std::remove_reference<Axis>::type,
                std::tuple<histogram<N, StoreType>>
            >::value)
        >::type>
			histogram(Axis axis, Axes... axes):
			data(new StoreType[1]),
			underflow(0),overflow(0){
				static_assert(sizeof...(axes)+1==dimensions,
				              "Number of axes used to construct a histogram "
				              "must be equal to its number of dimensions");
				//TODO: this is wrong because it disallows passing axis*s obtained from axis::copy()
				//static_assert(all_args_convertible<axis&,typename std::add_lvalue_reference<Axes>::type...>::value,
				//              "all arguments to histogram(Axes... axes) "
				//			  "must be axis objects");
				setAxes_impl<0>(axis, axes...);
			}
			
			///Copy construct a histogram from another of the same type.
			///All axes and data are copied. 
			histogram<N,StoreType>(const histogram<N,StoreType>& other):
			BaseType(other.initialized,other.useContentScaling),
			data(new StoreType[std::max(1U,std::accumulate(other.count,other.count+N,1U,std::multiplies<unsigned int>()))]),
			underflow(other.underflow),
			overflow(other.overflow){
				if(!other.initialized){
					//other may have some initialized axes, so we need to copy those
					for(unsigned int i=0; i<N; i++){
						if(other.axes[i])
							setAxis(i,other.axes[i]->copy());
						else
							axes[i]=nullptr;
					}
					
					return;
				}
				std::copy(&other.count[0],&other.count[N],&count[0]);
				unsigned long totalSize=1;
				for(unsigned int i=0; i<N; i++){
					axes[i]=other.axes[i]->copy();
					totalSize*=count[i];
					if(count[i]>0)
						detail::histogram_access::appendBins(*axes[i],(count[i]));
				}
				std::copy(&other.data[0],&other.data[totalSize],&data[0]);
			}
			
			///Move construct a histogram from another of the same type.
			///All axes and data are taken from the other histogram and moved into the newly
			///constructed object, leaving the original empty and unusable.
			histogram<N,StoreType>(histogram<N,StoreType>&& other):
			BaseType(other.initialized,other.useContentScaling),
			data(other.data),
			underflow(other.underflow),
			overflow(other.overflow){
				other.data=new StoreType[1];
				other.initialized=false;
				other.underflow=StoreType();
				other.overflow=StoreType();
				std::copy(&other.count[0],&other.count[N],&count[0]);
				std::fill_n(&other.count[0],N,0);
				for(unsigned int i=0; i<N; i++){
					axes[i]=other.axes[i];
					other.axes[i]=nullptr;
				}
			}
			
			///Copy from another histogram of the same type.
			///All axes and data are copied.
			histogram<N,StoreType>& operator=(const histogram<N,StoreType>& other){
				if(&other==this)
					return(*this);
				delete[] data;
				for(auto a : axes)
					delete a;
				
				this->useContentScaling=other.useContentScaling;
				this->initialized=other.initialized;
				underflow=other.underflow;
				overflow=other.overflow;
				std::copy(&other.count[0],&other.count[N],&count[0]);
				for(unsigned int i=0; i<N; i++){
					axes[i]=other.axes[i]->copy();
					if(count[i]>0)
						detail::histogram_access::appendBins(*axes[i],count[i]);
						//axes[i]->appendBins(count[i]);
				}
				unsigned long totalSize=std::accumulate(other.count,other.count+N,1U,std::multiplies<unsigned int>());
				data=new StoreType[totalSize];
				std::copy(&other.data[0],&other.data[totalSize],&data[0]);
				
				return(*this);
			}
			
			///Move from another histogram of the same type.
			///All axes and data are taken from the other histogram and moved into this
			///histogram, leaving the other empty an unusable.
			histogram<N,StoreType>& operator=(histogram<N,StoreType>&& other){
				if(&other==this)
					return(*this);
				delete[] data;
				for(auto a : axes)
					delete a;
				data=other.data;
				other.data=new StoreType[1];
				this->useContentScaling=other.useContentScaling;
				other.useContentScaling=false;
				this->initialized=other.initialized;
				other.initialized=false;
				underflow=other.underflow;
				other.underflow=StoreType(0);
				overflow=other.overflow;
				other.overflow=StoreType(0);
				std::copy(&other.count[0],&other.count[N],&count[0]);
				std::fill_n(&other.count[0],N,0);
				for(unsigned int i=0; i<N; i++){
					axes[i]=other.axes[i];
					other.axes[i]=nullptr;
				}
				return(*this);
			}
			
			~histogram(){
				delete[] data;
				for(auto a : axes)
					delete a;
			}
			
		private:
			//internal functions for setting axes
			template<unsigned int index>
			void setAxes_impl() const{
				static_assert(index==N,"number of axes passed to histogram<>::setAxes() must match the dimension of the histogram");
			}
			
			template<unsigned int index,typename AxisType,typename... Args>
			void setAxes_impl(AxisType a, Args... args){
				setAxis(index,new AxisType(a));
				setAxes_impl<index+1>(args...);
			}
			
			template<unsigned int index,typename... Args>
			void setAxes_impl(axis* a, Args... args){
				setAxis(index,a);
				setAxes_impl<index+1>(args...);
			}
			
			template<unsigned int index,typename AxisType,typename... Args>
			void setAxes_impl(const AxisType* a, Args... args){
				setAxis(index,a->copy());
				setAxes_impl<index+1>(args...);
			}
		public:	
			///Set all of the axes of the histogram
			template<typename... Axes>
			void setAxes(Axes... axes){
				static_assert(sizeof...(axes)==N,"number of axes passed to histogram<>::setAxes() must match the dimension of the histogram");
				static_assert(all_args_convertible<axis&,typename std::add_lvalue_reference<Axes>::type...>::value,
				              "all arguments to histogram(Axes... axes) "
							  "must be axis objects");
				setAxes_impl<0>(axes...);
			}
			
			///Assign the given value to every bin in the histogram
			void fill(const StoreType& v){
				unsigned int n=std::accumulate(count,count+dimensions,1U,std::multiplies<unsigned int>());
				std::fill_n(&data[0],n,v);
			}
			
			//enable printing!
			template<int, typename>
			friend struct detail::histogramPrinter;
			
			template<int N_, typename StoreType_>
			friend histogram<N_,StoreType_> operator/(const histogram<N_,StoreType_>& h1, const histogram<N_,StoreType_>& h2);
			
		private:
			//internal machinery for adding and accessing data
			
			//base case for default amount (1)
			template<unsigned int index>
			void addImpl(externalCoordinate (&coordinates)[dimensions]){
				addCore(coordinates,detail::histogramTraits<StoreType>::unit());
			}
			//base case for specified amount
			template<unsigned int index>
			void addImpl(externalCoordinate (&coordinates)[dimensions], amount a){
				addCore(coordinates,a);
			}
			template<unsigned int index,typename... Args>
			void addImpl(externalCoordinate (&coordinates)[dimensions], externalCoordinate nextCoord, Args... args){
				assert(index<dimensions);
				coordinates[index]=nextCoord;
				addImpl<index+1>(coordinates,args...);
			}
			
			template<unsigned int index>
			StoreType subscriptImpl(internalCoordinate (&coordinates)[dimensions], double scale) const{
				StoreType temp=getBin(coordinates);
				return(temp*scale);
			}
			template<unsigned int index,typename... Args>
			StoreType subscriptImpl(internalCoordinate (&coordinates)[dimensions], double scale, internalCoordinate nextCoord, Args... args) const{
				coordinates[index]=nextCoord;
				if(this->useContentScaling)
					scale/=axes[index]->binWidth(nextCoord);
				StoreType temp=subscriptImpl<index+1>(coordinates,scale,args...);
				return(temp);
			}
			
			//non-const for raw access
			template<unsigned int index>
			modificationProxy subscriptImplRaw(const internalCoordinate (&coordinates)[dimensions], double scale){
				return(modificationProxy(getBin(coordinates),scale));
			}
			template<unsigned int index,typename... Args>
			modificationProxy subscriptImplRaw(internalCoordinate (&coordinates)[dimensions], double scale, internalCoordinate nextCoord, Args... args){
				coordinates[index]=nextCoord;
				if(this->useContentScaling)
					scale/=axes[index]->binWidth(nextCoord);
				return(subscriptImplRaw<index+1>(coordinates,scale,args...));
			}
			
			///Throw an exception when rank of the histogram does not match a serialized state being read from a file
			void mismatchedDimensionsOnRead(const std::string& name, size_t dim){
				throw std::runtime_error("Dimension mismatch for histogram '"+name+
				                         "' (expected "+std::to_string(this->getDimensions())+
				                         ", got "+std::to_string(dim)+")");
			}
			
		public:
			//public interface for manipulating the data
			
			///\brief Insert a count into the histogram.
			///
			///The first N arguments must be the (external) coordinates at which the count is to be inserted,
			///and the optional last argument is the amount to insert, if not a unit count.
			///
			///Usage example:
			///\code
			/// histogram<2> h(LinearAxis(0,1),LinearAxis(0,1));
			/// h.add(1,3); //add a unit count at (1,3)
			/// h.add(2,4,histogram<2>::amount(5)); //add 5 units at (2,4)
			/// cout << h << endl;
			///\endcode
			///which prints:
			///\verbatim
			///1 3 1
			///2 3 0
			///
			///1 4 0
			///2 4 5
			///\endverbatim
			template<typename... Args, typename = typename std::enable_if<not first_is_pointer<Args...>::value>::type>
			void add(Args... args){
				//still need asserts later to verify that in the N case there were N external coordinates
				//or that in the N+1 case there were N external coordinates followed by one amount, but this 
				//will help catch most errors early
				/*static_assert((sizeof...(args)==dimensions && all_args_are<externalCoordinate,Args...>::value)
							  || (sizeof...(args)==dimensions+1 && all_args_except_last_are<externalCoordinate,Args...>::value && last_arg_is<amount,Args...>::value),
							  "incorrect number of arguments to histogram<>::add()");*/
				static_assert(sizeof...(args)==dimensions || sizeof...(args)==dimensions+1, "incorrect number of arguments to histogram<>::add()");
				static_assert((sizeof...(args)==dimensions && all_args_convertible<externalCoordinate,Args...>::value) || sizeof...(args)==dimensions+1,
				              "incorrect type of arguments to histogram<>::add(); when using N arguments they must all be external coordinates");
				static_assert((sizeof...(args)==dimensions+1 && all_args_except_last_convertible<externalCoordinate,Args...>::value) || sizeof...(args)==dimensions,
							  "incorrect type of arguments to histogram<>::add(); when using N+1 arguments all but the last must be external coordinates");
				//static_assert((sizeof...(args)==dimensions+1 && last_arg_convertible<amount,Args...>::value) || sizeof...(args)==dimensions,
				static_assert((sizeof...(args)==dimensions+1 && detail::last_arg_is_amount<Args...>::value) || sizeof...(args)==dimensions,
							  "incorrect type of arguments to histogram<>::add(); when using N+1 arguments the last must be an amount");
				externalCoordinate coordinates[dimensions];
				addImpl<0>(coordinates,args...);
			}
			
			///Add a count to the histogram at a set of external coordinates passed in an array.
			void add(const externalCoordinate coordinates[dimensions], amount a=1.0){
				addCore(coordinates,a);
			}
			
			///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates
			template<typename... Args>
			StoreType operator()(Args... args) const{
				static_assert(sizeof...(args)==dimensions,"incorrect number of arguments to histogram<>::operator()");
				static_assert(all_args_convertible<internalCoordinate,Args...>::value,"all arguments to histogram<>::operator() must be internal coordinates");
				unsigned int coordinates[dimensions];
				StoreType temp=subscriptImpl<0>(coordinates,1.0,args...);
				return(temp);
			}
			///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates, passed in an array
			StoreType operator()(const internalCoordinate (&coordinates)[dimensions]) const{
				double scale=1.0;
				if(this->useContentScaling)
					for(unsigned int i=0; i<dimensions; i++)
						scale/=axes[i]->binWidth(coordinates[i]);
				return(getBin(coordinates)*scale);
			}
			
			///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates
			template<typename... Args>
			modificationProxy operator()(Args... args){
				static_assert(sizeof...(args)==N,"incorrect number of arguments to histogram<>::operator()");
				static_assert(all_args_convertible<internalCoordinate,Args...>::value,"all arguments to histogram<>::operator() must be internal coordinates");
				unsigned int coordinates[N];
				return(subscriptImplRaw<0>(coordinates,1.0,args...));
			}
			///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates, passed in an array
			modificationProxy operator()(const internalCoordinate (&coordinates)[dimensions]){
				double scale=1.0;
				if(this->useContentScaling)
					for(unsigned int i=0; i<N; i++)
						scale/=axes[i]->binWidth(coordinates[i]);
				return(modificationProxy(getBin(coordinates),scale));
			}
			
			///Get the sum of the values in all bins of the histograms,
			///optionally including the underflow and overflow bins.
			StoreType integral(bool includeOverflow=false) const{
				StoreType total=std::accumulate(cbegin(),cend(),StoreType());
				if(includeOverflow)
					total+=underflow+overflow;
				return(total);
			}
			///Get the contents of the the underflow bin
			const StoreType& getUnderflow() const{
				return(underflow);
			}
			///Set the contents of the the underflow bin
			void setUnderflow(const StoreType& value){
				underflow=value;
			}
			///Get the contents of the the overflow bin
			const StoreType& getOverflow() const{
				return(overflow);
			}
			///Set the contents of the the overflow bin
			void setOverflow(const StoreType& value){
				overflow=value;
			}
			
			///Type for forward iterators
			///\param HistType the type of the histogram to which this iterator points
			///\param T the type pointed to by the iterator (the type stored in the histogram)
			///\param DerefType the type yielded by dereferencing the iterator
			///\param PointerType the type used as a pointer to the object to which the iterator points
			template<typename HistType, typename T, typename DerefType=T, typename PointerType=DerefType*>
			struct iteratorTempl : public std::iterator<std::bidirectional_iterator_tag, DerefType, std::ptrdiff_t, PointerType>{
			private:
				HistType* h;
				internalCoordinate coordinates[N];
				T value;
				bool isEnd;
				
				iteratorTempl(HistType& h, bool e):h(&h),isEnd(e){
					std::fill_n(&coordinates[0],N,0);
				}
				
				void fetchValue(){
					if(!isEnd)
						value=(*h)(coordinates);
				}
				
				friend HistType;
				static iteratorTempl make_begin(HistType& h){
					iteratorTempl it(h,std::accumulate(h.count,h.count+dimensions,1U,std::multiplies<internalCoordinate>())==0);
					it.fetchValue();
					return(it);
				}
				static iteratorTempl make_end(HistType& h){
					iteratorTempl it(h,true);
					return(it);
				}
			public:
				typedef HistType HistogramType;
				
				iteratorTempl(const iteratorTempl& it):
				h(it.h),value(it.value),isEnd(it.isEnd){
					if(!isEnd)
						std::copy(&it.coordinates[0],&it.coordinates[0]+N,&coordinates[0]);
				}
				
				iteratorTempl& operator=(const iteratorTempl& it){
					h=it.h;
					isEnd=it.isEnd;
					std::copy(&it.coordinates[0],&it.coordinates[0]+N,&coordinates[0]);
					if(!isEnd)
						value=it.value;
					return(*this);
				}
				
				///prefix increment
				iteratorTempl& operator++(){
					int i;
					for(i=0; i<N; i++){
						coordinates[i]++;
						if(coordinates[i]>=h->count[i])
							coordinates[i]=0;
						else
							break;
					}
					if(i==N)
						isEnd=true;
					else
						fetchValue();
					return(*this);
				}
				
				///postfix increment
				iteratorTempl operator++(int){
					iteratorTempl storeIt(*this);
					++(*this);
					return(storeIt);
				}
				
				///prefix decrement
				iteratorTempl& operator--(){
					if(isEnd)
						isEnd=false;
					int i;
					for(i=0; i<N; i++){
						if(coordinates[i]==0)
							coordinates[i]=h->count[i]-1;
						else{
							coordinates[i]--;
							break;
						}
					}
					fetchValue();
					return(*this);
				}
				
				///postfix decrement
				iteratorTempl operator--(int){
					iteratorTempl storeIt(*this);
					--(*this);
					return(storeIt);
				}
				
				bool operator==(const iteratorTempl& rhs) const{
					if(isEnd!=rhs.isEnd) //an end iterator cannot equal a non-end iterator
						return(false);
					if(isEnd) //all end iterators are equal
						return(true);
					if(h!=rhs.h) //otherwise, iterators to different histograms cannot be equal
						return(false);
					for(int i=0; i<N; i++){
						if(coordinates[i]!=rhs.coordinates[i])
							return(false);
					}
					return(true);
				}
				
				bool operator!=(const iteratorTempl& rhs) const{
					if(isEnd!=rhs.isEnd) //any end iterator is unequal to any non-end iterator
						return(true);
					if(isEnd) //no two end iterators are unequal
						return(false);
					if(h!=rhs.h) //otherwise, iterators to different histograms cannot be equal
						return(true);
					for(int i=0; i<N; i++){
						if(coordinates[i]!=rhs.coordinates[i])
							return(true);
					}
					return(false);
				}
				
				DerefType operator*(){
					return(value);
				}
				PointerType operator->(){
					return(&value);
				}
				
				///Get the internal coordinates within the histogram (bin indices) to which this iterator refers
				const internalCoordinate (&getCoordinates() const)[N]{
					return(coordinates);
				}
				///Get the internal coordinate in dimension dim within the histogram (bin index) to which this iterator refers
				internalCoordinate getCoordinate(unsigned int dim) const{
					assert(dim<N);
					return(coordinates[dim]);
				}
				
				///Get the smaller edge in dimension dim of the bin to which this iterator refers
				externalCoordinate getBinEdge(unsigned int dim) const{
					return(h->getBinEdge(dim,coordinates[dim]));
				}
				///Get the width in dimension dim of the bin to which this iterator refers
				externalCoordinate getBinWidth(unsigned int dim) const{
					return(h->getBinWidth(dim,coordinates[dim]));
				}
				///Get the center in dimension dim of the bin to which this iterator refers
				externalCoordinate getBinCenter(unsigned int dim) const{
					return(h->getBinCenter(dim,coordinates[dim]));
				}
				
				HistType& histogram(){ return(*h); }
				const HistType& histogram() const{ return(*h); }
			};
				
			///Type for reverse iterators
			///\param HistType the type of the histogram to which this iterator points
			///\param T the type pointed to by the iterator (the type stored in the histogram)
			///\param DerefType the type yielded by dereferencing the iterator
			///\param PointerType the type used as a pointer to the object to which the iterator points
			template<typename HistType, typename T, typename DerefType=T, typename PointerType=DerefType*>
			struct rIteratorTempl : public std::iterator<std::bidirectional_iterator_tag, DerefType, std::ptrdiff_t, PointerType>{
			private:
				HistType* h;
				internalCoordinate coordinates[N];
				T value;
				bool isEnd;
				
				rIteratorTempl(HistType& h, bool e):h(&h),isEnd(e){
					std::transform(&h.count[0],&h.count[0]+N,&coordinates[0],
								   [](internalCoordinate c){ return(c-1); });
				}
				
				void fetchValue(){
					if(!isEnd)
						value=(*h)(coordinates);
				}
				
				friend HistType;
				static rIteratorTempl make_begin(HistType& h){
					rIteratorTempl it(h,std::accumulate(h.count,h.count+dimensions,1U,std::multiplies<internalCoordinate>())==0);
					it.fetchValue();
					return(it);
				}
				static rIteratorTempl make_end(HistType& h){
					rIteratorTempl it(h,true);
					return(it);
				}
			public:
				typedef HistType HistogramType;
				
				rIteratorTempl(const rIteratorTempl& it):
				h(it.h),value(it.value),isEnd(it.isEnd){
					std::copy(&it.coordinates[0],&it.coordinates[0]+N,&coordinates[0]);
				}
				
				rIteratorTempl& operator=(const rIteratorTempl& it){
					h=it.h;
					isEnd=it.isEnd;
					std::copy(&it.coordinates[0],&it.coordinates[0]+N,&coordinates[0]);
					if(!isEnd)
						value=it.value;
					return(*this);
				}
				
				///prefix increment
				rIteratorTempl& operator++(){
					int i;
					for(i=0; i<N; i++){
						if(coordinates[i]==0)
							coordinates[i]=h->count[i]-1;
						else{
							coordinates[i]--;
							break;
						}
					}
					if(i==N)
						isEnd=true;
					else
						fetchValue();
					return(*this);
				}
				
				///postfix increment
				rIteratorTempl operator++(int){
					rIteratorTempl storeIt(*this);
					++(*this);
					return(storeIt);
				}
				
				///prefix decrement
				rIteratorTempl& operator--(){
					if(isEnd)
						isEnd=false;
					int i;
					for(i=0; i<N; i++){
						coordinates[i]++;
						if(coordinates[i]>=h->count[i])
							coordinates[i]=0;
						else
							break;
					}
					fetchValue();
					return(*this);
				}
				
				///postfix decrement
				rIteratorTempl operator--(int){
					rIteratorTempl storeIt(*this);
					--(*this);
					return(storeIt);
				}
				
				bool operator==(const rIteratorTempl& rhs) const{
					if(isEnd!=rhs.isEnd) //an end iterator cannot equal a non-end iterator
						return(false);
					if(isEnd) //all end iterators are equal
						return(true);
					for(int i=0; i<N; i++){
						if(coordinates[i]!=rhs.coordinates[i])
							return(false);
					}
					return(true);
				}
				
				bool operator!=(const rIteratorTempl& rhs) const{
					if(isEnd!=rhs.isEnd) //any end iterator is unequal to any non-end iterator
						return(true);
					if(isEnd) //no two end iterators are unequal
						return(false);
					for(int i=0; i<N; i++){
						if(coordinates[i]!=rhs.coordinates[i])
							return(true);
					}
					return(false);
				}
				
				DerefType operator*(){
					return(value);
				}
				PointerType operator->(){
					return(&value);
				}
				
				///Get the internal coordinates within the histogram (bin indices) to which this iterator refers
				const internalCoordinate (&getCoordinates() const)[N]{
					return(coordinates);
				}
				///Get the internal coordinate in dimension dim within the histogram (bin index) to which this iterator refers
				internalCoordinate getCoordinate(unsigned int dim) const{
					assert(dim<N);
					return(coordinates[dim]);
				}
				
				///Get the smaller edge in dimension dim of the bin to which this iterator refers
				externalCoordinate getBinEdge(unsigned int dim) const{
					return(h->getBinEdge(dim,coordinates[dim]));
				}
				///Get the width in dimension dim of the bin to which this iterator refers
				externalCoordinate getBinWidth(unsigned int dim) const{
					return(h->getBinWidth(dim,coordinates[dim]));
				}
				///Get the center in dimension dim of the bin to which this iterator refers
				externalCoordinate getBinCenter(unsigned int dim) const{
					return(h->getBinCenter(dim,coordinates[dim]));
				}
				
				HistType& histogram(){ return(h); }
				const HistType& histogram() const{ return(*h); }
			};
			
			///A bidirectional iterator
			typedef iteratorTempl<histogram,modificationProxy,modificationProxy&,modificationProxy*> iterator;
			///A const bidirectional iterator
			typedef iteratorTempl<const histogram,StoreType,const StoreType&,const StoreType*> const_iterator;
			friend iterator;
			friend const_iterator;
			
			///A bidirectional iterator which runs in reverse
			typedef rIteratorTempl<histogram,modificationProxy,modificationProxy,modificationProxy*> reverse_iterator;
			///A const bidirectional iterator which runs in reverse
			typedef rIteratorTempl<const histogram,StoreType,StoreType&,const StoreType*> const_reverse_iterator;
			friend reverse_iterator;
			friend const_reverse_iterator;
			
			///Get an iterator to the first bin in the histogram
			iterator begin(){
				return(iterator::make_begin(*this));
			}
			///Get an iterator past the last bin in the histogram
			iterator end(){
				return(iterator::make_end(*this));
			}
			///Get an iterator to the first bin in the histogram
			const_iterator begin() const{
				return(const_iterator::make_begin(*this));
			}
			///Get an iterator past the last bin in the histogram
			const_iterator end() const{
				return(const_iterator::make_end(*this));
			}
			///Get a const iterator to the first bin in the histogram
			const_iterator cbegin() const{
				return(const_iterator::make_begin(*this));
			}
			///Get a const iterator past the last bin in the histogram
			const_iterator cend() const{
				return(const_iterator::make_end(*this));
			}
			///Get a reverse iterator to the last bin of the histogram
			reverse_iterator rbegin(){
				return(reverse_iterator::make_begin(*this));
			}
			///Get a reverse iterator before the first bin of the histogram
			reverse_iterator rend(){
				return(reverse_iterator::make_end(*this));
			}
			///Get a reverse iterator to the last bin of the histogram
			const_reverse_iterator rbegin() const{
				return(const_reverse_iterator::make_begin(*this));
			}
			///Get a reverse iterator before the first bin of the histogram
			const_reverse_iterator rend() const{
				return(const_reverse_iterator::make_end(*this));
			}
			///Get a const reverse iterator to the last bin of the histogram
			const_reverse_iterator crbegin() const{
				return(const_reverse_iterator::make_begin(*this));
			}
			///Get a const reverse iterator before the first bin of the histogram
			const_reverse_iterator crend() const{
				return(const_reverse_iterator::make_end(*this));
			}
			
		private:
			//base case
			template<unsigned int index>
			void findBinImpl(iterator& it){
				it.fetchValue();
			}
				
			//recursive case
			template<unsigned int index, typename... Args>
			void findBinImpl(iterator& it, externalCoordinate nextCoord, Args... args){
				auto result=axes[index]->findBin(nextCoord);
				if(result.first!=axis::FOUND_BIN){
					it.isEnd=true;
					return;
				}
				it.coordinates[index]=result.second;
				findBinImpl<index+1>(it,args...);
			}
			
			//base case
			template<unsigned int index>
			void findBinImpl(const_iterator& it) const{
				it.fetchValue();
			}
			
			//recursive case
			template<unsigned int index, typename... Args>
			void findBinImpl(const_iterator& it, externalCoordinate nextCoord, Args... args) const{
				auto result=axes[index]->findBin(nextCoord);
				if(result.first!=axis::FOUND_BIN){
					it.isEnd=true;
					return;
				}
				it.coordinates[index]=result.second;
				findBinImpl<index+1>(it,args...);
			}
			
		public:
			///Get an iterator to the bin in which an entry at the given coordinates would fall,
			///or an end iterator if such an entry is outside the current extent of the histogram
			template<typename... Args, typename = typename std::enable_if<not first_is_pointer<Args...>::value>::type>
			iterator findBin(Args... coordinates){
				static_assert(sizeof...(coordinates)==N,"incorrect number of arguments to histogram<>::findBin()");
				iterator it(*this,false); //assume that it will not be the end iterator
				findBinImpl<0>(it,coordinates...);
				return(it);
			}
			
			///Get an iterator to the bin in which an entry at the given coordinates would fall,
			///or an end iterator if such an entry is outside the current extent of the histogram
			template<typename... Args, typename = typename std::enable_if<not first_is_pointer<Args...>::value>::type>
			const_iterator findBin(Args... coordinates) const{
				static_assert(sizeof...(coordinates)==N,"incorrect number of arguments to histogram<>::findBin()");
				const_iterator it(*this,false); //assume that it will not be the end iterator
				findBinImpl<0>(it,coordinates...);
				return(it);
			}
			
			///Get an iterator to the bin in which an entry at the coordinates given in an array would fall
			///or an end iterator if such an entry is outside the current extent of the histogram
			iterator findBin(externalCoordinate coordinates[]){
				iterator it(*this,false); //assume that it will not be the end iterator
				for(unsigned int i=0; i<dimensions; i++){
					auto result=axes[i]->findBin(coordinates[i]);
					if(result.first!=axis::FOUND_BIN){
						it.isEnd=true;
						break;
					}
					it.coordinates[i]=result.second;
				}
				it.fetchValue();
				return(it);
			}
			///Get an iterator to the bin in which an entry at the coordinates given in an array would fall
			///or an end iterator if such an entry is outside the current extent of the histogram
			const_iterator findBin(externalCoordinate coordinates[]) const{
				iterator it(*this,false); //assume that it will not be the end iterator
				for(unsigned int i=0; i<dimensions; i++){
					auto result=axes[i]->findBin(coordinates[i]);
					if(result.first!=axis::FOUND_BIN){
						it.isEnd=true;
						break;
					}
					it.coordinates[i]=result.second;
				}
				it.fetchValue();
				return(it);
			}
			
			///Find the bin in which the entry at the center of the bin pointed to by it (which may be an
			///iterator referring to another histogram) would fall in this histogram.
			template<typename IteratorType>
			iterator findBinIterator(IteratorType it){
				static_assert(dimensions==IteratorType::HistogramType::dimensions,"Histogram dimensions must match");
				externalCoordinate loc[dimensions];
				iterator ret(*this,false); //assume that it will not be the end iterator
				for(unsigned int i=0; i<dimensions; i++){
					loc[i]=it.histogram().getBinCenter(i,it.getCoordinate(i));
					auto result=getAxis(i)->findBin(loc[i]);
					if(result.first!=axis::FOUND_BIN){
						ret.isEnd=true;
						break;
					}
					ret.coordinates[i]=result.second;
				}
				ret.fetchValue();
				return(ret);
			}
			
			///Find the bin in which the entry at the center of the bin pointed to by it (which may be an
			///iterator referring to another histogram) would fall in this histogram.
			template<typename IteratorType>
			const_iterator findBinIterator(IteratorType it) const{
				static_assert(dimensions==IteratorType::HistogramType::dimensions,"Histogram dimensions must match");
				double loc;
				const_iterator ret(*this,false); //assume that it will not be the end iterator
				for(unsigned int i=0; i<dimensions; i++){
					loc=it.histogram().getBinCenter(i,it.getCoordinate(i));
					auto result=getAxis(i)->findBin(loc);
					if(result.first!=axis::FOUND_BIN){
						ret.isEnd=true;
						break;
					}
					ret.coordinates[i]=result.second;
				}
				ret.fetchValue();
				return(ret);
			}
			
			///Multiply all bin contents by a
			template<typename U>
			histogram& operator*=(U a){
				unsigned long totalSize=std::accumulate(count,count+N,1UL,std::multiplies<unsigned long>());
				for(unsigned long i=0; i<totalSize; i++)
					data[i]*=a;
				underflow*=a;
				overflow*=a;
				return(*this);
			}
			///Divide all bin contents by a
			template<typename U>
			histogram& operator/=(U a){
				unsigned long totalSize=std::accumulate(count,count+N,1UL,std::multiplies<unsigned long>());
				for(unsigned long i=0; i<totalSize; i++)
					data[i]/=a;
				underflow/=a;
				overflow/=a;
				return(*this);
			}
			
			///Construct a new histogram with one fewer dimension with the projection of this histogram's data
			///\param dim the dimension to be projected out
			histogram<N-1,StoreType> project(unsigned int dim) const{
				static_assert(N>1,"Projecting out a dimension is only defined for histograms with more than one dimension");
				assert(dim<N);
				histogram<N-1,StoreType> h;
				//copy all of the axes except the dim'th into h
				for(unsigned int i=0,j=0;i<N; i++){
					if(i==dim)
						continue;
					h.setAxis(j,getAxis(i)->copy());
					j++; //note that this is skipped on the i==dim iteration
				}
				h.setUseContentScaling(false);
				
				//stretch h to have the correct range
				for(unsigned int i=0,j=0;i<N; i++){
					if(i==dim)
						continue;
					//strange little sidestep here to avoid appendBins messing up the axis range
					h.appendBins(j,1);
					h.appendBins(j,getBinCount(i)-1);
					//std::cout << "Dimension " << i << " (" << j << ")\n";
					//std::cout << " Old range: [" << range(i).first << ',' << range(i).second << "]\n";
					//std::cout << " New range: [" << h.range(j).first << ',' << h.range(j).second << "]\n";
					j++; //note that this is skipped on the i==dim iteration
				}
				
				//actually copy over the data, projecting out dimension dim
				internalCoordinate coordinates[N]={0/*All zeros*/};
				internalCoordinate projectedCoordinates[N-1]={0/*All zeros*/};
				const internalCoordinate n=std::accumulate(count,count+N,1U,std::multiplies<internalCoordinate>());
				for(internalCoordinate i=0; i<n; i++){
					//std::cout << "projected coordinates: (";
					//for(unsigned int k=0; k<N-1; k++)
					//	std::cout << (k?",":"") << projectedCoordinates[k];
					//std::cout << ")\n";
					h.getBin(projectedCoordinates)+=data[i];
					//std::cout << " adding " << data[i] << ", sum now " << h.getBin(projectedCoordinates) << std::endl;
					
					for(unsigned int j=0; j<N; j++){
						coordinates[j]++;
						if(j<dim)
							projectedCoordinates[j]++;
						//do nothing to projectedCoordinates when j==dim
						else if(j>dim)
							projectedCoordinates[j-1]++;
						
						if(coordinates[j]>=count[j]){
							coordinates[j]=0;
							if(j<dim)
								projectedCoordinates[j]=0;
							//do nothing to projectedCoordinates when j==dim
							else if(j>dim)
								projectedCoordinates[j-1]=0;
						}
						else
							break;
					}
				}
				return(h);
			}
			
			///Extract a histogram with one fewer dimension by dropping all bins in dimension dim except bin
			histogram<N-1,StoreType> extractSlice(unsigned int dim, unsigned int bin) const{
				static_assert(N>1,"Extracting a slice is only defined for histograms with more than one dimension");
				assert(dim<N);
				histogram<N-1,StoreType> h;
				//copy all of the axes except the dim'th into h
				for(unsigned int i=0,j=0;i<N; i++){
					if(i==dim)
						continue;
					h.setAxis(j,getAxis(i)->copy());
					j++; //note that this is skipped on the i==dim iteration
				}
				//stretch h to have the correct range
				for(unsigned int i=0,j=0;i<N; i++){
					if(i==dim)
						continue;
					//strange little sidestep here to avoid appendBins messing up the axis range
					h.appendBins(j,1);
					h.appendBins(j,getBinCount(i)-1);
					j++; //note that this is skipped on the i==dim iteration
				}
				h.setUseContentScaling(false);
				//actually copy over the data, removing all but one bin in dimension dim
				internalCoordinate coordinates[N]={0/*All zeros*/};
				internalCoordinate projectedCoordinates[N-1]={0/*All zeros*/};
				const internalCoordinate n=std::accumulate(count,count+N,1U,std::multiplies<internalCoordinate>());
				for(internalCoordinate i=0; i<n; i++){
					/*std::cout << "coordinates: (";
					for(unsigned int k=0; k<N; k++)
						std::cout << (k?",":"") << coordinates[k];
					std::cout << ") ";
					std::cout << "projected coordinates: (";
					for(unsigned int k=0; k<N-1; k++)
						std::cout << (k?",":"") << projectedCoordinates[k];
					std::cout << ")\n";*/
					if(coordinates[dim]==bin){
						//std::cout << " adding " << data[i] << std::endl;
						h.getBin(projectedCoordinates)+=data[i];
						//std::cout << "---\n" << h << "---\n";
					}
					
					for(unsigned int j=0; j<N; j++){
						coordinates[j]++;
						if(j<dim)
							projectedCoordinates[j]++;
						//do nothing to projectedCoordinates when j==dim
						else if(j>dim)
							projectedCoordinates[j-1]++;
						
						if(coordinates[j]>=count[j]){
							coordinates[j]=0;
							if(j<dim)
								projectedCoordinates[j]=0;
							//do nothing to projectedCoordinates when j==dim
							else if(j>dim)
								projectedCoordinates[j-1]=0;
						}
						else
							break;
					}
				}
				h.setUseContentScaling(this->useContentScaling);
				return(h);
			}
		};
		
		namespace detail{
			///A helper class which sets the runtime dimensions of a histogram only if doing so is possible
			template<int N, typename StoreType>
			struct conditionalDimensionSetter{
				static void setDimensions(histogram<N,StoreType>& h, unsigned int dim){ /*Do nothing*/ }
			};
			
			template<typename StoreType>
			struct conditionalDimensionSetter<Dynamic,StoreType>{
				static void setDimensions(histogram<Dynamic,StoreType>& h, unsigned int dim){ h.setDimensions(dim); }
			};
		}
				
		///Creates an empty duplicate of a histogram
		template<typename HistogramType>
		HistogramType makeEmptyHistogramCopy(const HistogramType& h){
			HistogramType n;
			using HistDimExt = detail::dimensionExtractor<histogram,HistogramType::dimensions,typename HistogramType::dataType>;
			const unsigned int actualDim=HistDimExt::getDim(&h);
			//In case the histogram dimension is dynamic it will need to be set
			detail::conditionalDimensionSetter<HistogramType::dimensions,typename HistogramType::dataType>::setDimensions(n,actualDim);
			//All of the axes need to be copied
			for(unsigned int i=0; i<actualDim; i++)
				n.setAxis(i, h.getAxis(i)->copy());
			return(n);
		}
		
		namespace detail{
			///A helper type which prints histograms
			template<int N, typename StoreType>
			struct histogramPrinter{
				histogramPrinter(std::ostream& os, const histogram<N,StoreType>& h){
					if(!h.initialized){
						os << std::endl;
						return;
					}
					const unsigned int dimensions=detail::dimensionExtractor<histogram,N,StoreType>::getDim(&h);
					const unsigned int n=std::accumulate(h.count,h.count+dimensions,1U,std::multiplies<unsigned int>());
					for(unsigned int i=0; i<n; i++){
						unsigned int fact=1;
						double scale=1.0;
						for(unsigned int j=0; j<dimensions; j++){
							unsigned int bin=(i/fact)%h.count[j];
							os << h.axes[j]->binEdge(bin) << ' ';
							if(h.useContentScaling)
								scale/=h.axes[j]->binWidth(bin);
							fact*=h.count[j];
						}
						os << scale*h.data[i] << std::endl;
						if((i%h.count[0])==(h.count[0]-1))
							os << std::endl;
					}
				}
			};
		} //namespace detail
		
		template<int N, typename StoreType>
		std::ostream& operator<<(std::ostream& os, const histogram<N,StoreType>& h){
			detail::histogramPrinter<N,StoreType>(os,h);
			return(os);
		}
		
		namespace detail{
			///\param e1 the first operand; must be ignored if x1 is false
			///\param e2 the second operand; must be ignored if x2 is false
			///\param x1 whether the first operand exists (is valid for use)
			///\param x2 whether the second operand exists (is valid for use)
			template<typename StoreType>
			StoreType divideEntries(const StoreType& e1, const StoreType& e2, bool x1, bool x2){
				if(!x2)
					return(std::numeric_limits<StoreType>::quiet_NaN());
				if(!x1){
					if(e2==0)
						return(std::numeric_limits<StoreType>::quiet_NaN());
					return(0);
				}
				//finally, the case where we can actually do some real division
				//TODO: how will this play together with content scaling?
				return(e1/e2);
			}
			template<typename StoreType>
			StoreType addEntries(const StoreType& e1, const StoreType& e2, bool x1, bool x2){
				if(!x1 && !x2)
					return(0);
				if(!x2)
					return(e1);
				if(!x1)
					return(e2);
				return(e1+e2);
			}
			template<typename StoreType>
			StoreType subtractEntries(const StoreType& e1, const StoreType& e2, bool x1, bool x2){
				if(!x1 && !x2)
					return(0);
				if(!x2)
					return(e1);
				if(!x1)
					return(StoreType(0)-e2);
				return(e1-e2);
			}
			template<typename StoreType>
			StoreType multiplyEntries(const StoreType& e1, const StoreType& e2, bool x1, bool x2){
				if(!x2 || !x1)
					return(0);
				return(e1*e2);
			}
		}
			
		//TODO: this fails messily if one histogram is truly empty (has never allocated storage)
		//TODO: underflow and overflow are treated, but simply, such that results will be likely
		//      be useless if the input histograms have different limits. Perhaps the user could
		//      be warned in this situation?
		template<int N, typename StoreType, typename OpType>
		histogram<N,StoreType> applyHistogramBinaryOperation(OpType op,
											  const histogram<N,StoreType>& h1,
											  const histogram<N,StoreType>& h2){
			const unsigned int dimensions=detail::dimensionExtractor<histogram,N,StoreType>::getDim(&h1);
			const unsigned int dimensions2=detail::dimensionExtractor<histogram,N,StoreType>::getDim(&h2);
			assert(dimensions==dimensions2 && "Histograms used as arguments to a binary operator must have the same dimensions");
			
			//If one input histogram is not initialized, return a copy of the other.
			//If neither is initialized we return a copy of the second, although this should not be relied upon.
			if(!h1.isInitialized())
				return(h2);
			if(!h2.isInitialized())
				return(h1);
			
			histogram<N,StoreType> r;
			detail::conditionalDimensionSetter<N,StoreType>::setDimensions(r,dimensions);
			//instantiate axes for all of r's dimesions and make sure that their limits are wide enough to include all data from h1 and h2
			for(unsigned int i=0;i<dimensions; i++){
				axis* a=h1.getAxis(i)->copy();
				if(h2.getAxis(i)->getLowerLimit()<a->getLowerLimit())
					a->setLowerLimit(h2.getAxis(i)->getLowerLimit());
				if(h2.getAxis(i)->getUpperLimit()>a->getUpperLimit())
					a->setUpperLimit(h2.getAxis(i)->getUpperLimit());
				r.setAxis(i,a);
			}
			
			//std::cout << "Stretching result histogram" << std::endl;
			//stretch r to actually hold all data
			//in the process, figure out whether h1 or h2 extends furth in each direction in each dimension
			
			typename histogram<N,StoreType>::externalCoordinate loc1[dimensions]; //location of the 'smallest' corner of the new histogram
			typename histogram<N,StoreType>::externalCoordinate loc2[dimensions]; //location of the 'largest' corner of the new histogram
			typename histogram<N,StoreType>::internalCoordinate edges[dimensions][2]; //the amount by which the larger histogram is larger on the given edge
			bool edgeDirs[dimensions][2]; //entry is true if h1 extends farther than h2 on that edge
			//zero, but of the right type to put in the result histogram
			auto zero=typename histogram<N,StoreType>::amount(decltype(std::declval<typename histogram<N,StoreType>::amount>().value)(0));
			for(unsigned int i=0;i<dimensions; i++){
				//DEBUG std::cout << "Dimension " << i << ":\n";
				
				if(h1.getBinCount(i)==0 && h2.getBinCount(i)==0){
					//DEBUG std::cout << " Both histograms empty\n";
					//both histograms contain nothing, nothing to do
					return(r);
				}
				else if(h1.getBinCount(i)==0){
					//DEBUG std::cout << " Histogram 1 empty\n";
					edgeDirs[i][0]=edgeDirs[i][1]=false;
					loc1[i]=h2.getBinCenter(i,0);
					loc2[i]=h2.getBinCenter(i,h2.getBinCount(i)-1);
					edges[i][0]=h2.getBinCount(i);
					edges[i][1]=0;
				}
				else if(h2.getBinCount(i)==0){
					//DEBUG std::cout << " Histogram 2 empty\n";
					edgeDirs[i][0]=edgeDirs[i][1]=true;
					loc1[i]=h1.getBinCenter(i,0);
					loc2[i]=h1.getBinCenter(i,h1.getBinCount(i)-1);
					edges[i][0]=h1.getBinCount(i);
					edges[i][1]=0;
				}
				else{
					//DEBUG std::cout << " Both histograms non-trivial\n";
					std::pair<axis::lookupResult,axis::internalCoordinate> extent;
					//left
					edgeDirs[i][0]=h1.range(i).first<h2.range(i).first;
					//DEBUG std::cout << " edgeDirs[" << i << "][0]=" << edgeDirs[i][0] << std::endl;
					if(edgeDirs[i][0]){
						loc1[i]=h1.getBinCenter(i,0);
						extent=h2.getAxis(i)->findBin(h1.getBinCenter(i,0));
					}
					else{
						loc1[i]=h2.getBinCenter(i,0);
						extent=h1.getAxis(i)->findBin(h2.getBinCenter(i,0));
					}

					switch(extent.first){
						case axis::FOUND_BIN:
							edges[i][0]=0;
							break;
						case axis::PREPEND_BINS:
						case axis::APPEND_BINS:
							edges[i][0]=extent.second;
							break;
						default:
							throw std::logic_error("Incompatible input histograms for binary operation");
					}

					//right
					edgeDirs[i][1]=h1.range(i).second>h2.range(i).second;
					//DEBUG std::cout << " edgeDirs[" << i << "][1]=" << edgeDirs[i][1] << std::endl;
					if(edgeDirs[i][1]){
						loc2[i]=h1.getBinCenter(i,h1.getBinCount(i)-1);
						extent=h2.getAxis(i)->findBin(h1.getBinCenter(i,h1.getBinCount(i)-1));
					}
					else{
						loc2[i]=h2.getBinCenter(i,h2.getBinCount(i)-1);
						extent=h1.getAxis(i)->findBin(h2.getBinCenter(i,h2.getBinCount(i)-1));
					}

					switch(extent.first){
						case axis::FOUND_BIN:
							edges[i][1]=0;
							break;
						case axis::PREPEND_BINS:
						case axis::APPEND_BINS:
							edges[i][1]=extent.second;
							break;
						default:
							throw std::logic_error("Incompatible input histograms for binary operation");
					}
					//DEBUG std::cout << " loc1[" << i << "]=" << loc1[i] << std::endl;
					//DEBUG std::cout << " loc2[" << i << "]=" << loc2[i] << std::endl;
					//DEBUG std::cout << " edges[" << i << "][0]=" << edges[i][0] << std::endl;
					//DEBUG std::cout << " edges[" << i << "][1]=" << edges[i][1] << std::endl;
				}
			}
			//insert corner points to stretch
			r.add(loc1,zero);
			r.add(loc2,zero);
			
			//std::cout << "Applying binary operation" << std::endl;
			bool cs=h1.getUseContentScaling();
			r.setUseContentScaling(cs);
			//compute the contents of r, paying attention to the boundary region where it may extend beyond h1 or h2, which are then implicitly zero
			typename histogram<N,StoreType>::const_iterator h1it=h1.begin(), h2it=h2.begin();
			for(typename histogram<N,StoreType>::iterator rit=r.begin(), end=r.end(); rit!=end; rit++){
				//determine whether the current bin is outside one or both of the input histograms
				bool outside1=false, outside2=false;
				for(unsigned int i=0;i<dimensions; i++){
					if(rit.getCoordinate(i)<edges[i][0])
						(edgeDirs[i][0]?outside2:outside1)=true;
					if(rit.getCoordinate(i)>=(r.getBinCount(i)-edges[i][1]))
						(edgeDirs[i][1]?outside2:outside1)=true;
					if(outside1 && outside2)
						break;
				}
				
				*rit=op((StoreType)*h1it,(StoreType)*h2it,!outside1,!outside2);
				
				//advance each input iterator only if we were at a point inside the corresponding histogram
				if(!outside1)
					h1it++;
				if(!outside2)
					h2it++;
			}
			
			r.setUnderflow(op(h1.getUnderflow(),h2.getUnderflow(),true,true));
			r.setOverflow(op(h1.getOverflow(),h2.getOverflow(),true,true));
			
			return(r);
		}
		
		///Add two histograms (binwise).
		///\pre h1 and h2 must have the same binning, although not necessarily the same extents
		template<int N, typename StoreType>
		histogram<N,StoreType> operator+(const histogram<N,StoreType>& h1, const histogram<N,StoreType>& h2){
			return(applyHistogramBinaryOperation(&detail::addEntries<StoreType>,h1,h2));
		}
		///Subtract two histograms (binwise).
		///\pre h1 and h2 must have the same binning, although not necessarily the same extents
		template<int N, typename StoreType>
		histogram<N,StoreType> operator-(const histogram<N,StoreType>& h1, const histogram<N,StoreType>& h2){
			return(applyHistogramBinaryOperation(&detail::subtractEntries<StoreType>,h1,h2));
		}
		///Multiply two histograms (binwise).
		///\pre h1 and h2 must have the same binning, although not necessarily the same extents
		template<int N, typename StoreType>
		histogram<N,StoreType> operator*(const histogram<N,StoreType>& h1, const histogram<N,StoreType>& h2){
			return(applyHistogramBinaryOperation(&detail::multiplyEntries<StoreType>,h1,h2));
		}
		///Divide two histograms (binwise).
		///\pre h1 and h2 must have the same binning, although not necessarily the same extents
		template<int N, typename StoreType>
		histogram<N,StoreType> operator/(const histogram<N,StoreType>& h1, const histogram<N,StoreType>& h2){
			return(applyHistogramBinaryOperation(&detail::divideEntries<StoreType>,h1,h2));
		}
		
		///Multiply a histogram by a scalar
		template<int N, typename StoreType, typename U>
		histogram<N,StoreType> operator*(const histogram<N,StoreType>& h, U a){
			return(histogram<N,StoreType>(h)*=a);
		}
		///Divide a histogram by a scalar
		template<int N, typename StoreType, typename U>
		histogram<N,StoreType> operator/(const histogram<N,StoreType>& h, U a){
			return(histogram<N,StoreType>(h)/=a);
		}
		///Multiply a histogram by a scalar
		template<int N, typename StoreType, typename U>
		histogram<N,StoreType> operator*(U a, const histogram<N,StoreType>& h){
			return(histogram<N,StoreType>(h)*=a);
		}
	}
}

#include "dynamic_histogram.h"

#endif
