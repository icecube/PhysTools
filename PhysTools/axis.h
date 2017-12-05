///\file axis.h
///This file defines the base type for histogram axes, and also the built-in concrete axis types

#ifndef PHYSTOOLS_AXIS_H
#define PHYSTOOLS_AXIS_H

#include <vector>

#include "PhysTools/detail/advance_copy.h"
#include "PhysTools/hdf5_serialization.h"

namespace phys_tools{
	namespace histograms{
		using namespace phys_tools::detail;
	
		namespace detail{ class histogram_access; }
		
		///\brief The base type for all histogram axes.
		///The histogram just stores a block of data, and relies on axis objects to define what the
		///entries in the block mean by giving the mapping from external (user) coordinates to
		///internal bin indices.
		class axis{
		public:
			typedef int internalCoordinate;
			typedef double externalCoordinate;
		protected:
			///the base point at which binning begins
			externalCoordinate base;
			///the stride which defines the bin size
			externalCoordinate stride;
			///the offset, in bins, of the current base point from the fundamental base point (base)
			internalCoordinate offset;
			///the current number of bins
			internalCoordinate count;
			///the minimum external coordinate values which does not underflow
			externalCoordinate lowerLimit;
			///the maximum external coordinate values which does not overflow
			externalCoordinate upperLimit;
			
			///constructor for use by derived classes
			axis(externalCoordinate b, externalCoordinate s, internalCoordinate o=0, internalCoordinate c=0,
				 externalCoordinate lowerLimit=-std::numeric_limits<externalCoordinate>::infinity(),
				 externalCoordinate upperLimit=std::numeric_limits<externalCoordinate>::infinity()):
			base(b),stride(s),offset(o),count(c),lowerLimit(lowerLimit),upperLimit(upperLimit){}
		public:
			///describes qualitatively the result of transforming an external coordinate to an internal coordinate
			enum lookupResult{FOUND_BIN,PREPEND_BINS,APPEND_BINS,UNDERFLOW_,OVERFLOW_};
			
			virtual ~axis(){}
			
			///All derived classes are expected to implement this polymorphic copy function to facilitate copying histograms
			///\return an axis object allocated with new, owned by the caller
			virtual axis* copy() const=0;
			
			///get the number of bins currently used on this axis
			virtual internalCoordinate numBins() const{
				return(count);
			}
			///performs the conversion of an external coordinate to an internal coordinate
			///\return a pair which tells whether a bin was found, an overflow or underflow occurred,
			///        or bins must be added, and in which bin the given external coordinate fell.
			///        If the first element is FOUND_BIN, the second is simply the bin index (internal coordinate).
			///        If the first element is PREPEND_BINS (APPEND_BINS) the second is the number of bins which
			///        must be prepended (appended); in this case the bin in which the external coordinate falls
			///        is always the newly added first (last) bin.
			///        If the fist element is UNDERFLOW_ or UNDERFLOW_ the second element is unspecified.
			virtual std::pair<lookupResult,internalCoordinate> findBin(externalCoordinate x) const=0;
			///gets the smallest external coordinate which falls in the bin with the given internal coordinate
			virtual externalCoordinate binEdge(internalCoordinate bin) const=0;
			///gets the range of external coordinates which fall in the bin with the given internal coordinate
			virtual externalCoordinate binWidth(internalCoordinate bin) const=0;
			///gets the central external coordinate which falls in the bin with the given internal coordinate
			virtual externalCoordinate binCenter(internalCoordinate bin) const=0;
			
			///gets guidance from the axis about whether the entries in histogram bins should be divided by the bin volumes.
			///For non-linear axes this scaling is often desirable in order to display data densities intuitively.
			virtual bool adviseContentScaling() const=0;
			
		protected:
			//these functions are intended for use only by the histogram class, via histogram_access
			friend class detail::histogram_access;
			
			///Inserts bins on the 'left' side of the axis.
			///It is expected that this function will usually be called as a result of findBin()
			/// returning a PREPEND_BINS recommendation.
			virtual void prependBins(unsigned int bins)=0;
			///Inserts bins on the 'right' side of the axis.
			///It is expected that this function will usually be called as a result of findBin()
			/// returning a APPEND_BINS recommendation.
			virtual void appendBins(unsigned int bins)=0;
			///Moves the base point of this axis to the right by bins, which may be negative
			virtual void shiftBase(int bins)=0;
			
			void writeBasicPropertiesHDF5(hid_t dataset_id) const;
			
		public:
			///get the smallest external coordinate which is not considered to underflow
			virtual externalCoordinate getLowerLimit() const{
				return(lowerLimit);
			}
			///set the smallest external coordinate which is not considered to underflow to be lim
			virtual void setLowerLimit(externalCoordinate lim){
				if(lim>upperLimit)
					throw std::runtime_error("Attempt to set histogram axis lower limit above upper limit");
				if(count>0 && lim>(base+offset*stride))
					throw std::runtime_error("Attempt to set histogram axis lower limit above current base");
				lowerLimit=lim;
			}
			
			///get the smallest external coordinate which is not considered to overflow
			virtual externalCoordinate getUpperLimit() const{
				return(upperLimit);
			}
			///set the smallest external coordinate which is not considered to overflow to be lim
			virtual void setUpperLimit(externalCoordinate lim){
				if(lim<lowerLimit)
					throw std::runtime_error("Attempt to set histogram axis upper limit below lower limit");
				if(count>0 && lim<(base+stride*(offset+count)))
					throw std::runtime_error("Attempt to set histogram axis upper limit below extent");
				upperLimit=lim;
			}
			
			///Serializes this axis to an HDF5 file
			///
			///\param container the location within the HDF5 file within which the axis should be written
			///\param idx index number of the axis which should be appended as a suffix to its name
			virtual void writeHDF(hid_t container, unsigned int idx) const=0;
			///All subclasses must implement
			///\code
			///static axis* readHDF(hid_t axis_id);
			///\endcode
			///To read seialized axis objects back into memory.
			///This function must call H5Dclose on axis_id before returning.
		};
		
		axis* deserializeAxis(hid_t container, const std::string& name);
		
		namespace detail{
			class histogram_access{
			public:
				template<typename AxisType>
				static void prependBins(AxisType& a, unsigned int bins){ a.prependBins(bins); }
				template<typename AxisType>
				static void appendBins(AxisType& a, unsigned int bins){ a.appendBins(bins); }
				template<typename AxisType>
				static void shiftBase(AxisType& a, int bins){ a.shiftBase(bins); }
			};
		}
		
		///An axis which has constant size binning and allows bins to be added dynamically, as needed
		class LinearAxis : public axis{
		protected:
			LinearAxis(externalCoordinate b, externalCoordinate s, internalCoordinate o, internalCoordinate c, externalCoordinate l, externalCoordinate u):
			axis(b,s,o,c,l,u){}
		public:
			///\param b the 'base' point for the axis, a bin edge relative to which all other bins will be defined
			///\param s the 'stride' of the axis, the width of the bins
			LinearAxis(externalCoordinate b, externalCoordinate s):
			axis(b,s){}
			
			virtual ~LinearAxis(){}
			
			virtual axis* copy() const{
				return(new LinearAxis(base,stride,offset,0,lowerLimit,upperLimit));
			}
			
			virtual std::pair<lookupResult,internalCoordinate> findBin(externalCoordinate x) const{
				if(x<lowerLimit || (std::isinf(lowerLimit) && x==lowerLimit))
					return(std::make_pair(UNDERFLOW_,0));
				if(x>=upperLimit || std::isnan(x)) //call NaNs overflow
					return(std::make_pair(OVERFLOW_,0));
				externalCoordinate effectiveBase=base+offset*stride;
				if(x<effectiveBase)
					return(std::make_pair(PREPEND_BINS,ceil((effectiveBase-x)/stride)));
				internalCoordinate idx=(x-effectiveBase)/stride;
				if(idx>=count)
					return(std::make_pair(APPEND_BINS,idx-count+1));
				return(std::make_pair(FOUND_BIN,idx));
			}
			virtual externalCoordinate binEdge(internalCoordinate bin) const{
				externalCoordinate edge=base+(bin+offset)*stride;
				if(edge<lowerLimit)
					return(lowerLimit);
				return(edge);
			}
			virtual externalCoordinate binWidth(internalCoordinate bin) const{
				externalCoordinate edge=base+(bin+offset)*stride;
				if(edge+stride>upperLimit)
					return(upperLimit-edge);
				return(stride);
			}
			virtual externalCoordinate binCenter(internalCoordinate bin) const{
				return(binEdge(bin)+binWidth(bin)/2);
			}
			virtual bool adviseContentScaling() const{
				return(false);
			}
			virtual void writeHDF(hid_t container, unsigned int idx) const;
			static axis* readHDF(hid_t axis_id);
		protected:
			virtual void prependBins(unsigned int bins){
				offset-=bins;
				count+=bins;
			}
			virtual void appendBins(unsigned int bins){
				count+=bins;
			}
			virtual void shiftBase(int bins){
				offset+=bins;
			}
		};
		
		///An axis which has logarithmically sized binning and allows bins to be added dynamically, as needed.
		///Main axis parameters (base and stride, in particular) are considered logarithmic,
		///so to create an axis with its base at 1 and 4 bins per decade, use:
		/// LogarithmicAxis(0.0,.25)
		class LogarithmicAxis : public axis{
		protected:
			LogarithmicAxis(externalCoordinate b, externalCoordinate s, internalCoordinate o, internalCoordinate c, externalCoordinate l, externalCoordinate u):
			axis(b,s,o,c,l,u){}
		public:
			///\param b the 'base' point for the axis, a bin edge relative to which all other bins will be defined
			///\param s the 'stride' of the axis, the width of the bins
			LogarithmicAxis(externalCoordinate b, externalCoordinate s):
			axis(b,s,/*offset*/0,/*count*/0,0.0){}
			
			virtual ~LogarithmicAxis(){}
			
			virtual axis* copy() const{
				return(new LogarithmicAxis(base,stride,offset,0,lowerLimit,upperLimit));
			}
			
			virtual std::pair<lookupResult,internalCoordinate> findBin(externalCoordinate x) const{
				if(x<lowerLimit || (lowerLimit==0 && x==lowerLimit))
					return(std::make_pair(UNDERFLOW_,0));
				if(x>=upperLimit || std::isnan(x)) //call NaNs overflow
					return(std::make_pair(OVERFLOW_,0));
				x = log10(x);
				externalCoordinate effectiveBase=base+offset*stride;
				if(x<effectiveBase)
					return(std::make_pair(PREPEND_BINS,ceil((effectiveBase-x)/stride)));
				internalCoordinate idx=(x-effectiveBase)/stride;
				if(idx>=count)
					return(std::make_pair(APPEND_BINS,idx-count+1));
				return(std::make_pair(FOUND_BIN,idx));
			}
			virtual externalCoordinate binEdge(internalCoordinate bin) const{
				externalCoordinate edge=pow((externalCoordinate)10,base+(bin+offset)*stride);
				if(edge<lowerLimit)
					return(lowerLimit);
				return(edge);
			}
			virtual externalCoordinate binWidth(internalCoordinate bin) const{
				externalCoordinate lowerEdge=binEdge(bin);
				externalCoordinate upperEdge=pow((externalCoordinate)10,base+(bin+offset+1)*stride);
				if(upperEdge>upperLimit)
					upperEdge=upperLimit;
				return(upperEdge-lowerEdge);
			}
			virtual externalCoordinate binCenter(internalCoordinate bin) const{
				externalCoordinate lowerEdge=binEdge(bin);
				externalCoordinate upperEdge=pow((externalCoordinate)10,base+(bin+offset+1)*stride);
				if(upperEdge>upperLimit)
					upperEdge=upperLimit;
				return((lowerEdge+upperEdge)/2);
			}
			virtual bool adviseContentScaling() const{
				return(true);
			}
			virtual void writeHDF(hid_t container, unsigned int idx) const;
			static axis* readHDF(hid_t axis_id);
		protected:
			virtual void prependBins(unsigned int bins){
				offset-=bins;
				count+=bins;
			}
			virtual void appendBins(unsigned int bins){
				count+=bins;
			}
			virtual void shiftBase(int bins){
				base+=bins*stride;
			}
		public:
			///Set the lower limit allowed for this axis.
			///Note that this value is not treated logarithmically.
			virtual void setLowerLimit(externalCoordinate lim){
				if(lim<0.0)
					throw std::runtime_error("Attempt to set logarithmic histogram axis upper limit below zero");
				if(lim>upperLimit)
					throw std::runtime_error("Attempt to set histogram axis lower limit above upper limit");
				if(count>0 && log10(lim)>base+offset*stride)
					throw std::runtime_error("Attempt to set histogram axis lower limit above current base");
				lowerLimit=lim;
			}
			///Set the upper limit allowed for this axis.
			///Note that this value is not treated logarithmically.
			virtual void setUpperLimit(externalCoordinate lim){
				if(lim<lowerLimit)
					throw std::runtime_error("Attempt to set histogram axis upper limit below lower limit");
				if(count>0 && log10(lim)<(base+stride*(count+offset)))
					throw std::runtime_error("Attempt to set histogram axis upper limit below extent");
				upperLimit=lim;
			}
		};
		
		///An axis which has constant size binning and a fixed number (and range) of bins
		class FixedLinearAxis : public LinearAxis{
		private:
			double max;
			unsigned int totalCount;
		protected:
			FixedLinearAxis(externalCoordinate b, externalCoordinate s, internalCoordinate o, internalCoordinate c, externalCoordinate l, externalCoordinate u):
			LinearAxis(b,s,o,c,l,u){}
		public:
			FixedLinearAxis(double min, double max, unsigned int count):
			//use of nextafter is critical to avoid allowing a bin containing only max
			LinearAxis(min, (max-min)/count, 0, 0, min, nextafter(max,-std::numeric_limits<double>::infinity())),
			max(max),totalCount(count){}
			
			virtual axis* copy() const{
				return(new FixedLinearAxis(lowerLimit,max,totalCount));
			}
			
			virtual void setLowerLimit(externalCoordinate lim){
				throw std::logic_error("The limits of a FixedLinearAxis cannot be altered after construction");
			}
			virtual void setUpperLimit(externalCoordinate lim){
				throw std::logic_error("The limits of a FixedLinearAxis cannot be altered after construction");
			}
			virtual void writeHDF(hid_t container, unsigned int idx) const;
			static axis* readHDF(hid_t axis_id);
		};
		
		///An axis which has logarithmic size binning and a fixed number (and range) of bins
		class FixedLogarithmicAxis : public LogarithmicAxis{
		private:
			double max;
			unsigned int totalCount;
		protected:
			FixedLogarithmicAxis(externalCoordinate b, externalCoordinate s, internalCoordinate o, internalCoordinate c, externalCoordinate l, externalCoordinate u):
			LogarithmicAxis(b,s,o,c,l,u){}
		public:
			FixedLogarithmicAxis(double min, double max, unsigned int count):
			//use of nextafter is critical to avoid allowing a bin containing only max
			LogarithmicAxis(log10(min), (log10(max)-log10(min))/count, 0, 0, min, nextafter(max,-std::numeric_limits<double>::infinity())),
			max(max),totalCount(count){}
			
			virtual axis* copy() const{
				return(new FixedLogarithmicAxis(lowerLimit,max,totalCount));
			}
			
			virtual void setLowerLimit(externalCoordinate lim){
				throw std::logic_error("The limits of a FixedLogarithmicAxis cannot be altered after construction");
			}
			virtual void setUpperLimit(externalCoordinate lim){
				throw std::logic_error("The limits of a FixedLogarithmicAxis cannot be altered after construction");
			}
			virtual void writeHDF(hid_t container, unsigned int idx) const;
			static axis* readHDF(hid_t axis_id);
		};
		
		///An axis with arbitrary, fixed binning, specified by the user as a set of bin edges.
		///Note that to define n contiguous bins, n+1 edges are required.
		///
		///Usage example:
		///\code
		/// histogram<1> h(FixedUserDefinedAxis{0,1,3,4});
		/// h.setUseContentScaling(false);
		/// h.add(.5); h.add(1.5); h.add(2.5); h.add(3.5);
		/// std::cout << h << std::endl;
		///\endcode
		///which prints:
		///\verbatim
		/// 0 1
		/// 1 2
		/// 3 1
		///\endverbatim
		///Note the use of setUseContentScaling to supress division of the middle bin's content
		///by its greater width, without which all bins would have been displayed with equal
		///contents, reflecting the uniformity of the input distribution.
		class FixedUserDefinedAxis : public axis{
		private:
			std::vector<externalCoordinate> edges;
			
		public:
			///Construct an axis with a set of edges contained in an iterator range
			///of sorted edges.
			///(This is not terribly efficient for non-RandomAccess iterators, as it
			///will make several passes through them.)
			template<typename ForwardIterator>
			FixedUserDefinedAxis(ForwardIterator edgesStart, ForwardIterator edgesEnd):
			axis(*edgesStart,1.0,0,0,*edgesStart,*advance_copy(edgesStart,std::distance(edgesStart,edgesEnd)-1)),
			edges(edgesStart,edgesEnd){
				assert(std::is_sorted(edgesStart,edgesEnd) && "Edges must be sorted!");
			}
			
			///Construct an axis with a set of edges contained in an initializer_list
			///of sorted edges.
			FixedUserDefinedAxis(std::initializer_list<externalCoordinate> edges):
			axis(*edges.begin(),1.0,0,0,*edges.begin(),*(edges.end()-1)),
			edges(edges){
				assert(std::is_sorted(edges.begin(),edges.end()) && "Edges must be sorted!");
			}
			
			FixedUserDefinedAxis(const axis& other):edges(other.numBins()+1),
			axis(other.binEdge(0),1.0,0,other.numBins(),other.binEdge(0),other.binEdge(other.numBins())){
				for(unsigned int i=0; i<=other.numBins(); i++)
					edges[i]=other.binEdge(i);
			}
			
			virtual axis* copy() const{
				return(new FixedUserDefinedAxis(edges.begin(),edges.end()));
			}
			
			virtual std::pair<lookupResult,internalCoordinate> findBin(externalCoordinate x) const{
				if(x<lowerLimit || (std::isinf(lowerLimit) && x==lowerLimit))
					return(std::make_pair(UNDERFLOW_,0));
				if(x>=upperLimit || std::isnan(x)) //call NaNs overflow
					return(std::make_pair(OVERFLOW_,0));
				auto it=std::upper_bound(edges.begin(),edges.end(),x);
				assert(it!=edges.begin() && it!=edges.end() && "Internal Error: unable to find correct bin");
				//this is the true bin index
				internalCoordinate idx=std::distance(edges.begin(),it)-1;
				//now twiddle the bin index if necessary to correct for the current 'view' of the true bins
				idx-=offset;
				if(idx<0)
					return(std::make_pair(PREPEND_BINS,-idx));
				else if(idx>=count)
					return(std::make_pair(APPEND_BINS,idx-count+1));
				return(std::make_pair(FOUND_BIN,idx));
			}
			virtual externalCoordinate binEdge(internalCoordinate bin) const{
				externalCoordinate edge=edges[bin+offset];
				if(edge<lowerLimit)
					return(lowerLimit);
				return(edge);
			}
			virtual externalCoordinate binWidth(internalCoordinate bin) const{
				externalCoordinate lowerEdge=binEdge(bin+offset);
				externalCoordinate upperEdge=edges[bin+offset+1];
				if(upperEdge>upperLimit)
					upperEdge=upperLimit;
				return(upperEdge-lowerEdge);
			}
			virtual externalCoordinate binCenter(internalCoordinate bin) const{
				externalCoordinate lowerEdge=binEdge(bin+offset);
				externalCoordinate upperEdge=edges[bin+offset+1];
				if(upperEdge>upperLimit)
					upperEdge=upperLimit;
				return((lowerEdge+upperEdge)/2);
			}
			virtual bool adviseContentScaling() const{
				//TODO: this answer is conservative, but not always desirable. Allow user to specify on construction?
				return(true);
			}
			virtual void prependBins(unsigned int bins){
				offset-=bins;
				assert(offset>=0 && "Internal Error: offset should never become negative");
				count+=bins;
			}
			virtual void appendBins(unsigned int bins){
				count+=bins;
			}
			virtual void shiftBase(int bins){
				offset+=bins;
				assert(offset>=0 && "Internal Error: offset should never become negative");
			}
			
			virtual void setLowerLimit(externalCoordinate lim){
				throw std::logic_error("The limits of a FixedUserDefinedAxis cannot be altered after construction");
			}
			virtual void setUpperLimit(externalCoordinate lim){
				throw std::logic_error("The limits of a FixedUserDefinedAxis cannot be altered after construction");
			}
			virtual void writeHDF(hid_t container, unsigned int idx) const;
			static axis* readHDF(hid_t axis_id);
		};
		
	} //namespace histograms
} //namespace phys_tools

#endif