#ifndef PHYSTOOLS_DYNAMIC_HISTOGRAM_H
#define PHYSTOOLS_DYNAMIC_HISTOGRAM_H

namespace phys_tools{
	namespace histograms{
///\brief A specialization of the histogram template for histograms whose dimensions are not known until runtime.
///
///\warning This code is incomplete and largely untested!
template<typename StoreType>
class histogram<Dynamic,StoreType> : public histogramBase<Dynamic,StoreType,unsigned int,double,typename detail::histogramTraits<StoreType>::amount>{
public:
	typedef StoreType dataType;
	typedef unsigned int internalCoordinate;
	typedef double externalCoordinate;
	typedef typename detail::histogramTraits<StoreType>::amount amount;
	typedef histogramBase<Dynamic,StoreType,internalCoordinate,externalCoordinate,amount> BaseType;
	
	using BaseType::initialized;
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
	using typename BaseType::modificationProxy;
	
	//the base needs access to private members
	friend BaseType;
	//all histograms are friends
	template<int M, typename OtherStoreType>
	friend class histogram;
	
	template <template<int,class> class H, int D, class S>
	friend struct detail::dimensionExtractor;
	
	template<int N, typename OtherStoreType, typename OpType>
	friend histogram<N,StoreType> applyHistogramBinaryOperation(OpType op,const histogram<N,OtherStoreType>& h1,const histogram<N,OtherStoreType>& h2);
private:
	unsigned int dimensions; //stores the number of dimensions this histgrom instance actually has, since the type is just 'Dynamic' (=-1)
	axis** axes;
	internalCoordinate* count;
	StoreType* data;
	StoreType underflow, overflow;
	
public:
	///Construct a histogram without axes or a specific number of dimensions.
	///\post nothing else may be done with the newly constructed histogram until setDimensions has been called
	histogram():
	dimensions(0),
	axes(NULL),
	count(NULL),
	data(new StoreType[1]),
	underflow(0),overflow(0)
	{}
	
	///Construct a histogram without axes
	///\param dim the number of dimensions the histogram will have
	///\pre dim>0
	histogram(unsigned int dim):
	dimensions(dim),
	axes(new axis*[dim]),
	count(new internalCoordinate[dim]),
	data(new StoreType[1]),
	underflow(0),overflow(0)
	{
		std::fill_n(&count[0],dim,0);
		std::fill_n(&axes[0],dim,nullptr);
		assert(dim>0 && "Histograms must have at least one dimension");
	}
	
	///Construct a histogram from a set of axes.
	///The histogram's dimensionality will be assumed from the number of axes supplied
	template<typename... Axes, typename = typename std::enable_if<not first_arg_convertible<unsigned int,Axes...>::value>::type>
	histogram(Axes... axes):
	dimensions(sizeof...(axes)),
	axes(new axis*[dimensions]),
	count(new internalCoordinate[dimensions]),
	data(new StoreType[1]),
	underflow(0),overflow(0){
		static_assert(all_args_convertible<axis&,typename std::add_lvalue_reference<Axes>::type...>::value,
					  "all arguments to histogram(Axes... axes) "
					  "must be axis objects");
		std::fill_n(&count[0],dimensions,0);
		std::fill_n(&this->axes[0],dimensions,nullptr);
		setAxes_impl<0>(axes...);
	}
	
	///Copy construct a histogram from another of the same type.
	///All axes and data are copied.
	histogram(const histogram& other):
	BaseType(other.initialized),
	dimensions(other.dimensions),
	axes(new axis*[dimensions]),
	count(new internalCoordinate[dimensions]),
	data(new StoreType[std::accumulate(other.count,other.count+other.dimensions,1U,std::multiplies<unsigned int>())]),
	underflow(other.underflow),
	overflow(other.overflow){
		std::copy(&other.count[0],&other.count[dimensions],&count[0]);
		unsigned long totalSize=1;
		for(unsigned int i=0; i<dimensions; i++){
			axes[i]=other.axes[i]->copy();
			totalSize*=count[i];
			if(count[i]>0)
				detail::histogram_access::appendBins(*axes[i],count[i]);
		}
		std::copy(&other.data[0],&other.data[totalSize],&data[0]);
	}
	
	///Move construct a histogram from another of the same type.
	///All axes and data are taken from the other histogram and moved into the newly
	///constructed object, leaving the original empty and unusable.
	histogram(histogram&& other):
	BaseType(other.initialized,other.useContentScaling),
	dimensions(other.dimensions),
	axes(other.axes),
	count(other.count),
	data(other.data),
	underflow(other.underflow),
	overflow(other.overflow){
		other.axes=nullptr;
		other.count=nullptr;
		other.data=new StoreType[1];
		other.initialized=false;
		other.underflow=StoreType();
		other.overflow=StoreType();
		other.count=nullptr;
		other.axes=nullptr;
	}
	
	~histogram(){
		if(axes)
			for(unsigned int i=0; i<dimensions; i++)
				delete axes[i];
		delete[] axes;
		delete[] count;
		delete[] data;
	}
	
	///Copy from another histogram of the same type.
	///All axes and data are copied.
	histogram& operator=(const histogram& other){
		if(&other==this)
			return(*this);
		for(unsigned int i=0; i<dimensions; i++)
			delete axes[i];
		delete[] axes;
		delete[] count;
		delete[] data;
		dimensions=other.dimensions;
		initialized=other.initialized;
		underflow=other.underflow;
		overflow=other.overflow;
		axes=new axis*[dimensions];
		count=new internalCoordinate[dimensions];
		std::copy(&other.count[0],&other.count[dimensions],&count[0]);
		for(unsigned int i=0; i<dimensions; i++){
			axes[i]=other.axes[i]->copy();
			if(count[i]>0)
				axes[i]->appendBins(count[i]);
		}
		unsigned long totalSize=std::accumulate(other.count,other.count+dimensions,1U,std::multiplies<unsigned int>());
		data=new StoreType[totalSize];
		std::copy(&other.data[0],&other.data[totalSize],&data[0]);
		return(*this);
	}
	
	///Move from another histogram of the same type.
	///All axes and data are taken from the other histogram and moved into this
	///histogram, leaving the other empty an unusable.
	histogram& operator=(histogram&& other){
		if(&other==this)
			return(*this);
		for(unsigned int i=0; i<dimensions; i++)
			delete axes[i];
		delete[] axes;
		delete[] count;
		delete[] data;
		dimensions=other.dimensions;
		data=other.data;
		other.data=new StoreType[1];
		initialized=other.initialized;
		other.initialized=false;
		underflow=other.underflow;
		other.underflow=StoreType(0);
		overflow=other.overflow;
		other.overflow=StoreType(0);
		count=other.count;
		other.count=nullptr;
		axes=other.axes;
		other.axes=nullptr;
		return(*this);
	}
	
	///Set the number of dimensions for a freshly default-constructed histogram
	///\param dim the number of dimensions the histogram will have
	///\pre dim>0
	///\pre the histogram was constructed with the default constructor and setDimensions has not been previously called
	void setDimensions(unsigned int dim){
		assert(!initialized && "Once initialized, a histogram's number of dimensions must not be changed");
		assert(dim>0 && "Histograms must have at least one dimension");
		for(unsigned int i=0; i<dimensions; i++)
			delete axes[i];
		delete[] axes;
		axes=new axis*[dim];
		delete[] count;
		count=new internalCoordinate[dim];
		dimensions=dim;
	}
	
	unsigned int getDimensions() const{
		return(dimensions);
	}
	
private:
	//internal functions for setting axes
	template<unsigned int index>
	void setAxes_impl() const{
		assert(index==dimensions && "incorrect number of arguments to histogram<Dynamic>::setAxes()");
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
	//internal machinery for adding and accessing data
	
	//base case for default amount (1)
	template<unsigned int index>
	void addImpl(externalCoordinate coordinates[]){
		assert(index==dimensions && "incorrect number of arguments to histogram<>::add()");
		addCore(coordinates,detail::histogramTraits<StoreType>::unit());
	}
	//base case for specified amount
	template<unsigned int index>
	void addImpl(externalCoordinate coordinates[], amount a){
		assert(index==dimensions && "incorrect number of arguments to histogram<>::add()");
		addCore(coordinates,a);
	}
	template<unsigned int index,typename... Args>
	void addImpl(externalCoordinate coordinates[], externalCoordinate nextCoord, Args... args){
		assert(index<dimensions && "incorrect number of arguments to histogram<>::add()");
		coordinates[index]=nextCoord;
		addImpl<index+1>(coordinates,args...);
	}
	
	template<unsigned int index>
	StoreType subscriptImpl(internalCoordinate* coordinates, double scale) const{
		StoreType temp=getBin(coordinates);
		return(scale*temp);
	}
	template<unsigned int index,typename... Args>
	StoreType subscriptImpl(internalCoordinate* coordinates, double scale, internalCoordinate nextCoord, Args... args) const{
		coordinates[index]=nextCoord;
		StoreType temp=subscriptImpl<index+1>(coordinates,scale,args...);
		return(temp);
	}
	
	//non-const for raw access
	template<unsigned int index>
	modificationProxy subscriptImplRaw(const internalCoordinate* coordinates, double scale){
		return(modificationProxy(getBin(coordinates),scale));
	}
	template<unsigned int index,typename... Args>
	modificationProxy subscriptImplRaw(internalCoordinate* coordinates, double scale, internalCoordinate nextCoord, Args... args){
		coordinates[index]=nextCoord;
		return(subscriptImplRaw<index+1>(coordinates,scale,args...));
	}
	
	///Change the rank of the histogram to match a serialized state being read from a file
	void mismatchedDimensionsOnRead(const std::string&, size_t dim){
		setDimensions(dim);
	}
	
public:
	///Set the axes for all dimensions
	template<typename... Axes>
	void setAxes(Axes... axes){
		static_assert(sizeof...(axes)==dimensions,
					  "incorrect number of arguments to histogram<Dynamic>::setAxes()");
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
	
	///\brief Insert a count into the histogram.
	///
	///The first N arguments must be the (external) coordinates at which the count is to be inserted,
	///and the optional last argument is the amount to insert, if not a unit count.
	///
	///Usage example:
	///\code
	/// histogram<Dynamic> h(LinearAxis(0,1),LinearAxis(0,1));
	/// h.add(1,3); //add a unit count at (1,3)
	/// h.add(2,4,histogram<Dynamic>::amount(5)); //add 5 units at (2,4)
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
		static_assert(all_args_convertible<externalCoordinate,Args...>::value ||
					  (all_args_except_last_convertible<externalCoordinate,Args...>::value && last_arg_convertible<amount,Args...>::value),
					  "incorrect type of arguments to histogram<>::add(): all arguments must be external coordinates, except the optional last which may be an amount");
		
		assert((sizeof...(args)==dimensions || sizeof...(args)==dimensions+1) && "incorrect number of arguments to histogram<Dynamic>::add()");
		assert(((sizeof...(args)==dimensions && all_args_convertible<externalCoordinate,Args...>::value) || sizeof...(args)==dimensions+1)
			   && "incorrect type of arguments to histogram<Dynamic>::add(); when using N arguments they must all be external coordinates");
		assert(((sizeof...(args)==dimensions+1 && all_args_except_last_convertible<externalCoordinate,Args...>::value) || sizeof...(args)==dimensions)
			   && "incorrect type of arguments to histogram<>::add(); when using N+1 arguments all but the last must be external coordinates");
		assert(((sizeof...(args)==dimensions+1 && last_arg_convertible<amount,Args...>::value) || sizeof...(args)==dimensions)
			   && "incorrect type of arguments to histogram<>::add(); when using N+1 arguments the last must be an amount");
		
		externalCoordinate coordinates[dimensions];
		addImpl<0>(coordinates,args...);
	}
	
	///Add a count to the histogram at a set of external coordinates passed in an array.
	void add(const externalCoordinate* coordinates, amount a=1.0){
		addCore(coordinates,a);
	}
	
	///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates
	template<typename... Args, typename = typename std::enable_if<not first_is_pointer<Args...>::value>::type>
	StoreType operator()(Args... args) const{
		assert(sizeof...(args)==dimensions && "incorrect number of arguments to histogram<Dynamic>::operator()");
		static_assert(all_args_convertible<internalCoordinate,Args...>::value,"all arguments to histogram<Dynamic>::operator() must be internal coordinates");
		unsigned int coordinates[dimensions];
		StoreType temp=subscriptImpl<0>(coordinates,1.0,args...);
		return(temp);
	}
	///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates, passed in an array
	StoreType operator()(const internalCoordinate* coordinates) const{
		double scale=1.0;
		return(getBin(coordinates)*scale);
	}
	
	///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates
	template<typename... Args, typename = typename std::enable_if<not first_is_pointer<Args...>::value>::type>
	modificationProxy operator()(Args... args){
		assert(sizeof...(args)==dimensions && "incorrect number of arguments to histogram<Dynamic>::operator()");
		static_assert(all_args_convertible<internalCoordinate,Args...>::value,"all arguments to histogram<Dynamic>::operator() must be internal coordinates");
		unsigned int coordinates[dimensions];
		return(subscriptImplRaw<0>(coordinates,1.0,args...));
	}
	///Get the value stored in the bin corresponding to the passed set of N (internal) coordinates, passed in an array
	modificationProxy operator()(const internalCoordinate* coordinates){
		double scale=1.0;
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
		internalCoordinate* coordinates;
		T value;
		bool isEnd;
		
		iteratorTempl(HistType& h, bool e):
		h(&h),
		coordinates(new internalCoordinate[h.dimensions]),
		isEnd(e){
			std::fill_n(&coordinates[0],h.dimensions,0);
		}
		
		void fetchValue(){
			if(!isEnd)
				value=(*h)(coordinates);
		}
		
		friend HistType;
		static iteratorTempl make_begin(HistType& h){
			iteratorTempl it(h,std::accumulate(h.count,h.count+h.getDimensions(),1U,std::multiplies<internalCoordinate>())==0);
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
		h(it.h),coordinates(new internalCoordinate[h->dimensions]),
		value(it.value),isEnd(it.isEnd){
			std::copy(&it.coordinates[0],&it.coordinates[0]+h->dimensions,&coordinates[0]);
		}
		
		iteratorTempl(iteratorTempl&& it):
		h(it.h),coordinates(it.coordinates),value(it.value),isEnd(it.isEnd){
			it.coordinates=nullptr;
		}
		
		~iteratorTempl(){
			delete[] coordinates;
		}
		
		iteratorTempl& operator=(const iteratorTempl& it){
			if(this==&it)
				return(*this);
			h=it.h;
			isEnd=it.isEnd;
			delete[] coordinates;
			coordinates=new internalCoordinate[h->dimensions];
			std::copy(&it.coordinates[0],&it.coordinates[0]+h->dimensions,&coordinates[0]);
			if(!isEnd)
				value=it.value;
			else
				coordinates=nullptr;
			return(*this);
		}
		
		iteratorTempl& operator=(iteratorTempl&& it){
			if(this==&it)
				return(*this);
			h=it.h;
			isEnd=it.isEnd;
			coordinates=it.coordinates;
			it.coordinates=nullptr;
			if(!isEnd)
				value=it.value;
			return(*this);
		}
		
		///prefix increment
		iteratorTempl& operator++(){
			int i;
			const unsigned int N=h->dimensions;
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
			const unsigned int N=h->dimensions;
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
			for(int i=0; i<h->dimensions; i++){
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
			for(int i=0; i<h->dimensions; i++){
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
		const internalCoordinate* getCoordinates() const{
			return(coordinates);
		}
		///Get the internal coordinate in dimension dim within the histogram (bin index) to which this iterator refers
		internalCoordinate getCoordinate(unsigned int dim) const{
			assert(dim<h->dimensions);
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
		internalCoordinate* coordinates;
		T value;
		bool isEnd;
		
		rIteratorTempl(HistType& h, bool e):
		h(&h),
		coordinates(new internalCoordinate[h.dimensions]),
		isEnd(e){
			std::transform(&h.count[0],&h.count[0]+h.dimensions,&coordinates[0],
						   [](internalCoordinate c){ return(c-1); });
		}
		
		void fetchValue(){
			if(!isEnd)
				value=(*h)(coordinates);
		}
		
		friend HistType;
		static rIteratorTempl make_begin(HistType& h){
			rIteratorTempl it(h,std::accumulate(h.count,h.count+h.getDimensions(),1U,std::multiplies<internalCoordinate>())==0);
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
		h(it.h),coordinates(new internalCoordinate[h->dimensions]),
		value(it.value),isEnd(it.isEnd){
			std::copy(&it.coordinates[0],&it.coordinates[0]+h->dimensions,&coordinates[0]);
		}
		
		rIteratorTempl(rIteratorTempl&& it):
		h(it.h),coordinates(it.coordinates),value(it.value),isEnd(it.isEnd){
			it.coordinates=nullptr;
		}
		
		~rIteratorTempl(){
			delete[] coordinates;
		}
		
		rIteratorTempl& operator=(const rIteratorTempl& it){
			if(this==&it)
				return(*this);
			h=it.h;
			isEnd=it.isEnd;
			delete[] coordinates;
			coordinates=new internalCoordinate[h->dimensions];
			std::copy(&it.coordinates[0],&it.coordinates[0]+h->dimensions,&coordinates[0]);
			if(!isEnd)
				value=it.value;
			else
				coordinates=nullptr;
			return(*this);
		}
		
		rIteratorTempl& operator=(rIteratorTempl&& it){
			if(this==&it)
				return(*this);
			h=it.h;
			isEnd=it.isEnd;
			coordinates=it.coordinates;
			it.coordinates=nullptr;
			if(!isEnd)
				value=it.value;
			return(*this);
		}
		
		///prefix increment
		rIteratorTempl& operator++(){
			int i;
			const unsigned int N=h->dimensions;
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
			const unsigned int N=h->dimensions;
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
			if(h!=rhs.h) //otherwise, iterators to different histograms cannot be equal
				return(false);
			for(int i=0; i<h->dimensions; i++){
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
			if(h!=rhs.h) //otherwise, iterators to different histograms cannot be equal
				return(true);
			for(int i=0; i<h->dimensions; i++){
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
		const internalCoordinate* getCoordinates() const{
			return(coordinates);
		}
		///Get the internal coordinate in dimension dim within the histogram (bin index) to which this iterator refers
		internalCoordinate getCoordinate(unsigned int dim) const{
			assert(dim<h->dimensions);
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
	typedef iteratorTempl<histogram,modificationProxy,modificationProxy,modificationProxy*> iterator;
	///A const bidirectional iterator
	typedef iteratorTempl<const histogram,StoreType,StoreType&,const StoreType*> const_iterator;
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
	
public:
	
	///Get an iterator to the bin in which an entry at the given coordinates would fall,
	///or an end iterator if such an entry is outside the current extent of the histogram
	template<typename... Args, typename = typename std::enable_if<not first_is_pointer<Args...>::value>::type>
	iterator findBin(Args... coordinates){
		assert(sizeof...(coordinates)==dimensions && "incorrect number of arguments to histogram<Dynamic>::findBin()");
		iterator it(*this,false); //assume that it will not be the end iterator
		findBinImpl<0>(it,coordinates...);
		return(it);
	}
	
	///Get an iterator to the bin in which an entry at the given coordinates would fall,
	///or an end iterator if such an entry is outside the current extent of the histogram
	template<typename... Args, typename = typename std::enable_if<not first_is_pointer<Args...>::value>::type>
	const_iterator findBin(Args... coordinates) const{
		assert(sizeof...(coordinates)==dimensions && "incorrect number of arguments to histogram<Dynamic>::findBin()");
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
		for(auto bin : *this)
			bin*=a;
		underflow*=a;
		overflow*=a;
		return(*this);
	}
	///Divide all bin contents by a
	template<typename U>
	histogram& operator/=(U a){
		for(auto bin : *this)
			bin/=a;
		underflow/=a;
		overflow/=a;
		return(*this);
	}
	
	///Construct a new histogram with one fewer dimension with the projection of this histogram's data
	///\param dim the dimension to be projected out
	histogram<Dynamic,StoreType> project(unsigned int dim) const{
		assert(dim<dimensions);
		histogram<Dynamic,StoreType> h(dimensions-1);
		//copy all of the axes except the dim'th into h
		for(unsigned int i=0,j=0;i<dimensions; i++){
			if(i==dim)
				continue;
			h.setAxis(j,getAxis(i)->copy());
			j++; //note that this is skipped on the i==dim iteration
		}
		
		//stretch h to have the correct range
		for(unsigned int i=0,j=0;i<dimensions; i++){
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
		internalCoordinate coordinates[dimensions];
		std::fill_n(&coordinates[0],dimensions,0);
		internalCoordinate projectedCoordinates[dimensions-1];
		std::fill_n(&projectedCoordinates[0],dimensions-1,0);
		const internalCoordinate n=std::accumulate(count,count+dimensions,1U,std::multiplies<internalCoordinate>());
		for(internalCoordinate i=0; i<n; i++){
			//std::cout << "projected coordinates: (";
			//for(unsigned int k=0; k<N-1; k++)
			//	std::cout << (k?",":"") << projectedCoordinates[k];
			//std::cout << ")\n";
			h.getBin(projectedCoordinates)+=data[i];
			//std::cout << " adding " << data[i] << ", sum now " << h.getBin(projectedCoordinates) << std::endl;
			
			for(unsigned int j=0; j<dimensions; j++){
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
	histogram<Dynamic,StoreType> extractSlice(unsigned int dim, unsigned int bin) const{
		assert(dim<dimensions);
		histogram<Dynamic,StoreType> h(dimensions-1);
		//copy all of the axes except the dim'th into h
		for(unsigned int i=0,j=0;i<dimensions; i++){
			if(i==dim)
				continue;
			h.setAxis(j,getAxis(i)->copy());
			j++; //note that this is skipped on the i==dim iteration
		}
		//stretch h to have the correct range
		for(unsigned int i=0,j=0;i<dimensions; i++){
			if(i==dim)
				continue;
			//strange little sidestep here to avoid appendBins messing up the axis range
			h.appendBins(j,1);
			h.appendBins(j,getBinCount(i)-1);
			j++; //note that this is skipped on the i==dim iteration
		}
		//actually copy over the data, removing all but one bin in dimension dim
		internalCoordinate coordinates[dimensions];
		std::fill_n(&coordinates[0],dimensions,0);
		internalCoordinate projectedCoordinates[dimensions-1];
		std::fill_n(&projectedCoordinates[0],dimensions-1,0);
		const internalCoordinate n=std::accumulate(count,count+dimensions,1U,std::multiplies<internalCoordinate>());
		for(internalCoordinate i=0; i<n; i++){
			/*std::cout << "coordinates: (";
			for(unsigned int k=0; k<dimensions; k++)
				std::cout << (k?",":"") << coordinates[k];
			std::cout << ") ";
			std::cout << "projected coordinates: (";
			for(unsigned int k=0; k<dimensions-1; k++)
				std::cout << (k?",":"") << projectedCoordinates[k];
			std::cout << ")\n";*/
			if(coordinates[dim]==bin){
				//std::cout << " adding " << data[i] << std::endl;
				h.getBin(projectedCoordinates)+=data[i];
				//std::cout << "---\n" << h << "---\n";
			}
			
			for(unsigned int j=0; j<dimensions; j++){
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
	
	//enable printing!
	template<int, typename>
	friend struct detail::histogramPrinter;
};
} //namespace histograms
} //namespace phys_tools

#endif
