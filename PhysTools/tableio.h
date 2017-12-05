#ifndef PHYSTOOLS_TABLEIO_H
#define PHYSTOOLS_TABLEIO_H

#include <cassert>
#include <string>
#include <map>
#include <memory>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/string.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <PhysTools/hdf5_serialization.h>

namespace phys_tools{

//Slightly ugly: CTS is a macro, so it effectively lives in the global namespace, 
//but it expands to things which are in the phys_tools namespace, which may not
//be in scope. Putting these utilities in a nested namespace makes it easier for
//users to bring them into scope with `using namespace phys_tools::cts` rather
//than a much more drastic `using namespace phys_tools`.
namespace cts{
template <int N>
constexpr char at(char const(&s)[N], int i){
	return(i >= N ? '\0' : s[i]);
}

template <class S, char C, bool EOS>
struct push_back_if_c :
boost::mpl::push_back<S, boost::mpl::char_<C>>
{};

template <class S, char C>
struct push_back_if_c<S, C, true> : S {};

#ifndef CTS_STRING_MAX_LENGTH
#define CTS_STRING_MAX_LENGTH 48
#endif

#define CTS_PRE(z, n, u) push_back_if_c<
#define CTS_POST(z, n, s) , at(s, n), (n >= sizeof(s))>::type

#define CTS(s)                                   \
BOOST_PP_REPEAT(CTS_STRING_MAX_LENGTH, CTS_PRE, ~)    \
boost::mpl::string<>                       \
BOOST_PP_REPEAT(CTS_STRING_MAX_LENGTH, CTS_POST, s)

} //namespace cts

namespace tableio{
//Type is the type of the field, Name is an MPL sequence from which the
//field name can be derived by application of mpl::c_str
template<typename Type, typename Name, bool AutoConvert=true>
struct field{
	typedef Type type;
	typedef Name name; //note that this is a type, not a value!
	static constexpr bool autoConvert=AutoConvert;
	static std::string getName(){
		return(boost::mpl::c_str<Name>::value);
	}
	static hid_t getType(){
		return(hdf_interface::getHDFDatatype<type>().datatype);
	}
};

namespace detail{
	template<class... FieldTypes> struct hoboTuple;
	
	template<>
	struct hoboTuple<>{
		static const int nFields=0;
		
		static void getFieldSizes(size_t sizes[], unsigned int counter){}
		static void getFieldOffsets(size_t offsets[], unsigned int counter, size_t base){}
		static void getFieldNames(std::string names[], unsigned int counter){}
		static void getFieldTypes(hid_t types[], unsigned int counter){}
		static bool getAutoConvert(size_t){ return(false); }
	};
	
	template<unsigned int index, class FieldType, class... OtherFields>
	struct FieldType_{
		typedef typename FieldType_<index-1,OtherFields...>::type type;
	};
	
	template<class FieldType, class... OtherFields>
	struct FieldType_<0, FieldType, OtherFields...>{
		typedef typename FieldType::type type;
	};
	
	template<unsigned int index, class FieldType, class... OtherFields>
	struct getImpl{
		static typename FieldType_<index,FieldType,OtherFields...>::type get(const hoboTuple<FieldType,OtherFields...>& tup){
			return(getImpl<index-1,OtherFields...>::get(tup.suffix));
		}
		static typename FieldType_<index,FieldType,OtherFields...>::type& get(hoboTuple<FieldType,OtherFields...>& tup){
			return(getImpl<index-1,OtherFields...>::get(tup.suffix));
		}
	};
	
	template<class FieldType, class... OtherFields>
	struct getImpl<0,FieldType,OtherFields...>{
		static typename FieldType_<0,FieldType,OtherFields...>::type get(const hoboTuple<FieldType,OtherFields...>& tup){
			return(tup.field);
		}
		static typename FieldType_<0,FieldType,OtherFields...>::type& get(hoboTuple<FieldType,OtherFields...>& tup){
			return(tup.field);
		}
	};
	
	template<class FieldType, class... OtherFields>
	struct hoboTuple<FieldType,OtherFields...>{
		typedef hoboTuple<OtherFields...> SuffixType;
		
		typename FieldType::type field;
		SuffixType suffix;
		
		static const int nFields=SuffixType::nFields+1;
		
		hoboTuple() = default;
		
		template<typename Arg, typename... ArgTypes>
		hoboTuple(Arg arg, ArgTypes... args):field(arg),suffix(args...){}
		
		static void getFieldSizes(size_t sizes[], unsigned int counter){
			sizes[counter++]=sizeof(typename FieldType::type);
			SuffixType::getFieldSizes(sizes,counter);
		}
		
		static void getFieldOffsets(size_t offsets[], unsigned int counter, size_t base){
			typedef hoboTuple<FieldType,OtherFields...> ThisType;
			offsets[counter++]=base+offsetof(ThisType,field);
			SuffixType::getFieldOffsets(offsets,counter,base+offsetof(ThisType,suffix));
		}
		
		static void getFieldNames(std::string names[], unsigned int counter){
			names[counter++]=FieldType::getName();
			SuffixType::getFieldNames(names,counter);
		}

		static void getFieldTypes(hid_t types[], unsigned int counter){
			types[counter++]=FieldType::getType();
			SuffixType::getFieldTypes(types,counter);
		}

		static bool getAutoConvert(size_t index){
			if(index==0)
				return(FieldType::autoConvert);
			return(SuffixType::getAutoConvert(index-1));
		}
	};
	
	template<unsigned int index_, typename type_>
	struct IndexForNameHolder{
		static constexpr unsigned int index=index_;
		typedef type_ type;
	};
	
	//fwd decl
	template<unsigned int index_, typename Name, typename FieldType, typename... OtherFieldTypes>
	struct IndexForNameImpl;
	
	//base case
	template<unsigned int index_, typename Name, typename FieldType>
	struct IndexForNameImpl<index_,Name,FieldType>{
	private:
		static_assert(std::is_same<Name,typename FieldType::name>::value, "Name does not name a field in the current tuple");
		using underlying_type = typename std::conditional<std::is_same<Name,typename FieldType::name>::value,
		IndexForNameHolder<index_,typename FieldType::type>,
		void
		>::type;
	public:
		static constexpr unsigned int index=underlying_type::index;
		using type = typename underlying_type::type;
	};
	
	//recursive case
	template<unsigned int index_, typename Name, typename FieldType, typename... OtherFieldTypes>
	struct IndexForNameImpl{
	private:
		static_assert(std::is_same<Name,typename FieldType::name>::value || sizeof...(OtherFieldTypes)>0, "Name does not name a field in the current tuple");
		using underlying_type = typename std::conditional<std::is_same<Name,typename FieldType::name>::value,
		IndexForNameHolder<index_,typename FieldType::type>,
		IndexForNameImpl<index_+1,Name,OtherFieldTypes...>
		>::type;
	public:
		static constexpr unsigned int index=underlying_type::index;
		using type = typename underlying_type::type;
	};
	
	template<typename Name, typename... FieldTypes>
	struct IndexForName{
	private:
		using underlying_type = IndexForNameImpl<0,Name,FieldTypes...>;
	public:
		static constexpr unsigned int index=underlying_type::index;
		using type = typename underlying_type::type;
	};
	
}

struct RecordID{
	uint32_t run;
	uint32_t event;
	uint32_t subEvent;
    int subEventStream;
	bool exists;
	
	bool operator==(const RecordID& other) const{
		return(run==other.run && event==other.event && subEventStream==other.subEventStream && subEvent==other.subEvent);
	}
	
	bool operator<(const RecordID& other) const{
		if(run<other.run)
			return(true);
		if(run>other.run)
			return(false);
		if(event<other.event)
			return(true);
		if(event>other.event)
			return(false);
		if(subEventStream<other.subEventStream)
			return(true);
		if(subEventStream>other.subEventStream)
			return(false);
		if(subEvent<other.subEvent)
			return(true);
		return(false);
	}
};

template<class... FieldTypes>
struct TableRow{
	typedef detail::hoboTuple<FieldTypes...> FieldsType;
	
	RecordID id;
	FieldsType fields;
	
	static const int baseFields=5;
	static const int nFields=baseFields+FieldsType::nFields;
	
	TableRow() = default;
	
	template<typename... ArgTypes>
	TableRow(ArgTypes... args):fields(args...){}
	
	static void getFieldSizes(size_t sizes[]){
		sizes[0]=sizeof(uint32_t);
		sizes[1]=sizeof(uint32_t);
		sizes[2]=sizeof(uint32_t);
		sizes[3]=sizeof(int);
		sizes[4]=sizeof(bool);
		FieldsType::getFieldSizes(sizes,baseFields);
	}
	
	static void getFieldOffsets(size_t offsets[]){
		typedef TableRow<FieldTypes...> ThisType;
		offsets[0]=offsetof(ThisType,id.run);
		offsets[1]=offsetof(ThisType,id.event);
		offsets[2]=offsetof(ThisType,id.subEvent);
		offsets[3]=offsetof(ThisType,id.subEventStream);
		offsets[4]=offsetof(ThisType,id.exists);
		FieldsType::getFieldOffsets(offsets,baseFields,offsetof(ThisType,fields));
	}
	
	static void getFieldNames(std::string names[nFields]){
		names[0]="Run";
		names[1]="Event";
		names[2]="SubEvent";
		names[3]="SubEventStream";
		names[4]="exists";
		FieldsType::getFieldNames(names,baseFields);
	}

	static void getFieldTypes(hid_t types[nFields]){
		types[0]=H5T_NATIVE_UINT32;
		types[1]=H5T_NATIVE_UINT32;
		types[2]=H5T_NATIVE_UINT32;
		types[3]=H5T_NATIVE_INT32;
		types[4]=H5T_NATIVE_UCHAR;
		FieldsType::getFieldTypes(types,baseFields);
	}

	static bool getAutoConvert(size_t fieldIndex){
		if(fieldIndex<baseFields)
			return(false);
		return(FieldsType::getAutoConvert(fieldIndex-baseFields));
	}
	
	template<unsigned int index>
	typename detail::FieldType_<index,FieldTypes...>::type get() const{
		static_assert(index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<index,FieldTypes...>::get(fields));
	}
	
	template<unsigned int index>
	typename detail::FieldType_<index,FieldTypes...>::type& get(){
		static_assert(index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<index,FieldTypes...>::get(fields));
	}
	
	template<typename Name, typename Info=detail::IndexForName<Name,FieldTypes...>>
	const typename Info::type get() const{
		static_assert(Info::index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<Info::index,FieldTypes...>::get(fields));
	}
	
	template<typename Name, typename Info=detail::IndexForName<Name,FieldTypes...>>
	typename Info::type get(){
		static_assert(Info::index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<Info::index,FieldTypes...>::get(fields));
	}
};

struct H5File{
	hid_t file;
	H5File(std::string path, unsigned int flags = H5F_ACC_RDONLY, hid_t fapl_id = H5P_DEFAULT):
	file(H5Fopen(path.c_str(), flags, fapl_id)){}
	~H5File(){
		H5Fclose(file);
	}
	operator hid_t() const{ //implicit conversion
		return(file);
	}
	operator bool() const{ //conversion to bool, tests if file is valid
		return(file>=0);
	}
};

namespace detail{
	//this is supposed to live in readTableBlock, but ATM using a local struct for a deleter seems to make clang crash
	struct nameDeleter{
		unsigned int size;
		nameDeleter(unsigned int s):size(s){}
		void operator()(char** obj){
			for(unsigned int i=0; i<size; i++)
				delete[] (obj[i]);
			delete[] obj;
		}
	};

	//H5Tconvert has massive startup overhead, apparently due to very cautious sanity checking.
	//Unfortunately, while it has the ability to convert multiple items at a time, those items
	//must be located contiguously in memory. As a result, to convert one field of a number of
	//structs it's vastly more efficient to copy the scattered data aside into a contiguous
	//buffer, do the conversion, and copy it back.
	//The types being converted had beeter be POD for this to work safely, but that is almost 
	//certainly going to be the case.
	template<typename RecordType>
	void convertField(const hid_t sourceType, const hid_t targetType, const size_t fieldSize, const ptrdiff_t fieldOffset, std::vector<RecordType>& data){
		std::unique_ptr<char[]> buffer(new char[fieldSize*data.size()]); //TODO: worry about alignment
		char* buf_ptr=buffer.get();
		for(size_t i=0, n=data.size(); i<n; i++)
			memcpy(buf_ptr+i*fieldSize, (char*)&data[i]+fieldOffset, fieldSize);

		H5Tconvert(sourceType,targetType,data.size(),buf_ptr,NULL,H5P_DEFAULT);	

		for(size_t i=0, n=data.size(); i<n; i++)
			memcpy((char*)&data[i]+fieldOffset, buf_ptr+i*fieldSize, fieldSize);
	}
} //namespace detail

//read up to 1 MB of data at a time
template<typename T, typename Event, typename Callback>
void readTable(hid_t tableLoc, const std::string& tableName, std::map<RecordID,Event>& data, Callback callback){
	hsize_t nFields,nRecords;
	H5TBget_table_info(tableLoc,tableName.c_str(),&nFields,&nRecords);
	//TODO: check for errors
	assert(nFields>=T::nFields);
	std::unique_ptr<char*[],detail::nameDeleter> availableFieldNames(new char*[nFields],detail::nameDeleter(nFields));
	for(unsigned int i=0; i<nFields; i++)
		availableFieldNames[i]=new char[255];
	std::unique_ptr<size_t[]> availableFieldSizes(new size_t[nFields]);
	std::unique_ptr<size_t[]> availableFieldOffsets(new size_t[nFields]);
	size_t tableRowSize;
	H5TBget_field_info(tableLoc,tableName.c_str(),availableFieldNames.get(),availableFieldSizes.get(),availableFieldOffsets.get(),&tableRowSize);
	//TODO: check for errors
	
	std::unique_ptr<std::string[]> requestedFieldNames(new std::string[T::nFields]);
	T::getFieldNames(requestedFieldNames.get());
	
	int fieldIndices[T::nFields];
	size_t fieldSizes[T::nFields], fieldOffsets[T::nFields];
	hid_t expectedFieldTypes[T::nFields], actualFieldTypes[T::nFields];
	T::getFieldSizes(fieldSizes);
	T::getFieldOffsets(fieldOffsets);
	T::getFieldTypes(expectedFieldTypes);
	unsigned int lastIndex=0;

	for(unsigned int i=0; i<T::nFields; i++){
		for(; lastIndex<nFields; lastIndex++){
			if(availableFieldNames[lastIndex]==requestedFieldNames[i])
				break;
		}
		if(lastIndex==nFields)
			throw std::runtime_error("Field '"+requestedFieldNames[i]+
									 "' either does not exist in table '"+tableName+"' or was requested out of order");
		fieldIndices[i]=lastIndex;
		if(availableFieldSizes[lastIndex]!=fieldSizes[i])
			throw std::runtime_error("Field '"+requestedFieldNames[i]+
									 "' size does not match: target size is "+boost::lexical_cast<std::string>(fieldSizes[i])+
									 " but size in table is "+boost::lexical_cast<std::string>(availableFieldSizes[lastIndex]));
	}

	hid_t table=H5Dopen(tableLoc, tableName.c_str(),H5P_DEFAULT);
	hid_t tableType=H5Dget_type(table);
	for(unsigned int i=0; i<T::nFields; i++)
		actualFieldTypes[i]=H5Tget_member_type(tableType, fieldIndices[i]);
	H5Tclose(tableType);
	H5Dclose(table);
	struct typeListDeleter{
		hid_t* types;
		size_t count;
		typeListDeleter(hid_t* t, size_t c):types(t),count(c){}
		~typeListDeleter(){
			for(size_t i=0; i<count; i++)
				H5Tclose(types[i]);
		}
	} actualFieldTypesCleanup(actualFieldTypes,T::nFields);
	
	std::vector<unsigned int> indicesToConvert;
	for(unsigned int i=T::baseFields; i<T::nFields; i++){
		//if the on-disk and in-memory types don't match and the user has allowed auto-conversion
		//try to set it up
		if(H5Tequal(expectedFieldTypes[i],actualFieldTypes[i])==false && T::getAutoConvert(i)){
			indicesToConvert.push_back(i);
			if(H5Tcompiler_conv(actualFieldTypes[i],expectedFieldTypes[i])==-1){
				//no conversion available; try to figure out why
				std::ostringstream reason;
				if(H5Tget_class(actualFieldTypes[i])!=H5Tget_class(expectedFieldTypes[i]))
					reason << "Unable to convert datatypes for field " << i << " of table " << tableName << " which have different classes";
				else if(H5Tget_class(actualFieldTypes[i])==H5T_ENUM){
					//a common reason is that one enumerator in the source type is not in the destination type
					int n=H5Tget_nmembers(actualFieldTypes[i]);
					for(int j=0; j<n; j++){
						char* name=H5Tget_member_name(actualFieldTypes[i],j);
						int oldval, newidx, newval;
						H5Tget_member_value(actualFieldTypes[i],j,&oldval);
						newidx=H5Tget_member_index(expectedFieldTypes[i],name);
						if(newidx<0)
							reason << "Value '" << name << "' (=" << oldval << ") from source (on-disk) enum type not found in target (in-memory) type" << std::endl;
		
						free(name);
					}
				}
				if(!reason.str().empty())
					throw std::runtime_error(reason.str());

				throw std::runtime_error("Unable to convert datatypes for field "+std::to_string(i));
			}
		}
	}

	const unsigned int maxBlockSize=1<<20; //1 MB
	unsigned int recordsPerBlock=maxBlockSize/sizeof(T);
	if(!recordsPerBlock)
		recordsPerBlock=1;
	std::vector<T> intermediate(recordsPerBlock);
	
	hsize_t recordsRead=0;
	while(recordsRead<nRecords){
		hsize_t toRead=recordsPerBlock;
		if(recordsRead+toRead>nRecords)
			toRead=nRecords-recordsRead;
		intermediate.resize(toRead);
		T* dataPtr=&intermediate[0];
		if(H5TBread_fields_index(tableLoc, tableName.c_str(), T::nFields, fieldIndices, recordsRead, toRead, sizeof(T),  fieldOffsets, fieldSizes, dataPtr)<0)
			throw std::runtime_error("Read error");

		//do any necessary conversions
		for(auto index : indicesToConvert)
			detail::convertField<T>(actualFieldTypes[index],expectedFieldTypes[index],fieldSizes[index],fieldOffsets[index],intermediate);

		//send newly read data to the callback
		for(T& sourceItem : intermediate){
			Event& targetItem=data[sourceItem.id];
			callback(sourceItem,targetItem);
		}
		recordsRead+=toRead;
	}
}

} //namespace tableio
} //namespace phys_tools

#endif
