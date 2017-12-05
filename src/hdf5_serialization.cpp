#include "PhysTools/hdf5_serialization.h"

namespace phys_tools{
	namespace hdf_interface{
		
		HDFDatatype::~HDFDatatype(){
			if(disposable){
				//The HDF5 documentation says we must call H5Tclose on every type we allocate,
				//however, doing so yields an error:
				//	H5T.c line 1722 in H5Tclose(): not a datatype
				//at runtime with version 1.8.9 of the library. So, for now, we leak.
				//std::cout << "Deallocating HDF datatype " << datatype << std::endl;
				//H5Tclose(datatype);
			}
		}
		
		std::unordered_map<std::type_index,HDFDatatype> DatatypeRegistry;

		namespace detail{
			///Registers all standard atomic datatypes which we are likely to want
			class basic_serialization_traits_init{
			public:
				basic_serialization_traits_init(){
					registerHDFDatatype<char>(H5T_NATIVE_CHAR,false);
					registerHDFDatatype<signed char>(H5T_NATIVE_SCHAR,false);
					registerHDFDatatype<unsigned char>(H5T_NATIVE_UCHAR,false);
					
					registerHDFDatatype<short>(H5T_NATIVE_SHORT,false);
					registerHDFDatatype<unsigned short>(H5T_NATIVE_USHORT,false);
					
					registerHDFDatatype<int>(H5T_NATIVE_INT,false);
					registerHDFDatatype<unsigned int>(H5T_NATIVE_UINT,false);
					
					registerHDFDatatype<long>(H5T_NATIVE_LONG,false);
					registerHDFDatatype<unsigned long>(H5T_NATIVE_ULONG,false);
					
					registerHDFDatatype<long long>(H5T_NATIVE_LLONG,false);
					registerHDFDatatype<unsigned long long>(H5T_NATIVE_ULLONG,false);
					
					registerHDFDatatype<float>(H5T_NATIVE_FLOAT,false);
					registerHDFDatatype<double>(H5T_NATIVE_DOUBLE,false);
					if(H5_SIZEOF_LONG_DOUBLE!=0)
						registerHDFDatatype<long double>(H5T_NATIVE_LDOUBLE,false);
					
					registerHDFDatatype<int8_t>(H5T_NATIVE_INT8,false);
					registerHDFDatatype<uint8_t>(H5T_NATIVE_UINT8,false);
					registerHDFDatatype<int_least8_t>(H5T_NATIVE_INT_LEAST8,false);
					registerHDFDatatype<uint_least8_t>(H5T_NATIVE_UINT_LEAST8,false);
					registerHDFDatatype<int_fast8_t>(H5T_NATIVE_INT_FAST8,false);
					registerHDFDatatype<uint_fast8_t>(H5T_NATIVE_UINT_FAST8,false);
					
					registerHDFDatatype<int16_t>(H5T_NATIVE_INT16,false);
					registerHDFDatatype<uint16_t>(H5T_NATIVE_UINT16,false);
					registerHDFDatatype<int_least16_t>(H5T_NATIVE_INT_LEAST16,false);
					registerHDFDatatype<uint_least16_t>(H5T_NATIVE_UINT_LEAST16,false);
					registerHDFDatatype<int_fast16_t>(H5T_NATIVE_INT_FAST16,false);
					registerHDFDatatype<uint_fast16_t>(H5T_NATIVE_UINT_FAST16,false);
					
					registerHDFDatatype<int32_t>(H5T_NATIVE_INT32,false);
					registerHDFDatatype<uint32_t>(H5T_NATIVE_UINT32,false);
					registerHDFDatatype<int_least32_t>(H5T_NATIVE_INT_LEAST32,false);
					registerHDFDatatype<uint_least32_t>(H5T_NATIVE_UINT_LEAST32,false);
					registerHDFDatatype<int_fast32_t>(H5T_NATIVE_INT_FAST32,false);
					registerHDFDatatype<uint_fast32_t>(H5T_NATIVE_UINT_FAST32,false);
					
					registerHDFDatatype<int64_t>(H5T_NATIVE_INT64,false);
					registerHDFDatatype<uint64_t>(H5T_NATIVE_UINT64,false);
					registerHDFDatatype<int_least64_t>(H5T_NATIVE_INT_LEAST64,false);
					registerHDFDatatype<uint_least64_t>(H5T_NATIVE_UINT_LEAST64,false);
					registerHDFDatatype<int_fast64_t>(H5T_NATIVE_INT_FAST64,false);
					registerHDFDatatype<uint_fast64_t>(H5T_NATIVE_UINT_FAST64,false);
				}
			} initializeBasicDataTypes;
		}
			
		template<>
		void addAttribute<std::string>(hid_t object, std::string name, const std::string& contents){
			hid_t strtype = H5Tcopy(H5T_C_S1);
			H5Tset_size(strtype, contents.size());
			hsize_t dim=1;
			hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
			hid_t attribute_id = H5Acreate(object,name.c_str(),strtype,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
			H5Awrite(attribute_id, strtype, &contents[0]);
			H5Aclose(attribute_id);
			H5Sclose(dataspace_id);
		};
		
		template<>
		void readAttribute<std::string>(hid_t object, std::string name, std::string& dest){
			hid_t attribute_id = H5Aopen(object, name.c_str(), H5P_DEFAULT);
			hid_t actualType = H5Aget_type(attribute_id);
			H5T_class_t typeClass = H5Tget_class(actualType);
			if(typeClass!=H5T_STRING){
				H5Aclose(attribute_id);
				H5Tclose(actualType);
				throw std::runtime_error("Expected and actual data types for attribute '"+name+"' do not match");
			}
			
			size_t size = H5Tget_size(actualType);
			dest.resize(size);
			
			if(H5Aread(attribute_id, actualType, &dest[0])<0){
				H5Aclose(attribute_id);
				H5Tclose(actualType);
				throw std::runtime_error("Failed to read attribute '"+name+"'");
			}
			
			H5Aclose(attribute_id);
			H5Tclose(actualType);
		}
		
		herr_t collectGroupContents(hid_t group_id, const char* member_name, const H5L_info_t *info, void* operator_data){
			std::set<std::string>* items=static_cast<std::set<std::string>*>(operator_data);
			items->insert(member_name);
			return(0);
		}
		
		std::set<std::string> groupContents(hid_t group_id){
			hsize_t start=0;
			std::set<std::string> tables;
			H5Literate(group_id,H5_INDEX_NAME,H5_ITER_NATIVE,&start,&collectGroupContents,&tables);
			return(tables);
		}
	
	} //namespace hdf_interface
} //namespace phys_tools
