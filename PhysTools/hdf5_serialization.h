///\file histogram_serialization.h
///This file contains support code for (de)serialization of histograms to and from HDF5

#ifndef PHYSTOOLS_HDF5_SERIALIZATION_H
#define PHYSTOOLS_HDF5_SERIALIZATION_H

#include <cassert>
#include <numeric>
#include <set>
#include <type_traits>
#include <typeindex>
#include <unordered_map>
#include <vector>

#include <boost/serialization/nvp.hpp>
#include <boost/lexical_cast.hpp>

#include "hdf5.h"

#include <iostream>

namespace phys_tools{	
	namespace hdf_interface{
		
		///\brief A handle for an HDF5 datatype
		class HDFDatatype{
		public:
			const hid_t datatype;
			
			~HDFDatatype();
			//copies are not okay
			HDFDatatype(const HDFDatatype&)=delete;
			//moves are okay
			HDFDatatype(HDFDatatype&& other):
			datatype(other.datatype),
			disposable(other.disposable){
				other.disposable=false;
			}
			//no assigment
			HDFDatatype& operator=(const HDFDatatype&)=delete;
			HDFDatatype& operator=(HDFDatatype&&)=delete;
		private:
			HDFDatatype(hid_t t, bool d):datatype(t),disposable(d){}
			bool disposable;
			template<typename T>
			friend void registerHDFDatatype(hid_t,bool);
		};
		
		///Tracks all HDF datatypes which are in use
		extern std::unordered_map<std::type_index,HDFDatatype> DatatypeRegistry;
		
		///Registers an HDF datatype corresponding to a C++ type
		template<typename T>
		void registerHDFDatatype(hid_t t, bool d){
			HDFDatatype hd(t,d);
			auto result=DatatypeRegistry.emplace(typeid(T),std::move(hd));
			//result.first->second.disposable=d;
		}
		
		//fwd decl
		template<typename T>
		const HDFDatatype& synthesizeHDFDatatype();
		
		///Checks whether there is a registered HDF datatype for the given C++ type
		template<typename T>
		bool HDFDatatypeExists(){
			auto it=DatatypeRegistry.find(typeid(T));
			return(it!=DatatypeRegistry.end());
		}
		
		///Fetches the HDF datatype for the given C++ datatype, creating and registering it if necessary
		template<typename T>
		const HDFDatatype& getHDFDatatype(){
			auto it=DatatypeRegistry.find(typeid(T));
			if(it==DatatypeRegistry.end())
				return(synthesizeHDFDatatype<T>());
			return(it->second);
		}
		
		namespace detail{
			///\brief Creates an HDF datatype from a C++ type.
			///
			///This is accomplished by serializing an instance of the type through the boost::serialization
			///interface to discover its structure.
			///TODO: make sure that this works for all the alternate ways boost::serialization allows
			///serialization to be set up
			class HDFDatatypeSynthesizer{
				hid_t datatype;
				void* objBase;
				unsigned int anonymousCounter;
				std::vector<std::string> subobjectNames;
				std::string currentItemName() const{
					assert(!subobjectNames.empty());
					return(std::accumulate(subobjectNames.begin()+1,subobjectNames.end(),subobjectNames.front(),
										   [](const std::string& s1, const std::string& s2){ return(s1+"."+s2); }));
				}
			public:
				HDFDatatypeSynthesizer():anonymousCounter(0){}
				
				//Behave as a Saving Archive
				using is_saving=boost::mpl::bool_<true>;
				using is_loading=boost::mpl::bool_<true>;
				
				template<typename T>
				void analyze(T&& t){
					datatype = H5Tcreate(H5T_COMPOUND, sizeof(T));
					//std::cout << "Created datatype " << datatype << std::endl;
					objBase=&t;
					visitor<T>::visit(*this,t);
				}
				
				hid_t getDatatype() const{ return(datatype); }
				
				template<typename T, bool S=std::is_scalar<T>::value >
				struct visitor{
					static void visit(HDFDatatypeSynthesizer& a, T& t){
						t.serialize(a,0);
					}
				};
				template<typename T>
				struct visitor<T,true>{
					static void visit(HDFDatatypeSynthesizer& a, T& t){
						throw std::logic_error("visit should never be called for scalars");
					}
				};
				//TODO: array members?
				
				template<typename T>
				HDFDatatypeSynthesizer& operator&(T& t){
					subobjectNames.push_back("anon_field_"+boost::lexical_cast<std::string>(anonymousCounter++));
					//std::cout << " inserting field " << currentItemName() << " at offset " << ((char*)&t-(char*)objBase) << std::endl;
					if(!HDFDatatypeExists<T>())
						synthesizeHDFDatatype<T>();
					H5Tinsert(datatype, currentItemName().c_str(),
							  ((char*)&t-(char*)objBase),
							  getHDFDatatype<T>().datatype);
					subobjectNames.pop_back();
					return(*this);
				}
				
				template<typename T>
				HDFDatatypeSynthesizer& operator&(const boost::serialization::nvp<T>& t){
					subobjectNames.push_back(t.name());
					if(!HDFDatatypeExists<T>())
						synthesizeHDFDatatype<T>();
					H5Tinsert(datatype, currentItemName().c_str(),
							  ((char*)&t.value()-(char*)objBase),
							  getHDFDatatype<T>().datatype);
					subobjectNames.pop_back();
					return(*this);
				}
				
				template<typename T>
				HDFDatatypeSynthesizer& operator<<(const T& t){
					return(*this & t);
				}
				//Don't need operator>> since we never load
			};
		}
		
		template<typename T>
		const HDFDatatype& synthesizeHDFDatatype(){
			detail::HDFDatatypeSynthesizer s;
			s.analyze(T());
			//std::cout << "Allocated HDF datatype " << s.getDatatype() << " for C++ type " << typeid(T).name() << std::endl;
			registerHDFDatatype<T>(s.getDatatype(),true);
			return(getHDFDatatype<T>());
		}
			
		template<typename T>
		void addAttribute(hid_t object, std::string name, const T& contents){
			hid_t dtype = getHDFDatatype<T>().datatype;
			hsize_t dim=1;
			hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
			hid_t attribute_id = H5Acreate(object,name.c_str(),dtype,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
			H5Awrite(attribute_id, dtype, &contents);
			H5Aclose(attribute_id);
			H5Sclose(dataspace_id);
		};
		
		template<>
		void addAttribute<std::string>(hid_t object, std::string name, const std::string& contents);
		
		template<typename T>
		void readAttribute(hid_t object, std::string name, T& dest){
			hid_t attribute_id = H5Aopen(object, name.c_str(), H5P_DEFAULT);
			hid_t actualType = H5Aget_type(attribute_id);
			
			hid_t expectedType = getHDFDatatype<T>().datatype;
			
			//TODO: this may be too harsh; it ignores whether conversion is possible
			if(H5Tequal(actualType, expectedType)<=0){
				H5Aclose(attribute_id);
				H5Tclose(actualType);
				throw std::runtime_error("Expected and actual data types for attribute '"+name+"' do not match");
			}
			
			if(H5Aread(attribute_id, expectedType, &dest)<0){
				H5Aclose(attribute_id);
				H5Tclose(actualType);
				throw std::runtime_error("Failed to read attribute '"+name+"'");
			}
			
			H5Aclose(attribute_id);
			H5Tclose(actualType);
		}
		
		template<>
		void readAttribute<std::string>(hid_t object, std::string name, std::string& dest);
		
		std::set<std::string> groupContents(hid_t group_id);
		
	} //namespace hdf_interface
} //namespace phys_tools

#endif //PHYSTOOLS_HDF5_SERIALIZATION_H
