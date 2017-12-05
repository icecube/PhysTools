#include "PhysTools/axis.h"

#include <map>

namespace phys_tools{
	namespace histograms{

		namespace{
			std::map<std::string,axis* (*)(hid_t axis_id)> axisDeserializers;
		}
		
#define RegisterAxisType(axis)\
struct Register ## axis { \
	Register ## axis (){ \
		axisDeserializers.emplace(#axis,&axis::readHDF); \
	} \
} reg ## axis ## _s;

		void axis::writeBasicPropertiesHDF5(hid_t dataset_id) const{
			hdf_interface::addAttribute(dataset_id,"base",base);
			hdf_interface::addAttribute(dataset_id,"stride",stride);
			hdf_interface::addAttribute(dataset_id,"offset",offset);
			hdf_interface::addAttribute(dataset_id,"count",count);
			hdf_interface::addAttribute(dataset_id,"lowerLimit",lowerLimit);
			hdf_interface::addAttribute(dataset_id,"upperLimit",upperLimit);
		}
		
		void LinearAxis::writeHDF(hid_t container, unsigned int idx) const{
			uint8_t dummy=0;
			hsize_t dim=1;
			hid_t dType=hdf_interface::getHDFDatatype<uint8_t>().datatype;
			hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
			hid_t dataset_id = H5Dcreate(container, ("_h_axis_"+std::to_string(idx)).c_str(), dType, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			
			H5Dwrite(dataset_id, dType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummy);
			
			writeBasicPropertiesHDF5(dataset_id);
			hdf_interface::addAttribute(dataset_id,"AxisType",std::string("Linear"));
			
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		axis* LinearAxis::readHDF(hid_t axis_id){
			externalCoordinate base, stride, lowerLimit, upperLimit;
			internalCoordinate count, offset;
			hdf_interface::readAttribute(axis_id,"base",base);
			hdf_interface::readAttribute(axis_id,"stride",stride);
			hdf_interface::readAttribute(axis_id,"offset",offset);
			hdf_interface::readAttribute(axis_id,"count",count);
			hdf_interface::readAttribute(axis_id,"lowerLimit",lowerLimit);
			hdf_interface::readAttribute(axis_id,"upperLimit",upperLimit);
			H5Dclose(axis_id);
			return(new LinearAxis(base,stride,offset,count,lowerLimit,upperLimit));
		}
		
		RegisterAxisType(LinearAxis);
		
		void LogarithmicAxis::writeHDF(hid_t container, unsigned int idx) const{
			uint8_t dummy=0;
			hsize_t dim=1;
			hid_t dType=hdf_interface::getHDFDatatype<uint8_t>().datatype;
			hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
			hid_t dataset_id = H5Dcreate(container, ("_h_axis_"+std::to_string(idx)).c_str(), dType, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			
			H5Dwrite(dataset_id, dType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummy);
			hdf_interface::addAttribute(dataset_id,"AxisType",std::string("Logarithmic"));
			
			writeBasicPropertiesHDF5(dataset_id);
			
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		axis* LogarithmicAxis::readHDF(hid_t axis_id){
			externalCoordinate base, stride, lowerLimit, upperLimit;
			internalCoordinate count, offset;
			hdf_interface::readAttribute(axis_id,"base",base);
			hdf_interface::readAttribute(axis_id,"stride",stride);
			hdf_interface::readAttribute(axis_id,"offset",offset);
			hdf_interface::readAttribute(axis_id,"count",count);
			hdf_interface::readAttribute(axis_id,"lowerLimit",lowerLimit);
			hdf_interface::readAttribute(axis_id,"upperLimit",upperLimit);
			H5Dclose(axis_id);
			return(new LogarithmicAxis(base,stride,offset,count,lowerLimit,upperLimit));
		}
		
		RegisterAxisType(LogarithmicAxis);
		
		void FixedLinearAxis::writeHDF(hid_t container, unsigned int idx) const{
			uint8_t dummy=0;
			hsize_t dim=1;
			hid_t dType=hdf_interface::getHDFDatatype<uint8_t>().datatype;
			hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
			hid_t dataset_id = H5Dcreate(container, ("_h_axis_"+std::to_string(idx)).c_str(), dType, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			
			H5Dwrite(dataset_id, dType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummy);
			
			writeBasicPropertiesHDF5(dataset_id);
			hdf_interface::addAttribute(dataset_id,"max",max);
			hdf_interface::addAttribute(dataset_id,"totalCount",totalCount);
			hdf_interface::addAttribute(dataset_id,"AxisType",std::string("FixedLinear"));
			
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		axis* FixedLinearAxis::readHDF(hid_t axis_id){
			externalCoordinate base, stride, lowerLimit, upperLimit;
			internalCoordinate count, offset;
			hdf_interface::readAttribute(axis_id,"base",base);
			hdf_interface::readAttribute(axis_id,"stride",stride);
			hdf_interface::readAttribute(axis_id,"offset",offset);
			hdf_interface::readAttribute(axis_id,"count",count);
			hdf_interface::readAttribute(axis_id,"lowerLimit",lowerLimit);
			hdf_interface::readAttribute(axis_id,"upperLimit",upperLimit);
			H5Dclose(axis_id);
			return(new FixedLinearAxis(base,stride,offset,count,lowerLimit,upperLimit));
		}
		
		RegisterAxisType(FixedLinearAxis);
		
		void FixedLogarithmicAxis::writeHDF(hid_t container, unsigned int idx) const{
			uint8_t dummy=0;
			hsize_t dim=1;
			hid_t dType=hdf_interface::getHDFDatatype<uint8_t>().datatype;
			hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
			hid_t dataset_id = H5Dcreate(container, ("_h_axis_"+std::to_string(idx)).c_str(), dType, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			
			H5Dwrite(dataset_id, dType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummy);
			
			writeBasicPropertiesHDF5(dataset_id);
			hdf_interface::addAttribute(dataset_id,"max",max);
			hdf_interface::addAttribute(dataset_id,"totalCount",totalCount);
			hdf_interface::addAttribute(dataset_id,"AxisType",std::string("FixedLogarithmic"));
			
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		axis* FixedLogarithmicAxis::readHDF(hid_t axis_id){
			externalCoordinate base, stride, lowerLimit, upperLimit;
			internalCoordinate count, offset;
			hdf_interface::readAttribute(axis_id,"base",base);
			hdf_interface::readAttribute(axis_id,"stride",stride);
			hdf_interface::readAttribute(axis_id,"offset",offset);
			hdf_interface::readAttribute(axis_id,"count",count);
			hdf_interface::readAttribute(axis_id,"lowerLimit",lowerLimit);
			hdf_interface::readAttribute(axis_id,"upperLimit",upperLimit);
			H5Dclose(axis_id);
			return(new FixedLogarithmicAxis(base,stride,offset,count,lowerLimit,upperLimit));
		}
		
		RegisterAxisType(FixedLogarithmicAxis);
		
		void FixedUserDefinedAxis::writeHDF(hid_t container, unsigned int idx) const{
			hsize_t dim=edges.size();
			hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
			hid_t dType=hdf_interface::getHDFDatatype<externalCoordinate>().datatype;
			hid_t dataset_id = H5Dcreate(container, ("_h_binedges_"+std::to_string(idx)).c_str(), dType, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			
			H5Dwrite(dataset_id, dType, H5S_ALL, H5S_ALL, H5P_DEFAULT, edges.data());
			hdf_interface::addAttribute(dataset_id,"AxisType",std::string("FixedUserDefined"));
			hdf_interface::addAttribute(dataset_id,"offset",offset);
			hdf_interface::addAttribute(dataset_id,"count",count);
			
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		axis* FixedUserDefinedAxis::readHDF(hid_t axis_id){
			hid_t dataspace_id = H5Dget_space(axis_id);
			int axis_dim = H5Sget_simple_extent_ndims(dataspace_id);
			if(axis_dim!=1){
				H5Dclose(axis_id);
				H5Sclose(dataspace_id);
				throw std::runtime_error("Rank of axis data is not 1");
			}
			hsize_t extent;
			H5Sget_simple_extent_dims(dataspace_id,&extent,NULL);
			H5Sclose(dataspace_id);
			if(extent<2){
				H5Dclose(axis_id);
				throw std::runtime_error("Number of bin edges for axis is too small");
			}
			
			std::vector<double> edges(extent);
			herr_t status = H5Dread(axis_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &edges[0]);
			if(status<0){
				H5Dclose(axis_id);
				throw std::runtime_error("Failed to read edges for axis");
			}
			
			internalCoordinate count, offset;
			hdf_interface::readAttribute(axis_id,"offset",offset);
			hdf_interface::readAttribute(axis_id,"count",count);
			
			H5Dclose(axis_id);
			auto result=new FixedUserDefinedAxis(edges.begin(),edges.end());
			result->offset=offset;
			result->count=count;
			return(result);
		}
		
		RegisterAxisType(FixedUserDefinedAxis);
		
		axis* deserializeAxis(hid_t container, const std::string& name){
			hid_t axis_id = H5Dopen(container, name.c_str(), H5P_DEFAULT);
			if(axis_id<0)
				throw std::runtime_error("Failed to open axis");
			
			std::string axisType;
			try{
				hdf_interface::readAttribute(axis_id,"AxisType",axisType);
				axisType+="Axis";
			}catch(...){
				H5Dclose(axis_id);
				throw;
			}
			
			auto desIt=axisDeserializers.find(axisType);
			if(desIt==axisDeserializers.end()){
				throw std::runtime_error("Unrecognized axis type: '"+axisType+"'");
				H5Dclose(axis_id);
			}
			
			return((*desIt->second)(axis_id));
		}
		
	} //namespace histograms
} //namespace phys_tools