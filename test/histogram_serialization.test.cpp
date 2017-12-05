#include <iostream>
#include <PhysTools/hdf5_serialization.h>
#include <PhysTools/histogram.h>
#include <PhysTools/bin_types.h>

//#include <cxxabi.h>

/*struct thingy{
	double x;
	double y;
	struct{
		int a, b;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::make_nvp("a",a);
			ar & boost::serialization::make_nvp("b",b);
		}
	} z;
	
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        //ar & boost::serialization::make_nvp("x",x);
		ar & x;
        ar & boost::serialization::make_nvp("y",y);
        ar & boost::serialization::make_nvp("z",z);
    }
};*/



int main(){
	/*std::cout << "hash_code for char: " << typeid(std::declval<char>()).hash_code() << std::endl;
	
	for(const auto& item : phys_tools::hdf_interface::DatatypeRegistry){
		int status = 0;
		char* demangled = abi::__cxa_demangle(item.first.name(), 0, 0, &status);
		std::cout << "Datatype for type " << (status?item.first.name():demangled) << " is " << item.second.datatype << std::endl;
		if(demangled)
			free(demangled);
	}
	
	std::cout << phys_tools::hdf_interface::getHDFDatatype<char>().datatype << std::endl;
	hid_t thingy_t=phys_tools::hdf_interface::getHDFDatatype<thingy>().datatype;
	std::cout << thingy_t << std::endl;
	std::cout << H5Tget_nmembers(thingy_t) << std::endl;
	for(unsigned int i=0; i<H5Tget_nmembers(thingy_t); i++){
		hid_t mType=H5Tget_member_type(thingy_t,i);
		H5T_class_t mClass=H5Tget_member_class(thingy_t,i);
		std::cout << " field " << i << ": " << H5Tget_member_name(thingy_t,i)
		<< " of type " << mType << std::endl;
		H5Tclose(mType);
	}*/
	
	using namespace phys_tools::histograms;
	{
		histogram<2> h1(LinearAxis(0,1),LinearAxis(0,1));
		h1.add(0,0,amount(1));
		h1.add(0,1,amount(2));
		h1.add(1,0,amount(3));
		h1.add(1,1,amount(4));
		hid_t h5file=H5Fcreate("histogram_serialization_test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		//standard mode
		h1.write(h5file, "h1", /*compressLevel=*/0, /*dashiCompat=*/false);
		//compressed
		h1.write(h5file, "h1_compressed", /*compressLevel=*/1, /*dashiCompat=*/false);
		//Dashi-compatible
		h1.write(h5file, "h1_dashi", /*compressLevel=*/0, /*dashiCompat=*/true);
		H5Fclose(h5file);
	}
	{
		histogram<2> h1;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h1.read(h5file, "h1");
		H5Fclose(h5file);
		std::cout << h1 << std::endl;
	}
	{
		histogram<2> h1;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h1.read(h5file, "h1_compressed");
		H5Fclose(h5file);
		std::cout << h1 << std::endl;
	}
	{
		histogram<2> h1;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h1.read(h5file, "h1_dashi");
		H5Fclose(h5file);
		std::cout << h1 << std::endl;
	}
	
	{
		histogram<Dynamic> h2(LinearAxis(0,1),LinearAxis(0,1));
		h2.add(0,0,amount(11));
		h2.add(0,1,amount(12));
		h2.add(1,0,amount(13));
		h2.add(1,1,amount(14));
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDWR, H5P_DEFAULT);
		h2.write(h5file, "h2", /*compressLevel=*/0, /*dashiCompat=*/false);
		H5Fclose(h5file);
	}
	{
		histogram<Dynamic> h2;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h2.read(h5file, "h2");
		H5Fclose(h5file);
		std::cout << h2 << std::endl;
	}
	
	{
		histogram<1,entryCounter<double>> h3(LinearAxis(0,1));
		h3.add(0,amount(1));
		h3.add(0,amount(2));
		h3.add(0,amount(3));
		h3.add(1,amount(4));
		h3.add(1,amount(5));
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDWR, H5P_DEFAULT);
		h3.write(h5file, "h3", /*compressLevel=*/0, /*dashiCompat=*/false);
		H5Fclose(h5file);
	}
	{
		histogram<1,entryCounter<double>> h3;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h3.read(h5file, "h3");
		H5Fclose(h5file);
		std::cout << h3 << std::endl;
	}
	
	{
		histogram<1,meanVarTracker<double>> h4(LinearAxis(0,1));
		h4.add(0,amount(1));
		h4.add(0,amount(2));
		h4.add(0,amount(3));
		h4.add(1,amount(4));
		h4.add(1,amount(5));
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDWR, H5P_DEFAULT);
		h4.write(h5file, "h4", /*compressLevel=*/0, /*dashiCompat=*/false);
		H5Fclose(h5file);
	}
	{
		histogram<1,meanVarTracker<double>> h4;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h4.read(h5file, "h4");
		H5Fclose(h5file);
		std::cout << h4 << std::endl;
	}
	
	{
		histogram<1,sqErrorValue<double>> h5(LinearAxis(0,1));
		h5.add(0,amount(1));
		h5.add(0,amount(2));
		h5.add(0,amount(3));
		h5.add(1,amount(4));
		h5.add(1,amount(5));
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDWR, H5P_DEFAULT);
		h5.write(h5file, "h5", /*compressLevel=*/0, /*dashiCompat=*/false);
		H5Fclose(h5file);
	}
	{
		histogram<1,sqErrorValue<double>> h5;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h5.read(h5file, "h5");
		H5Fclose(h5file);
		std::cout << h5 << std::endl;
	}
	
	{
		histogram<1,fcErrorValue> h6(LinearAxis(0,1));
		h6.add(0,amount(1));
		h6.add(0,amount(2));
		h6.add(0,amount(3));
		h6.add(1,amount(4));
		h6.add(1,amount(5));
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDWR, H5P_DEFAULT);
		h6.write(h5file, "h6", /*compressLevel=*/0, /*dashiCompat=*/false);
		H5Fclose(h5file);
	}
	{
		histogram<1,fcErrorValue> h6;
		hid_t h5file=H5Fopen("histogram_serialization_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		h6.read(h5file, "h6");
		H5Fclose(h5file);
		std::cout << h6 << std::endl;
	}
}