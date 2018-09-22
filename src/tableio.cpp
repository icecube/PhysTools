#include <PhysTools/tableio.h>

namespace phys_tools{
namespace tableio{

namespace{
herr_t collectTableNames(hid_t group_id, const char* member_name, void* operator_data){
	std::set<std::string>* items=static_cast<std::set<std::string>*>(operator_data);
	items->insert(member_name);
	return(0);
}
}

std::set<std::string> getTables(H5File& file, std::string loc){
	std::set<std::string> tables;
	H5Giterate(file,"/",NULL,&collectTableNames,&tables);
	return(tables);
}

} //namespace tableio
} //namespace phys_tools
