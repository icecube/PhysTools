#include <PhysTools/memory_dump.h>

#include <fstream>

namespace phys_tools{

uint32_t getFileChecksum(const std::string& filename){
	const std::streamsize bufferSize=1u<<16; //1MB

	std::ifstream file(filename, std::ios_base::binary);
	if(!file)
		throw std::runtime_error("Unable to open "+filename+" for reading");

	boost::crc_32_type fileCRC;
	char buffer[bufferSize];
	do{
		file.read(buffer, bufferSize);
		fileCRC.process_bytes(buffer, file.gcount());
	} while(file);

	return(fileCRC.checksum());
}

} //namespace
