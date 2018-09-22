#ifndef PHYSTOOLS_MEMORYDUMP_H
#define PHYSTOOLS_MEMORYDUMP_H

#include <deque>
#include <istream>
#include <ostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <boost/crc.hpp>

namespace phys_tools{

///Compute a CRC32 checksum of the contents of the given file.
///Uses 1 MB of stack space as a buffer.
uint32_t getFileChecksum(const std::string& filename);

///Store a collection of events in raw (in memory) form to a file.
///Produces a file which is non-portable, but can be read very quickly by unsplatData.
///In case of changes to the Event structure, old files should not be read by
///different program versions. A checksum of the executable is stored in the
///file to prevent this.
///\param output the stream to which the data should be written
///\param progChecksum the checksum associated with the currently running code,
///                    which should be matched by any code attempting to load
///                    this data
///\param data the objects to be written
template<typename Event>
void dumpMemory(std::ostream& output, const uint32_t progChecksum,
			   const std::deque<Event>& data){
	size_t size;
	boost::crc_32_type fileCRC;

	output.write((char*)&progChecksum,sizeof(progChecksum));
	fileCRC.process_bytes((char*)&progChecksum,sizeof(progChecksum));

	size=data.size();
	output.write((char*)&size,sizeof(size));
	fileCRC.process_bytes((char*)&size,sizeof(size));
	for(const Event& e : data){
		output.write((char*)&e,sizeof(e));
		fileCRC.process_bytes((char*)&e,sizeof(e));
	}

	uint32_t checksum=fileCRC.checksum();
	output.write((char*)&checksum,sizeof(checksum));
}

///Read back a collection of data written by dumpMemory.
///The current execautable checksum is compared to the one stored in the file
///for safety.
///\param input the stream from which to read data
///\param expectedChecksum the checksum associated with the currently running
///                        code, which is checked against the checksum stored in
///                        the file
///\param data the container into which to load the data
///\param dangerous skip comparing checksums
template<typename Event>
void reloadMemory(std::istream& input, const uint32_t expectedChecksum,
				 std::deque<Event>& data, bool dangerous=false){
	boost::crc_32_type fileCRC;

	//read checksum; check it
	uint32_t storedChecksum;
	input.read((char*)&storedChecksum,sizeof(storedChecksum));
	fileCRC.process_bytes((char*)&storedChecksum,sizeof(storedChecksum));
	if(storedChecksum!=expectedChecksum){
		std::ostringstream ss;
		ss << "Stored program checksum, " << std::hex
		<< storedChecksum << ", does not match expected (current) checksum, " << expectedChecksum;
		throw std::runtime_error(ss.str());
	}

	size_t size;
	input.read((char*)&size,sizeof(size));
	fileCRC.process_bytes((char*)&size,sizeof(size));
	//data.resize(size);
	//for(Event& e : data){
	//	datafile.read((char*)&e,sizeof(e));
	//	fileCRC.process_bytes((char*)&e,sizeof(e));
	//}
	for(std::size_t i=0; i<size; i++){
		Event e;
		input.read((char*)&e,sizeof(e));
		fileCRC.process_bytes((char*)&e,sizeof(e));
		data.push_back(e);
	}

	input.read((char*)&storedChecksum,sizeof(storedChecksum));
	if(storedChecksum!=fileCRC.checksum()){
		std::ostringstream ss;
		ss << "Data appears to be corrupted: stored checksum, " << std::hex
		<< storedChecksum << ", does not match recomputed checksum, " << fileCRC.checksum();
		throw std::runtime_error(ss.str());
	}
}

} //namespace phys_tools

#endif
