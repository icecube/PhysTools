// gnuplot-normalize.cpp
// Access Gnuplot directly and easily from C++ code
//
// This file implements the normalization methods

#include "PhysTools/gnuplot.h"

#include <fstream>
#include <sstream>

namespace{
	struct range_result{
		bool valid;
		double min, max;
		range_result():valid(false){}
	};

	range_result parseRange(const std::string& raw){
		range_result result;
		std::istringstream ss(raw);
		char separator;
		ss >> result.min >> separator >> result.max;
		if(!ss.fail() && !ss.bad())
			result.valid=true;
		return(result);
	}

	range_result extractRange(const std::string& raw){
		std::size_t start_pos=0, end_pos;
		range_result result;
		while(true){
			start_pos=raw.find('[',start_pos);
			if(start_pos==std::string::npos)
				break;
			end_pos=raw.find(']',start_pos);
			if(end_pos==std::string::npos)
				break;
			result=parseRange(raw.substr(start_pos+1,end_pos-start_pos-1));
			if(result.valid)
				break;
			start_pos++;
		}
		return(result);
	}
}

namespace phys_tools {
namespace gnuplot {

Gnuplot::NormalizeResult
Gnuplot::normalize_once (const std::string& gnuplot_code) const
{
    // send Gnuplot the code
    commands << "set term dumb" << std::endl;
    commands << "set output '/dev/null'" << std::endl;
    commands << gnuplot_code;
    // obtain x min and max
    NormalizeResult normalization;
    commands << "set xrange restore" << std::endl;
    commands << "show xrange" << std::endl;
	range_result result;
	std::string          line;
    while (getline (errors, line).good ()) {
		if((result=extractRange(line)).valid){
            normalization.x_min=result.min;
			normalization.x_max=result.max;
			break;
		}
    }
    // obtain y min and max
    commands << "set yrange restore" << std::endl;
    commands << "show yrange" << std::endl;
    while (getline (errors, line).good ()) {
		if((result=extractRange(line)).valid){
            normalization.y_min=result.min;
			normalization.y_max=result.max;
			break;
		}
    }
    return normalization;
}

} // namespace gnuplot
} // namespace phys_tools
