#include "PhysTools/plottable_histograms.h"

namespace phys_tools{

	std::string Plottable1DHistogram::get_data_code() const{
		unsigned int count=getBinCount(0);
		std::ostringstream data;
		if(format().style()==phys_tools::gnuplot::Format::LEFT_STEPS
		   || format().style()==phys_tools::gnuplot::Format::RIGHT_STEPS)
			data << getBinEdge(0,0) << '\t' << 0.0 << '\n';
		for(unsigned int i=0; i<count; i++){
			double x=0.0;
			switch(format().style()){
				case phys_tools::gnuplot::Format::LEFT_STEPS:
					x=getBinEdge(0,i);
					break;
				case phys_tools::gnuplot::Format::RIGHT_STEPS:
					x=getBinEdge(0,i)+getBinWidth(0,i);
					break;
				case phys_tools::gnuplot::Format::CENTER_STEPS:
				case phys_tools::gnuplot::Format::LINES:
				case phys_tools::gnuplot::Format::POINTS:
				case phys_tools::gnuplot::Format::LINES_POINTS:
				case phys_tools::gnuplot::Format::IMPULSES:
					x=getBinCenter(0,i);
					break;
				default:
					throw std::runtime_error("Unsupported drawing format");
			}
			data << x << '\t' << (*this)(i) << '\n';
		}
		if(format().style()==phys_tools::gnuplot::Format::LEFT_STEPS
		   || format().style()==phys_tools::gnuplot::Format::RIGHT_STEPS)
			data << getBinEdge(0,count) << '\t' << 0.0 << '\n';
		data << 'e' << std::endl;
		return(data.str());
	}
	
	std::string Plottable1DHistogram::get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const{
		return(get_data_code());
	}
	
	std::string Plottable1DHistogram::get_plot_code () const{
        std::ostringstream plot_code;
		plot_code << "'-' title \"" << title_ << "\" ";
        plot_code << format_;
        return (plot_code.str());
    }

	std::string Plottable2DHistogram::get_data_code() const{
		std::ostringstream data;
		unsigned int xCount=getBinCount(0);
		unsigned int yCount=getBinCount(1);
		for(unsigned int j=0; j<yCount; j++){
			if(format().style()==gnuplot::Format::IMAGE || format().style()==gnuplot::Format::IMAGE_FAILSAFE){
				for(unsigned int i=0; i<xCount; i++)
					data << getBinCenter(0,i) << '\t' << getBinCenter(1,j) << '\t' << (*this)(i,j) << '\n';
			}
			else{
				for(unsigned int i=0; i<xCount; i++)
					data << getBinEdge(0,i) << '\t' << getBinEdge(1,j) << '\t' << (*this)(i,j) << '\n';
				data << getBinEdge(0,xCount) << '\t' << getBinEdge(1,j) << '\t' << (*this)(xCount-1,j) << "\n\n";
			}
		}
		if(yCount && format().style()!=gnuplot::Format::IMAGE && format().style()!=gnuplot::Format::IMAGE_FAILSAFE){
			for(unsigned int i=0; i<xCount; i++)
				data << getBinEdge(0,i) << '\t' << getBinEdge(1,yCount) << '\t' << (*this)(i,yCount-1) << '\n';
			if(xCount)
				data << getBinEdge(0,xCount) << '\t' << getBinEdge(1,yCount) << '\t' << (*this)(xCount-1,yCount-1) << "\n";
		}
		data << "e\n";
		return(data.str());
	}
	
	std::string Plottable2DHistogram::get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const{
		return(get_data_code());
	}
	
	std::string Plottable2DHistogram::get_plot_code () const{
        std::ostringstream plot_code;
        plot_code << "'-' title \"" << title_ << "\" ";
        plot_code << format_;
        return (plot_code.str());
    }

	std::string PlottableErrorBars::get_data_code() const{
		unsigned int count=values.getBinCount(0);
		std::ostringstream data;
		for(unsigned int i=0; i<count; i++){
			double x=values.getBinEdge(0,i)+values.getBinWidth(0,i)/2.;
			if(horizontal_offset)
				x+=horizontal_offset*values.getBinWidth(0,i);
			data << x << '\t' << values(i) << '\t' << errors(i) << '\n';
		}
		data << 'e' << std::endl;
		return(data.str());
	}
	
	std::string PlottableErrorBars::get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const{
		return(get_data_code());
	}
	
	std::string PlottableErrorBars::get_plot_code () const{
        std::ostringstream plot_code;
        plot_code << "'-' notitle ";
        plot_code << format_;
		plot_code << "pointtype -1 ";
        return (plot_code.str());
    }
	
	/*std::string Plottable1DErrorHistogram::get_data_code() const{
		unsigned int count=getBinCount(0);
		std::ostringstream data;
		if(format().style()==phys_tools::gnuplot::Format::LEFT_STEPS
		   || format().style()==phys_tools::gnuplot::Format::RIGHT_STEPS)
			data << getBinEdge(0,0) << '\t' << 0.0 << '\n';
		for(unsigned int i=0; i<count; i++){
			double x=0.0;
			switch(format().style()){
				case phys_tools::gnuplot::Format::LEFT_STEPS:
					x=getBinEdge(0,i);
					break;
				case phys_tools::gnuplot::Format::RIGHT_STEPS:
					x=getBinEdge(0,i)+getBinWidth(0,i);
					break;
				case phys_tools::gnuplot::Format::CENTER_STEPS:
				case phys_tools::gnuplot::Format::Y_ERROR_BARS:
				case phys_tools::gnuplot::Format::Y_ERROR_LINES:
					x=getBinEdge(0,i)+getBinWidth(0,i)/2.;
					break;
				default:
					throw std::runtime_error("Unsupported drawing format");
			}
			if(horizontal_offset)
				x+=horizontal_offset*getBinWidth(0,i);
			data << x << '\t' << (*this)(i).mean() << '\t' << (*this)(i).standardDeviation() << '\n';
		}
		if(format().style()==phys_tools::gnuplot::Format::LEFT_STEPS
		   || format().style()==phys_tools::gnuplot::Format::RIGHT_STEPS)
			data << getBinEdge(0,count) << '\t' << 0.0 << '\n';
		data << 'e' << std::endl;
		return(data.str());
	}
	
	std::string Plottable1DErrorHistogram::get_data_code (const double& x_min,
												   const double& x_max,
												   const double& y_min,
												   const double& y_max) const{
		return(get_data_code());
	}
	
	std::string Plottable1DErrorHistogram::get_plot_code () const{
        std::ostringstream plot_code;
        plot_code << "'-' title \"" << _title << "\" ";
        plot_code << _format;
        return (plot_code.str());
    }*/

} //namespace phys_tools