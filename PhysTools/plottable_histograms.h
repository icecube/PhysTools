#ifndef PHYS_TOOLS_PLOTTABLE_HISTOGRAMS_H
#define PHYS_TOOLS_PLOTTABLE_HISTOGRAMS_H

#include "PhysTools/histogram.h"
#include "PhysTools/bin_types.h"
#include "PhysTools/gnuplot.h"

namespace phys_tools{

class Plottable1DHistogram : public gnuplot::Plottable, public histograms::histogram<1>{
public:
	Plottable1DHistogram(const std::string& title,gnuplot::Format f):
	Plottable(title,f){}
	
	template<typename... Axes,
	         typename = typename std::enable_if<detail::all_args_convertible<histograms::axis&,typename std::add_lvalue_reference<Axes>::type...>::value>::type>
	Plottable1DHistogram(const std::string& title,gnuplot::Format f,Axes... axes):
	Plottable(title,f),histograms::histogram<1>(axes...){}
	
	Plottable1DHistogram(const std::string& title,gnuplot::Format f, const histograms::histogram<1>& h):
	Plottable(title,f),histograms::histogram<1>(h){}
		
	Plottable1DHistogram(const std::string& title,gnuplot::Format f, histograms::histogram<1>&& h):
	Plottable(title,f),histograms::histogram<1>(h){}
	
	virtual std::string get_data_code() const;
	
	virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const;
	
	virtual std::string get_plot_code () const;
	
	const std::string& getTitle() const{ return(title_); }
	void setTitle(const std::string& t){ title_=t; }
};

class Plottable2DHistogram : public gnuplot::Plottable3D, public histograms::histogram<2>{
public:
	template<typename... Axes>
	Plottable2DHistogram(const std::string& title,Axes... axes):
	Plottable3D(title,gnuplot::Format::PM3D),histograms::histogram<2>(axes...){}
	
	Plottable2DHistogram(const std::string& title,histograms::histogram<2>& h):
	Plottable3D(title,gnuplot::Format::PM3D),histograms::histogram<2>(h){}
	
	Plottable2DHistogram(const std::string& title,histograms::histogram<2>&& h):
	Plottable3D(title,gnuplot::Format::PM3D),histograms::histogram<2>(h){}
	
	virtual std::string get_data_code() const;
	
	virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const;
	
	virtual std::string get_plot_code () const;
	
	const std::string& getTitle() const{ return(title_); }
	void setTitle(const std::string& t){ title_=t; }
};

//TODO: This needs work
class PlottableErrorBars : public gnuplot::Plottable{
public:
	const histograms::histogram<1>& values;
	histograms::histogram<1> errors;

	PlottableErrorBars(const histograms::histogram<1>& v, histograms::histogram<1>&& e):
	Plottable("",gnuplot::Format::ERROR_BARS),values(v),errors(e){}
	
	virtual std::string get_data_code() const;
	
	virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const;
	
	virtual std::string get_plot_code () const;
	
	//When plotting multiple sets of points with error bars it is oftern useful to be able to
	//offset each data series horizontally so that all are clearly visible.
	//Strictly speaking this makes the plotted x values incorrect, but if this is tolerable
	//the plot can become much easier to read.
	//The offset specified here is treated as a fraction of the width of the bin.
	//(This helps for offsetting along a logarithmic axis.)
	void set_offset(double offset){ horizontal_offset=offset; }
private:
	double horizontal_offset;
};
	
template<typename T>
class ErrorValueTraits{
public:
	static double value(const T& t){ return(0); }
	static double error(const T& t){ return(0); }
};
	
template<typename T>
class ErrorValueTraits<histograms::meanVarTracker<T>>{
public:
	static double value(const histograms::meanVarTracker<T>& t){ return(t.mean()); }
	static double error(const histograms::meanVarTracker<T>& t){ return(t.standardDeviation()); }
};
	
template<typename T>
class ErrorValueTraits<histograms::sqErrorValue<T>>{
public:
	static double value(const histograms::sqErrorValue<T>& t){ return(t.value()); }
	static double error(const histograms::sqErrorValue<T>& t){ return(t.standardDeviation()); }
};
	
template<typename T>
class Plottable1DErrorHistogram : public gnuplot::Plottable, public histograms::histogram<1,T>{
public:
	template<typename... Axes>
	Plottable1DErrorHistogram(const std::string& title,gnuplot::Format f,Axes... axes):
	Plottable(title,f),histograms::histogram<1,T>(axes...),horizontal_offset(0){}
	
	Plottable1DErrorHistogram(const std::string& title,gnuplot::Format f, histograms::histogram<1,T>&& h):
	Plottable(title,f),histograms::histogram<1,T>(h),horizontal_offset(0){}
	
	using histograms::histogram<1,T>::getBinCount;
	using histograms::histogram<1,T>::getBinEdge;
	using histograms::histogram<1,T>::getBinWidth;
	
	virtual std::string get_data_code() const{
		unsigned int count=getBinCount(0);
		std::ostringstream data;
		if(format().style()==gnuplot::Format::LEFT_STEPS
		   || format().style()==gnuplot::Format::RIGHT_STEPS)
			data << getBinEdge(0,0) << '\t' << 0.0 << '\n';
		for(unsigned int i=0; i<count; i++){
			double x=0.0;
			switch(format().style()){
				case gnuplot::Format::LEFT_STEPS:
					x=getBinEdge(0,i);
					break;
				case gnuplot::Format::RIGHT_STEPS:
					x=getBinEdge(0,i)+getBinWidth(0,i);
					break;
				case gnuplot::Format::CENTER_STEPS:
				case gnuplot::Format::Y_ERROR_BARS:
				case gnuplot::Format::Y_ERROR_LINES:
					x=histograms::histogram<1,T>::getBinCenter(0,i);
					break;
				default:
					throw std::runtime_error("Unsupported drawing format");
			}
			if(horizontal_offset)
				x+=horizontal_offset*getBinWidth(0,i);
			data << x << '\t' << ErrorValueTraits<T>::value((*this)(i)) << '\t' << ErrorValueTraits<T>::error((*this)(i)) << '\n';
		}
		if(format().style()==gnuplot::Format::LEFT_STEPS
		   || format().style()==gnuplot::Format::RIGHT_STEPS)
			data << getBinEdge(0,count) << '\t' << 0.0 << '\n';
		data << 'e' << std::endl;
		return(data.str());
	}
	
	virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const{
		return(get_data_code());
	}
	
	virtual std::string get_plot_code () const{
        std::ostringstream plot_code;
        plot_code << "'-' title \"" << title_ << "\" ";
        plot_code << format_;
        return (plot_code.str());
    }
	
	const std::string& getTitle() const{ return(title_); }
	void setTitle(const std::string& t){ title_=t; }
	
	//When plotting multiple sets of points with error bars it is oftern useful to be able to
	//offset each data series horizontally so that all are clearly visible.
	//Strictly speaking this makes the plotted x values incorrect, but if this is tolerable
	//the plot can become much easier to read.
	//The offset specified here is treated as a fraction of the width of the bin.
	//(This helps for offsetting along a logarithmic axis.)
	void set_offset(double offset){ horizontal_offset=offset; }
private:
	double horizontal_offset;
};
	
template<typename BinType>
class FilledErrorHist : public histograms::histogram<1,BinType>, public gnuplot::Plottable{
public:
	class Format : public gnuplot::Format{
	private:
		Value<gnuplot::ColorSpec>  fill_color_;
		gnuplot::Format::Style error_style_;
		bool showFill_, showErrors_;
	public:
		Format():gnuplot::Format(FILLED_CURVES),error_style_(LEFT_STEPS),showFill_(true),showErrors_(true){}
		Format(gnuplot::Format::Style s):gnuplot::Format(s),error_style_(LEFT_STEPS),showFill_(true),showErrors_(true){}
		Format(gnuplot::Format f):gnuplot::Format(f),error_style_(LEFT_STEPS),showFill_(true),showErrors_(true){}
		Format(const Format& f):gnuplot::Format(f),fill_color_(f.fill_color_),error_style_(f.error_style_),showFill_(f.showFill_),showErrors_(f.showErrors_){}
		
		gnuplot::ColorSpec fill_color() const{ return(fill_color_); }
		Format& fill_color(gnuplot::ColorSpec fill_color){
			fill_color_=fill_color;
			return(*this);
		}
		
		gnuplot::Format::Style error_style() const{ return(error_style_); }
		Format& error_style(gnuplot::Format::Style es){
			error_style_=es;
			return(*this);
		}
		
		bool showFill() const{ return(showFill_); }
		Format& showFill(bool show){
			showFill_=show;
			return(*this);
		}
		
		bool showErrors() const{ return(showErrors_); }
		Format& showErrors(bool show){
			showErrors_=show;
			return(*this);
		}
		
		Format makeFillFormat() const{
			Format f(*this);
			if(style()==FILLED_CURVES && fill_color_.good())
				f.line_color(fill_color_);
			return(f);
		}
		Format makeErrorFormat() const{
			Format f(*this);
			f.style(error_style_).point_size(0);
			return(f);
		}
	};
private:
	Format format_;
public:
	using iterator=typename histograms::histogram<1,BinType>::iterator;
	using const_iterator=typename histograms::histogram<1,BinType>::const_iterator;

	template<typename AxisType, typename = typename std::is_convertible<typename std::add_lvalue_reference<AxisType>::type,histograms::axis&>::type>
	FilledErrorHist(const std::string& title,Format f,AxisType axis):
	Plottable(title),histograms::histogram<1,BinType>(axis),format_(f){}
	
	FilledErrorHist(const FilledErrorHist& other):
	Plottable(other.title()),
	histograms::histogram<1,BinType>((const histograms::histogram<1,BinType>&)other),
	format_(other.format_){}
	
	FilledErrorHist(histograms::histogram<1,BinType>&& other):histograms::histogram<1,BinType>(std::move(other)){}
	
	FilledErrorHist(const std::string& title, histograms::histogram<1,BinType>&& other):Plottable(title),histograms::histogram<1,BinType>(std::move(other)){}
	
	Format& format(){ return(format_); }
	const Format& format() const{ return(format_); }
	
	histograms::histogram<1> getCenter(){
		using namespace histograms;
		histogram<1> center;
		center.setAxis(0,this->getAxis(0)->copy());
		bool saved=this->getUseContentScaling();
		this->setUseContentScaling(false);
		center.setUseContentScaling(false);
		for(const_iterator it=this->cbegin(), end=this->cend(); it!=end; it++)
			center.add(it.getBinCenter(0),amount((double)*it));
		this->setUseContentScaling(saved);
		center.setUseContentScaling(saved);
		return(center);
	}
	
	histograms::histogram<1> getMin(){
		using namespace histograms;
		histogram<1> min;
		min.setAxis(0,this->getAxis(0)->copy());
		bool saved=this->getUseContentScaling();
		this->setUseContentScaling(false);
		min.setUseContentScaling(false);
		for(const_iterator it=this->cbegin(), end=this->cend(); it!=end; it++)
			min.add(it.getBinCenter(0),amount((*it).errorMin()));
		this->setUseContentScaling(saved);
		min.setUseContentScaling(saved);
		return(min);
	}
	
	histograms::histogram<1> getMax(){
		using namespace histograms;
		histogram<1> max;
		max.setAxis(0,this->getAxis(0)->copy());
		bool saved=this->getUseContentScaling();
		this->setUseContentScaling(false);
		max.setUseContentScaling(false);
		for(const_iterator it=this->cbegin(), end=this->cend(); it!=end; it++)
			max.add(it.getBinCenter(0),amount((*it).errorMax()));
		this->setUseContentScaling(saved);
		max.setUseContentScaling(saved);
		return(max);
	}
	
	virtual std::string get_plot_code() const{
        std::ostringstream plot_code;
		if(format_.style()==gnuplot::Format::FILLED_CURVES){
			if(format_.showFill())
				plot_code << "'-' title \"" << title_ << "\" " << format_.makeFillFormat();
			if(format_.showFill() && format_.showErrors())
				plot_code << ", ";
			if(format_.showErrors()){
				plot_code << "'-' notitle " << format_.makeErrorFormat();
				plot_code << ", '-' notitle " << format_.makeErrorFormat();
			}
		}
		else{
			if(format_.showFill())
				plot_code << "'-' title \"" << title_ << "\" " << format_.makeFillFormat();
			if(format_.showFill() && format_.showErrors())
				plot_code << ", ";
			if(format_.showErrors())
				plot_code << "'-' notitle " << format_.makeErrorFormat();
		}
        return (plot_code.str());
    }
	
private:
	std::string get_data_code_core(const double x_min=-std::numeric_limits<double>::infinity(),
									const double x_max=std::numeric_limits<double>::infinity(),
									const double y_min=-std::numeric_limits<double>::infinity(),
									const double y_max=std::numeric_limits<double>::infinity()) const{
		unsigned int count=this->getBinCount(0);
		std::ostringstream data;
		auto clamp=[y_min,y_max](double y){
			if(y<y_min)
				return(y_min);
			if(y>y_max)
				return(y_max);
			return(y);
		};
		
		if(format_.style()==gnuplot::Format::FILLED_CURVES){
			auto outputErrors=[&,this](int bin, int item){
				const double min=clamp((*this)(bin).errorMin());
				const double max=clamp((*this)(bin).errorMax());
				
				switch(item){
					case 0: data << min << '\t' << max; break;
					case 1: data << min; break;
					case 2: data << max; break;
					
				}
				data << '\n';
			};
			
			const bool doublePoints=format_.error_style()==gnuplot::Format::LEFT_STEPS;
			auto outputPoint=
			[&,this](int bin, int item){
				data << this->getBinEdge(0,bin) << '\t';
				outputErrors(bin,item);
				if(doublePoints){
					data << this->getBinEdge(0,bin)+this->getBinWidth(0,bin) << '\t';
					outputErrors(bin,item);
				}
			};
			if(format_.showFill()){
				for(unsigned int i=0; i<count; i++)
					outputPoint(i,0);
				data << 'e' << std::endl;
			}
			if(format_.showErrors()){
				for(int item=1; item<3; item++){
					for(unsigned int i=0; i<count; i++)
						outputPoint(i,item);
					data << 'e' << std::endl;
				}
			}
		}
		else{
			const auto cStyle=format_.style();
			auto outputData=[&](int item){
				for(unsigned int i=0; i<count; i++){
					double x=0.0;
					switch(cStyle){
						case gnuplot::Format::LEFT_STEPS:
							x=this->getBinEdge(0,i);
							break;
						case gnuplot::Format::RIGHT_STEPS:
							x=this->getBinEdge(0,i)+this->getBinWidth(0,i);
							break;
						case gnuplot::Format::CENTER_STEPS:
						case gnuplot::Format::LINES:
						case gnuplot::Format::POINTS:
						case gnuplot::Format::LINES_POINTS:
						case gnuplot::Format::IMPULSES:
						case gnuplot::Format::DOTS:
						case gnuplot::Format::ERROR_BARS:
						case gnuplot::Format::Y_ERROR_BARS:
							//case gnuplot::Format::XY_ERROR_BARS:
							x=this->getBinCenter(0,i);
							break;
						default:
							throw std::runtime_error("Unsupported drawing format: "+boost::lexical_cast<std::string>(format().style()));
					}
					data << x;
					if(item==0)
						data << '\t' << (double)(*this)(i);
					if(item)
						data << '\t' << clamp((double)(*this)(i)) << '\t' << clamp((*this)(i).errorMin()) << '\t' << clamp((*this)(i).errorMax());
					data << '\n';
				}
				data << 'e' << std::endl;
			};
			if(format_.showFill())
				outputData(0);
			if(format_.showErrors())
				outputData(1);
		}
		return(data.str());
	}
public:
	virtual std::string get_data_code() const{
		return(get_data_code_core());
	}
	
	virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const{
		return(get_data_code_core(x_min,x_max,y_min,y_max));
	}
};

} //namespace phys_tools

#endif