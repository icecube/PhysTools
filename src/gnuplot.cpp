// gnuplot.cpp
// Access Gnuplot directly and easily from C++ code


#include "PhysTools/gnuplot.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <boost/filesystem.hpp>

namespace phys_tools {
namespace gnuplot {

namespace fs = boost::filesystem;


// Range implementation ---------------------------------------------

Range::Range (const std::string& name)
:   name_ (name)
{
}

Range::Range (const std::string& name, const double min, const double max)
:   min_ (min), max_ (max), name_ (name)
{
}

Range::~Range ()
{
}

Range&
Range::auto_min ()
{
    min_.clear ();
    return *this;
}

Range&
Range::auto_max ()
{
    max_.clear ();
    return *this;
}

Range&
Range::min (const double min)
{
    min_ = min;
    return *this;
}

Value<double>
Range::min () const
{
    return min_;
}

Range&
Range::max (const double max)
{
    max_ = max;
    return *this;
}

Value<double>
Range::max () const
{
    return max_;
}

const std::string&
Range::name() const{
    return name_;
}

std::ostream&
operator<< (std::ostream& stream, const Range& range)
{
    stream << "set " << range.name_ << "range " << "[";
    if (range.min_.empty ()) {
        stream  << "*";
    }
    else {
        stream  << range.min_;
    }
    stream << ":";
    if (range.max_.empty ()) {
        stream  << "*";
    }
    else {
        if (!range.min_.empty () && range.min_ == range.max_) {
            stream  << range.max_ + .1;
        }
        else {
            stream  << range.max_;
        }
    }
    stream << "] writeback" << std::endl;
    return stream;
}


// Axis implementation ----------------------------------------------

Axis::Axis (const std::string& name, const std::string& label)
: Range (name), label_ (label), zeroaxis_ (false), logscale_(false)
{
}

Axis::Axis (const std::string& name, const std::string& label,
            const double& min, const double& max,
            const bool zeroaxis)
: Range (name, min, max), label_ (label), zeroaxis_ (zeroaxis), logscale_(false)
{
}


Axis::~Axis ()
{
}

Axis&
Axis::auto_min ()
{
    Range::auto_min ();
    return *this;
}

Axis&
Axis::auto_max ()
{
    Range::auto_max ();
    return *this;
}

Axis&
Axis::label (const std::string& label)
{
    label_ = label;
    return *this;
}

Axis&
Axis::logscale (const bool logscale)
{
    logscale_ = logscale;
    return *this;
}

bool
Axis::logscale () const
{
    return logscale_;
}

Axis&
Axis::min (const double& min)
{
    Range::min (min);
    return *this;
}

Value<double>
Axis::min () const
{
    return min_;
}

Axis&
Axis::max (const double& max)
{
    Range::max (max);
    return *this;
}

Value<double>
Axis::max () const
{
    return max_;
}
    
Axis&
Axis::tic_start(const double& start)
{
    tic_start_=start;
    return *this;
}

Axis&
Axis::tic_increment(const double& increment)
{
    tic_increment_=increment;
    return *this;
}

Axis&
Axis::tic_end(const double& end)
{
    tic_end_=end;
    return *this;
}
    
Axis&
Axis::add_tic(const double& tic)
{
    tics_.push_back(tic);
    return(*this);
}

Axis&
Axis::minor_tics(unsigned int intervals){
    if(!intervals)
        throw std::domain_error("gnuplot::Axis::minor_tics: number of intervals must be >0");
    minor_tics_=intervals;
    return *this;
}
	
Axis&
Axis::tic_scale(double scale)
{
    tic_scale_=scale;
    return(*this);
}

Axis&
Axis::minor_tic_scale(double scale)
{
    minor_tic_scale_=scale;
    return(*this);
}
    
Axis&
Axis::format(std::string fmt){
    format_ = fmt;
    return(*this);
}
    
Axis&
Axis::unset_format(){
    format_ = std::string("% g");
    return(*this);
}


Axis&
Axis::zeroaxis (const bool zeroaxis)
{
    zeroaxis_ = zeroaxis;
    return *this;
}

std::ostream&
operator<< (std::ostream& stream, const Axis& axis)
{
	if(axis.logscale_){
		if(!axis.min_.empty() && axis.min_<=0)
			throw std::runtime_error("Zero or negative minimum for "+axis.name_+" axis with log scale");
		if(!axis.max_.empty() && axis.max_<=0)
			throw std::runtime_error("Zero or negative maximum for "+axis.name_+" axis with log scale");
	}
	
    stream << static_cast<const Range&> (axis);
    stream  << "set " << axis.name_ << "label "
        << "\"" << axis.label_ << "\"" << std::endl;
    if (!axis.logscale_)
        stream << "un";
    stream << "set logscale " << axis.name_ << std::endl;
    if(!axis.format_.empty())
        stream << "set format " << axis.name_ << " \"" << axis.format_ << '"' << std::endl;
    if(!axis.tics_.empty()){
        stream << "set " << axis.name_ << "tics (";
        for(unsigned int i=0; i<axis.tics_.size(); i++){
            if(i)
                stream << ',';
            stream << axis.tics_[i];
        }
        stream << ')' << std::endl;
    }
    else if(!axis.tic_start_.empty() && !axis.tic_increment_.empty()){
        stream << "set " << axis.name_ << "tics " << axis.tic_start_ << ',' << axis.tic_increment_;
        if(!axis.tic_end_.empty())
            stream << ',' << axis.tic_end_;
        stream << std::endl;
    }
    else if(!axis.tic_increment_.empty())
        stream << "set " << axis.name_ << "tics " << axis.tic_increment_ << std::endl;
    if(!axis.minor_tics_.empty())
        stream << "set m" << axis.name_ << "tics " << axis.minor_tics_ << std::endl;
    if (axis.name_ == "z" || axis.name_ == "cb")    // older gnuplot doesn't accept zeroaxis
        return stream;
    if (!axis.zeroaxis_)
        stream << "un";
    stream  << "set " << axis.name_ << "zeroaxis" << std::endl;
    if (!axis.tic_scale_.empty()){
        stream << "set " << axis.name_ << "tics scale " << axis.tic_scale_;
        if(!axis.minor_tic_scale_.empty())
            stream << ',' << axis.minor_tic_scale_;
        stream << std::endl;
    }
    else if(!axis.minor_tic_scale_.empty())
        stream << "set " << axis.name_ << "tics scale 1," << axis.minor_tic_scale_ << std::endl;
    return stream;
}
	
	
// ColorSpec implementation -----------------------------------------

ColorSpec rgb_color(const std::string& colorname){
	return(ColorSpec("rgbcolor \""+colorname+"\""));
}

ColorSpec rgb_color(uint32_t colorcode){
	std::ostringstream ss;
	ss << std::setw(6) << std::setfill('0') << std::hex << colorcode;
	if(colorcode>0xFFFFFF)
		throw std::range_error("rgb_color: Color code "+ss.str()+" out of range (must be 24 bits for RRGGBB)");
	return(ColorSpec("rgbcolor \"#"+ss.str()+"\""));
}

ColorSpec argb_color(uint32_t colorcode){
	std::ostringstream ss;
	ss << std::setw(8) << std::setfill('0') << std::hex << colorcode;
	return(ColorSpec("rgbcolor \"#"+ss.str()+"\""));
}

ColorSpec frac_color(double value){
	if(!(value>=0 && value<=1)) //use this form to trap NaNs
		throw std::range_error("frac_color: Palette fraction value "+std::to_string(value)+" out of range (must be in [0,1])");
	return(ColorSpec("palette frac "+std::to_string(value)));
}

ColorSpec cb_color(double value){
	if(std::isinf(value) || std::isnan(value))
		throw std::range_error("frac_color: Palette fraction value "+std::to_string(value)+" out of range (must be finite)");
	return(ColorSpec("palette cb "+std::to_string(value)));
}

std::ostream& operator<< (std::ostream& stream, const ColorSpec& color){
	return(stream << color.representation_);
}

// Coordinate implementation ----------------------------------------

std::ostream& operator<<(std::ostream& os, const Coordinate::CoordinateSystem& cs){
	switch(cs){
		case Coordinate::first:     return(os << "first");
		case Coordinate::second:    return(os << "second");
		case Coordinate::graph:     return(os << "graph");
		case Coordinate::screen:    return(os << "screen");
		case Coordinate::character: return(os << "character");
	}
	return(os);
}

std::ostream& operator<<(std::ostream& os, const Coordinate& c){
	return(os << c.xs << ' ' << c.x << ", " << c.ys << ' ' << c.y << ", " << c.zs << ' ' << c.z);
}

// Format implementation --------------------------------------------

Format::Format (const Style style)
:   style_ (style), nohidden3d_ (false), nocontours_ (false)
{
}

Format&
Format::style (const Style style)
{
    style_ = style;
    return *this;
}
Format::Style
Format::style() const{
    return style_;
}

Format&
Format::line_type (const int line_type)
{
    line_type_ = line_type;
    return *this;
}

Format&
Format::line_width (const float line_width)
{
    line_width_ = line_width;
    return *this;
}

Format&
Format::line_color (ColorSpec line_color)
{
    line_color_ = line_color;
    return *this;
}

Format&
Format::point_type (const int point_type)
{
    point_type_ = point_type;
    return *this;
}

Format&
Format::point_size (const float point_size)
{
    point_size_ = point_size;
    return *this;
}

Format&
Format::fill_style (const std::string fill_style)
{
    fill_style_ = fill_style;
    return *this;
}

Format&
Format::nohidden3d (const bool nohidden3d)
{
    nohidden3d_ = nohidden3d;
    return *this;
}

Format&
Format::nocontours (const bool nocontours)
{
    nocontours_ = nocontours;
    return *this;
}

std::ostream&
operator<< (std::ostream& stream, const Format& format)
{
    stream << "with ";
    switch (format.style_) {
    case Format::LINES:
        stream << "lines ";
        break;
    case Format::POINTS:
        stream << "points ";
        break;
    case Format::LINES_POINTS:
        stream << "linespoints ";
        break;
    case Format::IMPULSES:
        stream << "impulses ";
        break;
    case Format::DOTS:
        stream << "dots ";
        break;
    case Format::LEFT_STEPS:
        stream << "steps ";
        break;
    case Format::RIGHT_STEPS:
        stream << "fsteps ";
        break;
    case Format::CENTER_STEPS:
        stream << "histeps ";
        break;
    case Format::ERROR_BARS:
        stream << "errorbars ";
        break;
    case Format::X_ERROR_BARS:
        stream << "xerrorbars ";
        break;
    case Format::Y_ERROR_BARS:
        stream << "yerrorbars ";
        break;
    case Format::XY_ERROR_BARS:
        stream << "xyerrorbars ";
        break;
    case Format::ERROR_LINES:
        stream << "errorlines ";
        break;
    case Format::X_ERROR_LINES:
        stream << "xerrorlines ";
        break;
    case Format::Y_ERROR_LINES:
        stream << "yerrorlines ";
        break;
    case Format::XY_ERROR_LINES:
        stream << "xyerrorlines ";
        break;
    case Format::BOXES:
        stream << "boxes ";
        break;
    case Format::PM3D:
        stream << "pm3d ";
        break;
    case Format::FILLED_CURVES:
        stream << "filledcurves ";
            break;
    case Format::IMAGE:
        stream << "image ";
		break;
	case Format::IMAGE_FAILSAFE:
		stream << "image failsafe ";
		break;
    }
    if (format.line_type_.good ()) {
        stream << "linetype " << format.line_type_ << " ";
    };
    if (format.line_width_.good ()) {
        stream << "linewidth " << format.line_width_ << " ";
    };
    if (format.line_color_.good ()) {
        stream << "linecolor " << format.line_color_ << " ";
    };
    if (format.point_type_.good ()) {
        stream << "pointtype " << format.point_type_ << " ";
    };
    if (format.point_size_.good ()) {
        stream << "pointsize " << format.point_size_ << " ";
    };
    if (format.fill_style_.good ()) {
        stream << "fill " << format.fill_style_ << " ";
    };
    if (format.nohidden3d_) {
        stream << "nohidden3d ";
    };
    if (format.nocontours_) {
        stream << "nocontours ";
    };
    return stream;
}


// LineStyle implementation --------------------------------------------

    
LineStyle::LineStyle():
cycle_(false){}

LineStyle::LineStyle(int lw, ColorSpec lc, int pt):
cycle_(false),width_(lw),color_(lc),point_type_(pt){}

LineStyle& LineStyle::cycle(bool should_cycle){
    cycle_=should_cycle;
    return(*this);
}
LineStyle& LineStyle::width(int lw){
    width_=lw;
    return(*this);
}
LineStyle& LineStyle::color(ColorSpec lc){
    color_=lc;
    return(*this);
}
LineStyle& LineStyle::point_type(int pt){
    point_type_=pt;
    return(*this);
}
    
std::ostream& operator<<(std::ostream& stream, const LineStyle& style){
	//TODO: Is this broken? Why is it disabled?
    //if(style.cycle_)
    //    stream << "cycle ";
    if(style.color_.good())
        stream << "lc " << style.color_ << ' ';
    if(style.width_.good())
        stream << "lw " << style.width_ << ' ';
    if(style.point_type_.good())
        stream << " pt " << style.point_type_ << ' ';
    return(stream);
}


// Plottable implementation -----------------------------------------

Plottable::Plottable (const std::string& title, const Format& format)
:   format_ (format), title_ (title)
{
}

Plottable::~Plottable ()
{
}

Format&
Plottable::format ()
{
    return format_;
}
    
const Format&
Plottable::format () const
{
    return format_;
}

std::string
Plottable::get_data_code () const
{
    return std::string ();
}

std::string
Plottable::get_data_code (const double& x_min, const double& x_max,
                          const double& y_min, const double& y_max) const
{
    return this->get_data_code ();
}


// Plottable implementation -----------------------------------------

Plottable3D::Plottable3D (const std::string& title, const Format& format)
:   Plottable (title, format)
{
}

Plottable3D::~Plottable3D ()
{
}

std::string
Plottable3D::get_data_code () const
{
    return std::string ();
}


// FunctionPlot implementation --------------------------------------

FunctionPlot::FunctionPlot (const std::string& function,
                            const std::string& title,
                            const Format& format)
:   Plottable (title, format), function_ (function)
{
}

FunctionPlot::~FunctionPlot ()
{
}

std::string
FunctionPlot::get_plot_code () const
{
    std::ostringstream plot_code;
    plot_code << function_ << " ";
    plot_code << "title \"" << title_ << "\" ";
    plot_code << format_;
    return plot_code.str ();
}


// VerticalLinePlot implementation ----------------------------------

VerticalLinePlot::VerticalLinePlot (const double& x,
                                    const std::string& title,
                                    const Format& format)
:   Plottable (title, format), x_ (x)
{
}

VerticalLinePlot::~VerticalLinePlot ()
{
}

std::string
VerticalLinePlot::get_plot_code () const
{
    std::ostringstream plot_code;
    plot_code << "'-' title \"" << title_ << "\" ";
    plot_code << format_;
    return plot_code.str ();
}

std::string
VerticalLinePlot::get_data_code () const
{
    return "e\n";
}

std::string
VerticalLinePlot::get_data_code (const double& x_min, const double& x_max,
                                 const double& y_min, const double& y_max) const
{
    std::ostringstream stream;
    stream << x_ << "\t" << y_min << std::endl;
    stream << x_ << "\t" << y_max << std::endl;
    stream << "e" << std::endl;
    return stream.str ();
}
    

// DataPlot implementation --------------------------------------
    
DataPlot::~DataPlot(){}

std::string
DataPlot::get_plot_code() const{
    std::ostringstream plot_code;
    plot_code << "'-' title \"" << title_ << "\" ";
    plot_code << format_;
    return plot_code.str ();
}
    
std::string
DataPlot::get_data_code() const{
    std::ostringstream data_code;
    for(const point& p : data)
        data_code << p.x << '\t' << p.y << std::endl;
    data_code << "e" << std::endl;
    return data_code.str ();
}
    
std::string
DataPlot::get_data_code(const double& x_min,
                        const double& x_max,
                        const double& y_min,
                        const double& y_max) const{
    return(get_data_code());
}


// FunctionPlot3D implementation --------------------------------------

FunctionPlot3D::FunctionPlot3D (const std::string& function,
                                const std::string& title,
                                const Format& format)
:   Plottable3D (title, format), function_ (function)
{
}

FunctionPlot3D::~FunctionPlot3D ()
{
}

std::string
FunctionPlot3D::get_plot_code () const
{
    std::ostringstream plot_code;
    plot_code << function_ << " ";
    plot_code << "title \"" << title_ << "\" ";
    plot_code << format_;
    return plot_code.str ();
}
    

// Terminal implementation

std::ostream& operator<<(std::ostream& os, const Terminal& term){
    os << "set terminal " << term.name_ << ' ';
    term.write_settings_code(os);
    os << std::endl;
    return(os);
}


// PostscriptTerminal implementation

PostscriptTerminal::PostscriptTerminal():
Terminal("postscript"),
mode_(eps),
enhanced_(true),
plexity_(0),
color_(true),
solid_(true),
dash_length_scale_(1.0),
line_width_scale_(1.0),
x_size_(5.0),
y_size_(3.5),
black_text_(true),
font_("Helvetica"),
font_size_(14)
{}
    
void PostscriptTerminal::write_settings_code(std::ostream& os) const{
    switch(mode_){
        case landscape:
            os << "landscape ";
            break;
        case portrait:
            os << "portrait ";
            break;
        case eps:
            os << "eps ";
            break;
    }
    os << (enhanced_?"enhanced ":"noenhanced ");
    /*switch(plexity_){
        case 0:
            os << "defaultplex ";
            break;
        case 1:
            os << "simplex ";
            break;
        case 2:
            os << "duplex ";
            break;
        default:
            throw std::logic_error("Invalid plexity for PostscriptTerminal");
    }*/
    os << (color_?"color ":"monochrome ");
    os << (solid_?"solid ":"dashed ");
    //os << "dashlength " << dash_length_scale_ << ' ';
    os << "linewidth " << line_width_scale_ << ' ';
    os << "size " << x_size_ << ',' << y_size_ << ' ';
    //os << (black_text_?"blacktext ":"colortext ");
    os << "font \"" << font_ << ',' << font_size_ << '"';
}

PostscriptTerminal& PostscriptTerminal::mode(Mode m){
    mode_=m;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::enhanced(bool e){
    enhanced_=e;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::plexity(unsigned int plexity){
    if(plexity>2)
        throw std::logic_error("Invalid plexity for PostscriptTerminal");
    plexity_=plexity;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::color(bool c){
    color_=c;
    return(*this);
}
    
PostscriptTerminal& PostscriptTerminal::solid(bool s){
    solid_=s;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::dash_length(float len){
    dash_length_scale_=len;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::line_width(float width){
    line_width_scale_=width;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::x_size(float x){
    x_size_=x;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::y_size(float y){
    y_size_=y;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::black_text(bool b){
    black_text_=b;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::font(std::string f){
    font_=f;
    return(*this);
}

PostscriptTerminal& PostscriptTerminal::font_size(unsigned int s){
    font_size_=s;
    return(*this);
}

// Label implementation
    
Label::Label(std::string text, Coordinate position):
text_(text),position_(position),front_(false){}
    
Label&
Label::text(std::string t){
    text_=t;
    return(*this);
}
std::string
Label::text() const{
    return(text_);
}

Label&
Label::at(Coordinate pos){
    position_=pos;
    return(*this);
}
Label&
Label::position(Coordinate pos){
    position_=pos;
    return(*this);
}
Coordinate
Label::at() const{
    return(position_);
}
Coordinate
Label::position() const{
    return(position_);
}
    
Label&
Label::align(Label::Alignment a){
    text_align_=a;
    return(*this);
}
Value<Label::Alignment>
Label::align() const{
    return(text_align_);
}
    
Label&
Label::rotate(double angle){
    if(angle==0)
        rotation_.clear();
    else
        rotation_=angle;
    return(*this);
}
double
Label::rotate() const{
    return(rotation_.empty()?0.0:rotation_.value());
}

Label& Label::front(){
    front_=true;
    return(*this);
}

Label& Label::back(){
    front_=false;
    return(*this);
}

Label&
Label::font(std::string name, Value<unsigned int> size){
    font_name_=name;
    font_size_=size;
    return(*this);
}
Value<std::string>
Label::font() const{
    return(font_name_);
}
    
Label&
Label::font_size(unsigned int size){
    font_size_=size;
    return(*this);
}
Value<unsigned int>
Label::font_size() const{
    return(font_size_);
}

Label&
Label::font_color(ColorSpec color){
    font_color_=color;
    return(*this);
}
Label&
Label::text_color(ColorSpec color){
    font_color_=color;
    return(*this);
}
Value<ColorSpec>
Label::font_color() const{
    return(font_color_);
}
Value<ColorSpec>
Label::text_color() const{
    return(font_color_);
}
    
std::ostream& operator<<(std::ostream& stream, const Label& label){
    stream << '"' << label.text_ << '"';
    stream << " at " << label.position_;
    if(!label.text_align_.empty()){
        switch(label.text_align_.value()){
            case Label::LEFT: stream << " left"; break;
            case Label::CENTER: stream << " center"; break;
            case Label::RIGHT: stream << " right"; break;
        }
    }
    if(!label.rotation_.empty())
        stream << " rotate by " << label.rotation_;
    if(!label.font_name_.empty()){
        stream << " font \"" << label.font_name_;
        if(label.font_size_.empty())
            stream << ',' << label.font_size_;
        stream << '"';
    }
    stream << (label.front_?" front":" back");
    if(!label.font_color_.empty())
        stream << " textcolor " << label.font_color_;
    return(stream);
}


// Arrow implementation ---------------------------------------------

std::ostream& operator<<(std::ostream& os, const Arrow::HeadType& ht){
	switch(ht){
		case Arrow::nohead:   return(os << "nohead");
		case Arrow::head:     return(os << "head");
		case Arrow::backhead: return(os << "backhead");
		case Arrow::heads:    return(os << "heads");
	}
	return(os);
}

std::ostream& operator<<(std::ostream& os, const Arrow::HeadFill& ht){
	switch(ht){
		case Arrow::filled:   return(os << "filled");
		case Arrow::empty:    return(os << "empty");
		case Arrow::nofilled: return(os << "nofilled");
	}
	return(os);
}

std::ostream& operator<<(std::ostream& os, const Arrow& a){
	os << "from " << a.from_ << " to " << a.to_ << ' ' << a.ht_;
	if(a.head_length_.good() || a.head_angle_.good() || a.head_back_angle_.good()){
		os << " size screen "; //FIXME: hard-coded coordinate system
		if(a.head_length_.good())
			os << a.head_length_;
		else
			os << 0.017;
		if(a.head_angle_.good() || a.head_back_angle_.good()){
			os << ',';
			if(a.head_angle_.good())
				os << a.head_angle_;
			else
				os << 20;
			if(a.head_back_angle_.good())
				os << ',' << a.head_back_angle_;
		}
	}
	os << ' ' << a.hf_ << (a.front_?" front":" back");
	if(a.width_.good())
		os << " lw " << a.width_;
	return(os);
}

// Gnuplot implementation -------------------------------------------

Gnuplot::Gnuplot (const /*fs::path&*/std::string& path_to_gnuplot)
:   x_axis_ ("x"),
    y_axis_ ("y"),
    z_axis_ ("z"),
    cb_axis_ ("cb"),
    table_output_ (false),
    base_contours_ (false),
    surface_contours_ (false),
    pipe (commands, path_to_gnuplot/*.string()*/, outputs, errors)
{
    //figure out what kind of terminal we're using now
    commands << "show terminal" << std::endl;
    std::string line;
	static const std::string label="terminal type is ";
    while(getline(errors, line).good()){
		std::size_t pos=line.find(label);
		if(pos!=std::string::npos){
			set_terminal<Terminal>(line.substr(pos+label.size()));
			break;
		}
    }
    errors.clear();
}

Gnuplot::~Gnuplot ()
{
    commands << "exit" << std::endl;
    pipe.wait ();
}

std::string
Gnuplot::title ()
{
    return title_;
}

void
Gnuplot::title (const std::string& title)
{
    title_ = title;
}
    
void Gnuplot::manual_settings(const std::string& settings){
    manual_settings_ = settings;
}

Axis&
Gnuplot::x_axis ()
{
    return x_axis_;
}

Axis&
Gnuplot::y_axis ()
{
    return y_axis_;
}

Axis&
Gnuplot::z_axis ()
{
    return z_axis_;
}
    
Axis&
Gnuplot::cb_axis ()
{
    return cb_axis_;
}
    
LineStyle&
Gnuplot::line_style(unsigned int num)
{
    if(num==0)
        throw std::logic_error("line style tags must be > 0");
    return(line_styles_[num]);
}
    
void
Gnuplot::unset_line_style(unsigned int num){
    if(num==0)
        throw std::logic_error("line style tags must be > 0");
    if(line_styles_.find(num)!=line_styles_.end())
        line_styles_.erase(line_styles_.find(num));
    commands << "unset linetype " << num << std::endl;
}
    
Label&
Gnuplot::label(unsigned int num)
{
    if(num==0)
        throw std::logic_error("label tags must be > 0");
    return(labels_[num]);
}

void
Gnuplot::unset_label(unsigned int num){
    if(num==0)
        throw std::logic_error("label tags must be > 0");
    if(labels_.find(num)!=labels_.end())
        labels_.erase(labels_.find(num));
    commands << "unset label " << num << std::endl;
}
    
unsigned int
Gnuplot::add_label(const Label& l){
    unsigned int idx=1;
    if(!labels_.empty()){
        if(labels_.rbegin()->first==labels_.size()) //used indices are dense
            idx=labels_.rbegin()->first+1;
        else //used indices are sparse; search for the first unused one the dumb way
            while(labels_.find(idx)!=labels_.end())
                idx++;
    }
    
    labels_[idx]=l;
    return(idx);
}
    
unsigned int
Gnuplot::add_label(Label&& l){
    unsigned int idx=1;
    if(!labels_.empty()){
        if(labels_.rbegin()->first==labels_.size()) //used indices are dense
            idx=labels_.rbegin()->first+1;
        else //used indices are sparse; search for the first unused one the dumb way
            while(labels_.find(idx)!=labels_.end())
                idx++;
    }
    
    labels_[idx]=std::move(l);
    return(idx);
}

Arrow&
Gnuplot::arrow(unsigned int num)
{
	if(num==0)
		throw std::logic_error("arrow tags must be > 0");
	return(arrows_[num]);
}

void
Gnuplot::unset_arrow(unsigned int num){
	if(num==0)
		throw std::logic_error("arrow tags must be > 0");
	if(arrows_.find(num)!=arrows_.end())
		arrows_.erase(arrows_.find(num));
	commands << "unset arrow " << num << std::endl;
}

unsigned int
Gnuplot::add_arrow(const Arrow& a){
	unsigned int idx=1;
	if(!arrows_.empty()){
		if(arrows_.rbegin()->first==arrows_.size()) //used indices are dense
			idx=arrows_.rbegin()->first+1;
		else //used indices are sparse; search for the first unused one the dumb way
			while(arrows_.find(idx)!=arrows_.end())
				idx++;
	}
	
	arrows_[idx]=a;
	return(idx);
}

unsigned int
Gnuplot::add_arrow(Arrow&& a){
	unsigned int idx=1;
	if(!arrows_.empty()){
		if(arrows_.rbegin()->first==arrows_.size()) //used indices are dense
			idx=arrows_.rbegin()->first+1;
		else //used indices are sparse; search for the first unused one the dumb way
			while(arrows_.find(idx)!=arrows_.end())
				idx++;
	}
	
	arrows_[idx]=std::move(a);
	return(idx);
}
    
void
Gnuplot::set_bar_size(double size){
    std::ostringstream ss;
    if(size<0)
        ss << "set bars fullwidth";
    else
        ss << "set bars " << size;
    error_bar_size=ss.str();
}

void
Gnuplot::set_bars_front(bool front){
    std::ostringstream ss;
    ss << "set bars " << (front?"front":"back");
    error_bar_front=ss.str();
}

void
Gnuplot::show_surface(bool show){
    show_surface_=show;
}
    
void
Gnuplot::show_grid(bool show){
    show_grid_=show;
}
    
void
Gnuplot::set_contour(bool base, bool surface){
    base_contours_=base;
    surface_contours_=surface;
}

/*void
Gnuplot::plot_eps (const Plottable& plottable,
                   const boost::filesystem::path filename,
                   const std::string fontname,
                   const unsigned fontsize,
                   const double x_size,
                   const double y_size,
                   const bool enable_color,
                   const bool solid_lines) const
{
    const Plottable* p = &plottable;
    plot_eps (&p, &p + 1,
              filename, fontname, fontsize, x_size, y_size,
              enable_color, solid_lines);
}
    
void
Gnuplot::plot_eps (std::initializer_list<Plottable*> plottables,
                   const boost::filesystem::path filename,
                   const std::string fontname,
                   const unsigned fontsize,
                   const double x_size,
                   const double y_size,
                   const bool enable_color,
                   const bool solid_lines) const
{
    plot_eps (plottables.begin(), plottables.end(),
              filename, fontname, fontsize, x_size, y_size,
              enable_color, solid_lines);
}

void
Gnuplot::plot_epslatex (const Plottable& plottable,
                        const boost::filesystem::path filename,
                        const bool standalone,
                        const double x_size,
                        const double y_size,
                        const bool enable_color,
                        const bool solid_lines) const
{
    const Plottable* p = &plottable;
    plot_epslatex (&p, &p + 1,
                   filename, standalone, x_size, y_size,
                   enable_color, solid_lines);
}

void
Gnuplot::plot_png (const Plottable& plottable,
                   const boost::filesystem::path filename,
                   const std::string fontname,
                   const unsigned fontsize,
                   const int width,
                   const int height) const
{
    const Plottable* p= &plottable;
    plot_png (&p, &p + 1, filename, fontname, fontsize, width, height);
}

void
Gnuplot::plot_wxt (const Plottable& plottable) const
{
    const Plottable* p= &plottable;
    plot_wxt (&p, &p + 1);
}

void
Gnuplot::plot_eps_3d (const Plottable3D& plottable,
                      const boost::filesystem::path filename,
                      bool useContours,
                      const std::string fontname,
                      const unsigned fontsize,
                      const double x_size,
                      const double y_size,
                      const bool enable_color,
                      const bool solid_lines) const
{
    const Plottable3D* p = &plottable;
    plot_eps_3d (&p, &p + 1,
                 filename, useContours, fontname, fontsize, x_size, y_size,
                 enable_color, solid_lines);
}

void
Gnuplot::plot_epslatex_3d (const Plottable3D& plottable,
                           const boost::filesystem::path filename,
                           const bool standalone,
                           const double x_size,
                           const double y_size,
                           const bool enable_color,
                           const bool solid_lines) const
{
    const Plottable3D* p = &plottable;
    plot_epslatex_3d (&p, &p + 1,
                   filename, standalone, x_size, y_size,
                   enable_color, solid_lines);
}

void
Gnuplot::plot_png_3d (const Plottable3D& plottable,
                      const boost::filesystem::path filename,
                      const std::string fontname,
                      const unsigned fontsize,
                      const int width,
                      const int height) const
{
    const Plottable3D* p= &plottable;
    plot_png_3d (&p, &p + 1, filename, fontname, fontsize, width, height);
}

void
Gnuplot::plot_wxt_3d (const Plottable3D& plottable) const
{
    const Plottable3D* p= &plottable;
    plot_wxt_3d (&p, &p + 1);
}

void
Gnuplot::plot_eps_contour (const Plottable3D& plottable,
                           const boost::filesystem::path filename,
                           const std::string fontname,
                           const unsigned fontsize,
                           const double x_size,
                           const double y_size,
                           const bool enable_color,
                           const bool solid_lines) const
{
    const Plottable3D* p= &plottable;
    std::ostringstream    terminal_code;
    terminal_code << "set term postscript eps enhanced ";
    if (enable_color) {
        terminal_code << "color ";
    }
    if (solid_lines) {
        terminal_code << "solid ";
    }
    terminal_code << "size " << x_size << ',' << y_size << ' ';
    terminal_code << '"' << fontname << "\" " << fontsize << ' ';
    terminal_code << std::endl;
    terminal_code << "set output \"" << filename.string() << '"';
    terminal_code << std::endl;
    plot_agnostic_3d (terminal_code.str (), &p, &p + 1, true);
}

void
Gnuplot::plot_epslatex_contour (const Plottable3D& plottable,
                                const boost::filesystem::path filename,
                                const bool standalone,
                                const double x_size,
                                const double y_size,
                                const bool enable_color,
                                const bool solid_lines) const
{
    const Plottable3D* p= &plottable;
    std::ostringstream    terminal_code;
    terminal_code << "set term epslatex ";
    if (standalone) {
        terminal_code << "standalone ";
    }
    if (enable_color) {
        terminal_code << "color ";
    }
    if (solid_lines) {
        terminal_code << "solid ";
    }
    terminal_code << "size " << x_size << ',' << y_size << ' ';
    terminal_code << std::endl;
    terminal_code << "set output \"" << filename.string() << '"';
    terminal_code << std::endl;
    plot_agnostic_3d (terminal_code.str (), &p, &p + 1, true);
}

void
Gnuplot::plot_png_contour (const Plottable3D& plottable,
                           const boost::filesystem::path filename,
                           const std::string fontname,
                           const unsigned fontsize,
                           const int width,
                           const int height) const
{
    const Plottable3D* p= &plottable;
    std::ostringstream    terminal_code;
    terminal_code << "set term png ";
    std::string basic_fonts[] = {"tiny", "small", "medium", "large", "giant"};
    if (find (basic_fonts, basic_fonts + 5, fontname) != basic_fonts + 5) {
        terminal_code << fontname << ' ';
    }
    else {
        terminal_code << "font \"" << fontname << "\" " << fontsize << " ";
    }
    terminal_code << "size " << width << ',' << height;
    terminal_code << std::endl;
    terminal_code << "set output \"" << filename.string() << '"';
    terminal_code << std::endl;
    plot_agnostic_3d (terminal_code.str (), &p, &p + 1, true);
}

void
Gnuplot::plot_wxt_contour (const Plottable3D& plottable) const
{
    const Plottable3D* p= &plottable;
    std::ostringstream    terminal_code;
    terminal_code << "set term wxt persist" << std::endl;
    plot_agnostic_3d (terminal_code.str (), &p, &p + 1, true);
}*/



std::string
Gnuplot::get_settings_code () const
{
    std::ostringstream ss;
    ss << "set title \"" << title_ << "\"" << std::endl;
    ss << x_axis_;
    ss << y_axis_;
    ss << z_axis_;
    ss << cb_axis_;
    ss << "set xyplane 0" << std::endl;
    for(const auto& lt : line_styles_){
        ss << "set linetype " << lt.first << ' ' << lt.second << std::endl;
        if(lt.second.cycle_)
            ss << "set linetype cycle " << lt.first << std::endl;
    }
    if(error_bar_size.good())
        ss << error_bar_size.value() << std::endl;
    if(error_bar_front.good())
        ss << error_bar_front.value() << std::endl;
    for(const auto& label : labels_)
        ss << "set label " << label.first << ' ' << label.second << std::endl;
    for(const auto& arrow : arrows_)
        ss << "set arrow " << arrow.first << ' ' << arrow.second << std::endl;
    if(show_surface_.good())
        ss << (show_surface_.value()?"set":"unset") << " surface" << std::endl;
    if(show_grid_.good())
        ss << (show_grid_.value()?"set":"unset") << " grid" << std::endl;
    if(base_contours_ || surface_contours_){
        ss << "set contour ";
        if(base_contours_ && surface_contours_)
            ss << "both";
        else if(base_contours_)
            ss << "base";
        else if(surface_contours_)
            ss << "surface";
        ss << std::endl;
        if(!contour_levels_.empty()){
            ss << "set cntrparam levels discrete ";
            for(auto it=contour_levels_.begin(), end=contour_levels_.end(); it!=end; it++){
                if(it!=contour_levels_.begin())
                    ss << ',';
                ss << *it;
            }
            ss << std::endl;
        }
    }
    ss << manual_settings_ << std::endl;
    return ss.str ();
}
    
void Gnuplot::begin_multiplot(unsigned int rows, unsigned int columns){
    if(!!rows != !!columns)
        throw std::logic_error("Both rows and columns or neither should be specified for multiplot");
    commands << "set multiplot";
    if(rows && columns)
        commands << " layout " << rows << ',' << columns;
    commands << std::endl;
}
    
void Gnuplot::end_multiplot(){
    commands << "unset multiplot" << std::endl;
}
    
void Gnuplot::interact(){
    std::string line;
    std::cout << "Interacting directly with gnuplot; type \"exit\" to end interactive mode." << std::endl;
    commands << "set print" << std::endl;
    while(true){
        std::cout << "gnuplot> ";
        std::cout.flush();
        getline(std::cin,line);
        if(line==".q" || line=="exit"){
            std::cout << "Ending interactive mode." << std::endl;
            break;
        }
        commands << line << " \n print '--mark--'" << std::endl;
        while(getline(errors,line) && line.find("--mark--")!=0)
            std::cout << line << std::endl;
    }
}
    
std::string Gnuplot::get_raw_output(){
    commands << "set print '-' \n print '--mark--' \n set print" << std::endl;
    std::string results,line;
    while(getline(outputs,line)){
        if(line.find("--mark--")==0)
            break;
        results+=line;
        results+="\n";
    }
    return(results);
}

std::string Gnuplot::get_raw_errors(){
    commands << "set print \n print '--mark--'" << std::endl;
    std::string results,line;
    while(getline(errors,line)){
        if(line.find("--mark--")==0)
            break;
        results+=line;
        results+="\n";
    }
    return(results);
}


std::ostream&
operator<< (std::ostream& stream, const Gnuplot::NormalizeResult& normalization)
{
    stream << "set xrange [" << normalization.x_min
        << ":" << normalization.x_max << "]" << std::endl;
    stream << "set yrange [" << normalization.y_min
        << ":" << normalization.y_max << "]" << std::endl;
    return stream;
}

} // namespace gnuplot
} // namespace phys_tools
