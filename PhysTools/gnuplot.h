/// @file gnuplot.hpp
/// @brief Access Gnuplot directly and easily from C++ code


#ifndef GNUPLOT_HPP
#define GNUPLOT_HPP 1

#include <string>
#include <sstream>
#include <map>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>

#include "pipe.h"
#include "value.h"

namespace phys_tools {

/// @ingroup System
/// @brief Namespace for the Gnuplot interface.
///
/// This namespace contains classes for interacting with Gnuplot.  It
/// tries to provide an intuitive interface with sensible defaults,
/// rather than a feature-complete interface to the program.
///
/// Just getting started?  Try Gnuplot::plot_wxt().
namespace gnuplot {

// 2D stuff -------------------------------------------------

/// @brief Represent a range.
///
/// Range represents a range in Gnuplot, such as r or t.  Perhaps more
/// importantly, it also is the base class for Axis.
///
/// The main task for Range is to encapsulate the details of ranges
/// which are determined partially or completely automatically by
/// Gnuplot at plot time.
class Range {
public:
    explicit Range (const std::string& name);
    Range (const std::string& name, double min, double max);
    virtual ~Range ();

    /// @brief Set the range lower bound to autoscale.
    /// @return A reference to the Range.
    Range&          auto_min ();

    /// @brief Set the range upper bound to autoscale.
    /// @return A reference to the Range.
    Range&          auto_max ();

    /// @brief Set the range lower bound.
    /// @return A reference to the Range.
    Range&          min (const double min);

    /// @brief Get the range lower bound.
    /// @return The lower bound of the Range.
    Value<double>   min () const;

    /// @brief Set the range upper bound.
    /// @return A reference to the Range.
    Range&          max (const double max);

    /// @brief Get the range lower bound.
    /// @return The lower bound of the Range.
    Value<double>   max () const;
    
    const std::string& name() const;

    /// @brief Write to @c stream the Gnuplot code corresponding to the Range.
    friend std::ostream& operator<< (std::ostream& stream, const Range& range);

protected:
    Value<double>   min_;
    Value<double>   max_;
    std::string     name_;
};


/// @brief Represent an axis.
///
/// Axis represents an axis in Gnuplot, such as x or y.  It derives
/// its ability to handle automatic ranges from Range, but it adds
/// axis-specific functionality like labels and the 'zeroaxis' property.
class Axis : Range {
    friend class Gnuplot;
public:
    Axis (const std::string& name, const std::string& label = "");
    Axis (const std::string& name, const std::string& label,
          const double& min, const double& max,
          const bool zeroaxis = false);
    virtual ~Axis ();

    /// @brief Set the axis lower bound to autoscale.
    /// @return A reference to the Axis.
    Axis&           auto_min ();

    /// @brief Set the axis upper bound to autoscale.
    /// @return A reference to the Axis.
    Axis&           auto_max ();

    /// @brief Set the axis label.
    /// @return A reference to the Axis.
    Axis&           label (const std::string& label);

    /// @brief Set whether the axis is 'logscale'.
    /// @return A reference to the Axis
    Axis&           logscale (const bool logscale);

    /// @brief Get whether the axis is 'logscale'.
    /// @return Whether the Axis is set to 'logscale'.
    bool            logscale () const;

    /// @brief Set the axis lower bound.
    /// @return A reference to the Axis.
    Axis&           min (const double& min);

    /// @brief Get the axis lower bound.
    /// @return The lower bound of the Axis.
    Value<double>   min () const;

    /// @brief Set the axis upper bound.
    /// @return A reference to the Axis.
    Axis&           max (const double& max);

    /// @brief Get the axis upper bound.
    /// @return The upper bound of the Axis.
    Value<double>   max () const;
    
    /// @brief Set the base point for the major tics.
    /// @param start The point at which the first tic mark should be placed.
    /// @return A reference to the Axis.
    Axis& tic_start(const double& start);
    /// @brief Set spacing of the major tics.
    /// @param increment The distance between major tics (or ratio on a log scale).
    /// @return A reference to the Axis.
    Axis& tic_increment(const double& increment);
    /// @brief Set the end point for the major tics.
    /// @param end The point at which the last tic mark should be placed.
    /// @return A reference to the Axis.
    Axis& tic_end(const double& end);
    /// @brief Add a manually specifed tic mark at an arbitrary location.
    /// @param tic The location of the new tic.
    /// @return A reference to the Axis.
    Axis& add_tic(const double& tic);
    /// @brief Set the frequency of minor tics.
    /// @param intervals The number of intervals between minor tics, which
    ///                  must be >0.
    /// @return A reference to the Axis.
    Axis& minor_tics(unsigned int intervals);
    /// @brief set the length of the tic marks
    ///
    /// Unless explicitly set, the minor tic scale will be half of the major tic
    /// scale.
    /// @param scale the factor by which the default tic mark length is scaled
    Axis& tic_scale(double scale);
    /// @brief set the length of the minor tic marks
    /// @param scale the factor by which the default tic mark length is scaled
    Axis& minor_tic_scale(double scale);
    
    Axis& format(std::string fmt);
    Axis& unset_format();

    /// @brief Set the zeroaxis property.
    /// @return A reference to the Axis.
    ///
    /// If this property is set, a visible line will be drawn for the
    /// axis, going through the origin.
    Axis&           zeroaxis (const bool zeroaxis=true);

    /// @brief Write to @c stream the Gnuplot code corresponding to the Axis.
    friend std::ostream& operator<< (std::ostream& stream, const Axis& axis);
private:
    std::string     label_;
    bool            zeroaxis_;
    bool            logscale_;
    Value<double>   tic_start_;
    Value<double>   tic_end_;
    Value<double>   tic_increment_;
    Value<unsigned int> minor_tics_;
    std::vector< Value<double> > tics_;
    Value<std::string> format_;
    Value<double>   tic_scale_;
    Value<double>   minor_tic_scale_;
};

/// @brief Encapsulate a Gnuplot color specification.
///
/// Gnuplot allows colors to be specified in a number of ways.
/// This class acts as an opaque container for these various
/// representations. Actual ColorSpec objects should be produced
/// by the associated factory functions (@ref rgb_color, @ref argb_color,
/// @ref frac_color, @ref cb_color, and @ref z_color) which deal with the
/// boilerplate needed by each representation and provide type
/// checking.
class ColorSpec{
private:
	ColorSpec()=default; //hack because Value does not have a proper, true empty state
	ColorSpec(std::string r):representation_(r){}
	std::string representation_;
	
	friend ColorSpec rgb_color(const std::string& colorname);
	friend ColorSpec rgb_color(uint32_t colorcode);
	friend ColorSpec argb_color(uint32_t colorcode);
	friend ColorSpec frac_color(double value);
	friend ColorSpec cb_color(double value);
	friend ColorSpec z_color();
	
	friend std::ostream& operator<< (std::ostream& stream, const ColorSpec& color);
	friend class Value<ColorSpec>;
};

/// @brief Get an RGB color by name.
///
/// @param colorname A colorname recongnized by Gnuplot. (The list of
///                  such names can be obtained using the
///                  "show colornames" command.)
ColorSpec rgb_color(const std::string& colorname);

/// @brief Get an RGB color by hexadecimal representation.
///
/// @param colorcode An integer whose value defines the color when
///                  interpreted as 'RRGGBB' in hexadecimal. This
///                  value may not be larger than 0xFFFFFF.
ColorSpec rgb_color(uint32_t colorcode);

/// @brief Get an ARGB color by hexadecimal representation.
///
/// @param colorcode An integer whose value defines the color when
///                  interpreted as 'AARRGGBB' in hexadecimal.
ColorSpec argb_color(uint32_t colorcode);
	
//TODO: figure out how rgb variable colors work

/// @brief Get a color from the 'fractional' palette.
///
/// This function provides a color which is looked up from the active
/// color palette by a fractional value in [0,1] onto the full range
/// of the color palette.
/// @param value The fractional position on the palette from which
///                  to select the color. This value must be in [0,1].
ColorSpec frac_color(double value);

/// @brief Get a color from the color axis palette.
///
/// This function provides a color which is looked up from the active
/// color palette according to where it falls on the current cbrange
/// (cbaxis).
/// @param value The coordinate on the colorbar axis from which to
///              look up the color.
ColorSpec cb_color(double value);

/// @brief Get a color which depends on the z coordinate of the plotted
///        object.
///
/// This function provides a special, variable color which varies with the z
/// coordinate of the plotted object. (That is, the color is determined by
/// looking up the current z coordinate on the colorbar axis.)
ColorSpec z_color();

//TODO: figure out how variable colors work
    
class Coordinate{
public:
    enum CoordinateSystem{first,second,graph,screen,character};
private:
    double x, y, z;
    CoordinateSystem xs, ys, zs;
public:
	Coordinate():
	x(0),y(0),z(0),xs(first),ys(first),zs(first){}
	
	Coordinate(CoordinateSystem xs, double x, CoordinateSystem ys, double y, CoordinateSystem zs=first, double z=0.):
	x(x),y(y),z(z),xs(xs),ys(ys),zs(zs){}
	
	Coordinate(double x, CoordinateSystem ys, double y, CoordinateSystem zs=first, double z=0.):
	x(x),y(y),z(z),xs(first),ys(ys),zs(zs){}
    
	Coordinate(CoordinateSystem xs, double x, double y, CoordinateSystem zs=first, double z=0.):
	x(x),y(y),z(z),xs(xs),ys(first),zs(zs){}
    
	Coordinate(double x, double y, CoordinateSystem zs=first, double z=0.):
	x(x),y(y),z(z),xs(first),ys(first),zs(zs){}
    
	Coordinate(double x, double y, double z):
	x(x),y(y),z(z),xs(first),ys(first),zs(first){}
	
	friend std::ostream& operator<<(std::ostream&, const Coordinate&);
};
	
std::ostream& operator<<(std::ostream&, const Coordinate::CoordinateSystem&);
std::ostream& operator<<(std::ostream&, const Coordinate&);

/// @brief Encapsulate a Gnuplot plot format specification.
///
/// Gnuplot enables very fine control over the formatting of a
/// function or data plot.  Some of this control is, unfortunately,
/// output-format- or system-dependent.  Format is a thin
/// abstraction layer that makes it fairly easy to give a plot format
/// specification, provided that you are familiar with the underlying
/// Gnuplot settings.
//
// @sa Try entering "help plot with" at the Gnuplot prompt for more
// information.
class Format {
public:

    /// @brief A Gnuplot line or point type specification.
    ///
    /// Most of these correspond exactly to the literals used directly
    /// in Gnuplot, modulo extra underscores for readibility.
    /// Exceptions to this rule include LEFT_STEPS (Gnuplot: "steps"),
    /// RIGHT_STEPS (Gnuplot "fsteps"), and CENTER_STEPS (Gnuplot:
    /// "histeps").  Also note that "labels", "financebars",
    /// "candlesticks", "vectors", and "rgbimage" are omitted.
    enum Style {
        LINES, POINTS, LINES_POINTS, IMPULSES, DOTS,
        LEFT_STEPS, RIGHT_STEPS, CENTER_STEPS,
        ERROR_BARS, X_ERROR_BARS, Y_ERROR_BARS, XY_ERROR_BARS,
        ERROR_LINES, X_ERROR_LINES, Y_ERROR_LINES, XY_ERROR_LINES,
        BOXES, PM3D, FILLED_CURVES, IMAGE, IMAGE_FAILSAFE
    };

    /// @brief Construct a Format object.
    ///
    /// This constructor only allows the specification of a line or
    /// point style  The other methods all return a reference to the
    /// Format, so method chaining can be used to set the other
    /// properties conveniently.  For example:
    /// @code
    /// Format f1 (Format::LINES_POINTS);
    /// f1.line_type (2).line_color ("rgbcolor 'red'");
    ///
    /// Format f2 = Format (Format::BOXES).line_width (3).fill_style ("solid");
    /// @endcode
    Format (const Style style = LINES);

    /// @brief Set the line or point style.
    /// @param style The line or point style.
    /// @return A reference to the Format.
    Format&             style (const Style style);
    Style style() const;

    /// @brief Set the line type.
    /// @param line_type The line type.
    /// @return A reference to the Format.
    Format&             line_type (const int line_type);

    /// @brief Set the line width.
    /// @param line_width The line width.
    /// @return A reference to the Format.
    Format&             line_width (const float line_width);

    /// @brief Set the line color
    /// @param line_color The line color specification.
    /// @return A reference to the Format.
    Format&             line_color (ColorSpec line_color);

    /// @brief Set the point type.
    /// @param point_type The point type.
    /// @return A reference to the Format.
    Format&             point_type (const int point_type);

    /// @brief Set the point size.
    /// @param point_size The point size.
    /// @return A reference to the Format.
    Format&             point_size (const float point_size);

    /// @brief Set the fill style.
    /// @param fill_style The fill style.
    /// @return A reference to the Format.
    ///
    /// Of the Gnuplot styles that can take a fill style, we only have
    /// implemented BOXES (Gnuplot: "boxes").  Try "help fillstyle" in
    /// Gnuplot for more information.
    Format&             fill_style (const std::string fill_style);

    /// @brief Set the 'nohidden3d' property.
    /// @param nohidden3d Whether to disable hidden3d.
    /// @return A reference to the Format.
    ///
    /// In Gnuplot 3D plots, you can specify hidden3d, which causes
    /// things that are behind other things to be hidden.  If
    /// 'nohidden3d' is set for a particular plot element, that element
    /// will "shine through", even if hidden3d is set and the element is
    /// behind something else.
    Format&             nohidden3d (const bool nohidden3d = true);

    /// @brief Set the 'nocontours' property.
    /// @param nocontours Whether to disable contours.
    /// @return A reference to the Format.
    ///
    /// In Gnuplot 3D plots, you can choose globally to plot contours
    /// for each plot element.  If 'nocontours' is set for a particular
    /// plot element, contouring will be turned off for that element
    /// even if contours are active globally.
    Format&             nocontours (const bool nocontours = true);

    /// @brief Write to @c stream the Gnuplot 'with' code
    /// corresponding to the Axis.
    /// @param stream A reference to a std::ostream.
    /// @param format A reference to a Format.
    ///
    /// This method assembles the Gnuplot 'with' code corresponding to
    /// the sum of the Format settings, including the word "with"
    /// itself.
    friend std::ostream& operator<<(std::ostream& stream, const Format& format);

private:
    Style               style_;
    Value<int>          line_type_;
    Value<float>        line_width_;
    Value<ColorSpec>    line_color_;
    Value<int>          point_type_;
    Value<float>        point_size_;
    Value<std::string>  fill_style_;
    bool                nohidden3d_;
    bool                nocontours_;

};
    
class LineStyle {
public:
    LineStyle();
    LineStyle(int lw, ColorSpec lc, int pt);
    LineStyle& cycle(bool should_cycle=true);
    LineStyle& width(int lw);
    LineStyle& color(ColorSpec lc);
    LineStyle& point_type(int pt);
private:
    bool cycle_;
    Value<int> width_;
    Value<ColorSpec> color_;
    Value<int> point_type_;
    
    friend std::ostream& operator<<(std::ostream& stream, const LineStyle& style);
    friend class Gnuplot;
};


/// @brief Base class for plottable 2D things.
///
/// Plottable is a general abstract base class for any objects
/// (hereafter referred to as plot elements) that can be plotted in 2D
/// (with "plot" rather than "splot").  For most applications, you
/// will not need to work directly with Plottable.  For plotting
/// functions, see FunctionPlot; for plotting data, see DataPlot; and
/// for plotting vertical lines (an unfortunate special case in
/// Gnuplot), see VerticalLinePlot.
class Plottable {
public:

    /// @brief Construct a Plottable with a given format and title.
    /// @param title The title for this plot element.
    /// @param format A Format.
    Plottable (const std::string& title = "",
               const Format& format = Format::LINES);

    virtual ~Plottable ();

    /// @brief Get the Format object for this plot element.
    virtual Format&         format ();
    virtual const Format&   format () const;
    
    /// @brief Get or set the title for this plot element.
    std::string& title(){ return(title_); }
    const std::string& title() const{ return(title_); }
    std::string& title(const std::string& t){
        title_=t;
        return(title_);
    }

    /// @brief Return any inlined data as Gnuplot code in a string.
    ///
    /// This method should write any inlined data the plot element
    /// needs Gnuplot to see.  That is, anything this plot element needs
    /// to write on extra lines after the "plot" command should be
    /// written by this method.  @sa Try "help datafile special-filenames"
    /// in Gnuplot for more information.
    virtual std::string get_data_code () const;

    /// @brief Return any inlined data as Gnuplot code in a string.
    ///
    /// This method should write any inlined data the plot element
    /// needs Gnuplot to see.  That is, anything this plot element needs
    /// to write on extra lines after the "plot" command should be
    /// written by this method.  This method differs from the
    /// zero-argument overload in that it has access to the established
    /// x and y ranges. @sa Try "help datafile special-filenames"
    /// in Gnuplot for more information.
    virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const;

    /// @brief Return the plot code.
    ///
    /// This method should write anything this plot element needs
    /// after the word 'plot' and before any commas.  @sa Try "help
    /// datafile special-filenames" in Gnuplot for more information.
    virtual std::string get_plot_code () const = 0;


protected:
    Format          format_;
    std::string     title_;

};


/// @brief Plot a function.
///
/// FunctionPlot allows you to plot a function given in Gnuplot's
/// syntax as a string.
class FunctionPlot : public Plottable {
public:

    /// @brief Construct a FunctionPlot
    /// @param function The function.
    /// @param title The title of the function to be displayed in the legend
    /// @param format The format for displaying the function.
    FunctionPlot (const std::string& function,
                  const std::string& title = "",
                  const Format& format = Format::LINES);
    virtual ~FunctionPlot ();

    virtual std::string get_plot_code () const;

private:
    std::string     function_;

};


/// @brief Plot a vertical line.
///
/// VerticalLinePlot allows you to plot a vertical line at a specified
/// x value.
class VerticalLinePlot : public Plottable {
public:

    /// @brief Construct a VerticalLinePlot
    /// @param x The location of the vertical line.
    /// @param title The title of the function to be displayed in the legend
    /// @param format The format for displaying the vertical line.
    ///
    /// Internally, VerticalLinePlot works by plotting the line's
    /// endpoints.  Therefore, the user should to choose an
    /// appropriate Format (i.e. with style == Format::LINES or
    /// Format::LINES_POINTS).
    VerticalLinePlot (const double& x,
                      const std::string& title = "",
                      const Format& format = Format::LINES);
    virtual ~VerticalLinePlot ();

    virtual std::string get_plot_code () const;
    virtual std::string get_data_code () const;
    virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const;

private:
    double  x_;

};


/// @brief Plot 2D data.
///
/// DataPlot is a general purpose class for plotting 2D data.  It
/// tries to anticipate common data structures in which your data might
/// be contained.
///
/// To use DataPlot, you must provide two iterators that each
/// dereference to either a std::pair; a numerical array of length at
/// least 2; or any type with public numerical @c x() and @c y() methods.
/// The first iterator is the beginning of the data range and the
/// second iterator is just past the end of the data range.  In any of
/// these cases, the first member will be interpreted as an 'x'
/// coordinate and the second member will be interpreted as a 'y'
/// coordinate.  Here is an example of DataPlot being used in three
/// different ways:
/// @code
/// vector<shared_ptr<Plottable> > elements;
///
/// // data plot with iterator to pair
/// typedef vector<pair<float, float> > points_list_type;
/// points_list_type points;
/// for (float x = -10; x <= 10; x += 1) {
///     points.push_back (make_pair (x, x * x * x));
/// }
/// typedef DataPlot<points_list_type::const_iterator> data_plot_type;
/// shared_ptr<data_plot_type> dpv = make_shared<data_plot_type> (
///     points.begin (), points.end (),
///     "Cubic from vector", Format::LINES_POINTS);
/// elements.push_back (dpv);
///
/// 
/// // data plot with iterator to array of length two or more
/// const int COUNT = 21;
/// const int N_POINTS = 3;
/// float array_points[COUNT][N_POINTS];
/// for (int i = 0; i < COUNT; ++i) {
///     array_points[i][0] = i - COUNT / 2;
///     array_points[i][1] = 0.66 * (array_points[i][0]
///                                 * array_points[i][0]
///                                 * array_points[i][0]);
/// }
/// typedef DataPlot<const float (*)[N_POINTS]> array_data_plot_type;
/// shared_ptr<array_data_plot_type> dpa = make_shared<array_data_plot_type> (
///     array_points, array_points + COUNT,
///     "Cubic from array", Format::LINES_POINTS);
/// elements.push_back (dpa);
///
/// // data plot with iterator to a type with members x and y
/// MyPoint struct_points[COUNT];
/// for (int i = 0; i < COUNT; ++i) {
///     struct_points[i].x () = i - COUNT / 2;
///     struct_points[i].y () = 0.33 * (struct_points[i].x ()
///                                 * struct_points[i].x ()
///                                 * struct_points[i].x ());
/// }
/// typedef DataPlot<const MyPoint*> struct_data_plot_type;
/// shared_ptr<struct_data_plot_type> dps = make_shared<struct_data_plot_type> (
///     struct_points, struct_points + COUNT,
///     "Cubic from struct array", Format::LINES_POINTS);
/// elements.push_back (dps);
/// 
/// Gnuplot g;
/// g.plot_wxt (elements.begin (), elements.end ());
/// @endcode
/*template <typename ForwardIterator>
class DataPlot : public Plottable {
public:

    /// @brief Construct a DataPlot for data in the range (first, last).
    /// @param first An iterator to the first data point.
    /// @param last An iterator just past the last data point.
    /// @param title The title of the function to be displayed in the legend
    /// @param format The format for displaying the data points.
    ///
    /// See the overview of DataPlot for more information on @c first
    /// and @c last.
    DataPlot (const ForwardIterator& first,
              const ForwardIterator& last,
              const std::string& title = "",
              const Format& format = Format::POINTS);
    virtual ~DataPlot ();

    virtual std::string get_plot_code () const;
    virtual std::string get_data_code () const;
    virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const;
private:

    template <typename T1, typename T2>
        std::string do_get_data_code (const std::pair<T1, T2>*) const;
    template <typename T, int count>
        std::string do_get_data_code (const T(*const)[count]) const;
    template <typename PointT>
        std::string do_get_data_code (const PointT*) const;

    const ForwardIterator first_;
    const ForwardIterator last_;
};*/
    
class DataPlot : public Plottable {
public:
    struct point{
        double x, y;
        
        point(const point& p):x(p.x),y(p.y){}
        point(double x,double y):x(x),y(y){}
        template<typename T1, typename T2>
        point(const std::pair<T1,T2>& p):x(p.first),y(p.second){}
        template<typename T, unsigned int size>
        point(const T(*const a)[size]):x(a[0]),y(a[1]){
            static_assert(size==2,"Points can only be constructed from arrays of size 2");
        }
        //other overloads?
    };
    
    //Construct empty; use add_point
    DataPlot(const std::string& title = "",
             const Format& format = Format::POINTS):
    Plottable (title, format){}
    
    //Contruct from an iterator range
    template<typename ForwardIterator>
    DataPlot(ForwardIterator first, ForwardIterator last,
             const std::string& title = "",
             const Format& format = Format::POINTS):
    Plottable (title, format){
        for(; first!=last; ++first)
            data.emplace_back(point{*first});
    }
    //construct from an iterable container
    template<typename Container>
    DataPlot(const Container& c,
             const std::string& title = "",
             const Format& format = Format::POINTS):
    Plottable (title, format){
        for(const auto& item : c)
            data.emplace_back(point{item});
    }
    
    virtual ~DataPlot ();
    
    virtual std::string get_plot_code () const;
    virtual std::string get_data_code () const;
    virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const;
    
    template<typename T>
    void add_data_point(const T& t){
        data.push_back(point{t});
    }
private:
    std::vector<point> data;
};


/// @brief Base class for plottable 3D things.
///
/// Plottable3D is a general abstract base class for plot elements
/// that can be plotted in 3D (with "splot" rather than "plot").  
class Plottable3D : public Plottable{
public:

    /// @brief Construct a Plottable with a given format and title.
    /// @param title The title for this plot element.
    /// @param format A Format.
    Plottable3D (const std::string& title = "",
                 const Format& format = Format::LINES);

    virtual ~Plottable3D ();

    /// @brief Return any inlined data as Gnuplot code in a string.
    ///
    /// This method should write any inlined data the plot element
    /// needs Gnuplot to see.  That is, anything this plot element needs
    /// to write on extra lines after the "plot" command should be
    /// written by this method.  @sa Try "help datafile special-filenames"
    /// in Gnuplot for more information.
    virtual std::string get_data_code () const;

    /// @brief Return the plot code.
    ///
    /// This method should write anything this plot element needs
    /// after the word 'plot' and before any commas.  @sa Try "help
    /// datafile special-filenames" in Gnuplot for more information.
    virtual std::string get_plot_code () const = 0;

};

/// @brief Plot a function in 3D.
///
/// FunctionPlot3D allows you to plot a function in 3D given in Gnuplot's
/// syntax as a string.
class FunctionPlot3D : public Plottable3D {
public:

    /// @brief Construct a FunctionPlot3D
    /// @param function The function.
    /// @param title The title of the function to be displayed in the legend
    /// @param format The format for displaying the function.
    FunctionPlot3D (const std::string& function,
                    const std::string& title = "",
                    const Format& format = Format::LINES);
    virtual ~FunctionPlot3D ();

    virtual std::string get_plot_code () const;

private:
    std::string     function_;

};

/// @brief Plot 3D data.
///
/// DataPlot3D is a general purpose class for plotting 3D data.  To use
/// DataPlot3D, you must provide two iterators that each dereference to
/// either a boost::tuple with 3 members; a numerical array of length
/// at least 3; or any type with public numerical @c x, @c y, and @c z
/// members.  The first iterator is the beginning of the data range and
/// the second iterator is just past the end of the data range.  In any
/// of these cases, the first member will be interpreted as an 'x'
/// coordinate, the second member will be interpreted as a 'y'
/// coordinate and the third member will be interpreted as a 'z'
/// coordinate.
template <typename ForwardIterator>
class DataPlot3D : public Plottable3D {
public:

    /// @brief Construct a DataPlot3D for data in the range (first, last).
    /// @param first An iterator to the first data point.
    /// @param last An iterator just past the last data point.
    /// @param title The title of the function to be displayed in the legend
    /// @param format The format for displaying the data points.
    ///
    /// See the overview of DataPlot3D for more information on @c first
    /// and @c last.
    DataPlot3D (const ForwardIterator& first,
                const ForwardIterator& last,
                const std::string& title = "",
                const Format& format = Format::POINTS);
    virtual ~DataPlot3D ();

    virtual std::string get_plot_code () const;
    virtual std::string get_data_code () const;

protected:

    template <typename T1, typename T2, typename T3>
        boost::tuple<double, double, double>
        extract_point (const boost::tuple<T1, T2, T3>& p) const;
    template <typename T, int count>
        boost::tuple<double, double, double>
        extract_point (const T p[count]) const;
    template <typename PointT>
        boost::tuple<double, double, double>
        extract_point (const PointT& p) const;

    const ForwardIterator first_;
    const ForwardIterator last_;
};

/// @brief Plot 3D data as a mesh.
///
/// MeshDataPlot3D is a class for plotting 3D data as vertical bars.
/// To use MeshDataPlot3D, you must provide two iterators that each
/// dereference to either a boost::tuple with 3 members; a numerical
/// array of length at least 3; or any type with public numerical @c x(),
/// @c y(), and @c z() methods.  The first iterator is the beginning of the
/// data range and the second iterator is just past the end of the data
/// range.  In any of these cases, the first member will be interpreted
/// as an 'x' coordinate, the second member will be interpreted as a
/// 'y' coordinate and the third member will be interpreted as a 'z'
/// coordinate.  You also must provide the mesh dimensions.
template <typename ForwardIterator>
class MeshDataPlot3D : public DataPlot3D<ForwardIterator> {
public:

    /// @brief Construct a MeshDataPlot3D for data in the range (first, last).
    /// @param first An iterator to the first data point.
    /// @param last An iterator just past the last data point.
    /// @param n_x The width of the mesh (number of "rows")
    /// @param n_y The depth of the mesh (number of "columns")
    /// @param title The title of the function to be displayed in the legend
    /// @param format The format for displaying the data points.
    ///
    /// See the overview of DataPlot3D for more information on @c first
    /// and @c last.
    MeshDataPlot3D (const ForwardIterator& first,
                    const ForwardIterator& last,
                    const int n_x,
                    const int n_y,
                    const std::string& title = "",
                    const Format& format = Format::POINTS);
    virtual ~MeshDataPlot3D ();

    virtual std::string get_data_code () const;

protected:
    const int   n_x_;
    const int   n_y_;
};
    
class Terminal{
protected:
    std::string name_;
public:
    Terminal(std::string name):name_(name){}
    virtual ~Terminal(){}
    virtual void write_settings_code(std::ostream& os) const{}
    friend std::ostream& operator<<(std::ostream&, const Terminal&);
};
    
class PostscriptTerminal : public Terminal{
public:
    enum Mode{landscape,portrait,eps};
private:
    Mode mode_;
    bool enhanced_;
    unsigned int plexity_; //0=default, 1=simplex, 2=duplex
    bool color_;
    //TODO: background
    bool solid_;
    float dash_length_scale_;
    float line_width_scale_;
    //TODO: rounded
    //TODO: clip
    //TODO: palfuncparam (!)
    float x_size_; //inches
    float y_size_; //inches
    bool black_text_;
    std::string font_;
    unsigned int font_size_;
public:
    PostscriptTerminal();
    virtual void write_settings_code(std::ostream& os) const;
    
    PostscriptTerminal& mode(Mode m);
    PostscriptTerminal& enhanced(bool e);
    PostscriptTerminal& plexity(unsigned int plexity);
    PostscriptTerminal& color(bool c);
    PostscriptTerminal& solid(bool s);
    PostscriptTerminal& dash_length(float len);
    PostscriptTerminal& line_width(float width);
    PostscriptTerminal& x_size(float x);
    PostscriptTerminal& y_size(float y);
    PostscriptTerminal& black_text(bool b);
    PostscriptTerminal& font(std::string f);
    PostscriptTerminal& font_size(unsigned int s);
};

class Label{
public:
    ///Possible values for label text alignment
    enum Alignment{LEFT,CENTER,RIGHT};
    
    ///Construct a label
    /// @param text The text of the label.
    /// @param position The location of the label.
    explicit Label(std::string text="", Coordinate position=Coordinate());
    
    /// @brief Set the label text.
    /// @param t The label text.
    /// @return A reference to the Label
    Label& text(std::string t);
    /// @brief Get the label text.
    /// @return The text of the label.
    std::string text() const;
    
    /// @brief Set the label position.
    /// @param pos A gnuplot position (coordinates with optional coordinate
    ////           system specification).
    /// @return A reference to the Label.
    Label& at(Coordinate pos);
    /// @brief Set the label position.
    /// @param pos A gnuplot position (coordinates with optional coordinate
    ////           system specification).
    /// @return A reference to the Label.
    Label& position(Coordinate pos);
    /// @brief Get the label position.
    /// @return The position of the label.
    Coordinate at() const;
    /// @brief Get the label position.
    /// @return The position of the label.
    Coordinate position() const;
    
    /// @brief Set the label text lignment.
    /// @param a The alignment to use.
    /// @return A reference to the Label.
    Label& align(Alignment a);
    /// @brief Get the label text alignment.
    /// @return The label text alignment.
    Value<Alignment> align() const;
    
    /// @brief Set the label rotation angle.
    /// @param angle The rotation angle (in degrees).
    /// @return A reference to the Label.
    Label& rotate(double angle);
    /// @brief Get the label rotation angle.
    /// @return The rotation angle (in degrees).
    double rotate() const;
    
    /// @brief Set the label to be displayed in front of plotted data.
    Label& front();
    /// @brief Set the label to be displayed behind plotted data.
    Label& back();
    
    /// @brief Set the label font.
    /// @param name The name of the font to use.
    /// @param size The (optional) font size.
    /// @return A reference to the Label.
    Label& font(std::string name, Value<unsigned int> size=Value<unsigned int>());
    /// @brief Get the label font.
    /// @return The name of the label font.
    Value<std::string> font() const;
    
    /// @brief Set the label font size.
    /// @param size The font size.
    /// @return A reference to the Label.
    Label& font_size(unsigned int size);
    /// @brief Get the label font size.
    /// @return The label font size.
    Value<unsigned int> font_size() const;
    
    /// @brief Set the label text color.
    /// @param color A gnuplot colorspec.
    /// @return A reference to the Label.
    Label& font_color(ColorSpec color);
    /// @brief Set the label text color.
    /// @param color A gnuplot colorspec.
    /// @return A reference to the Label.
    Label& text_color(ColorSpec color);
    /// @brief Get the label text color.
    /// @return The label text color.
    Value<ColorSpec> font_color() const;
    /// @brief Get the label text color.
    /// @return The label text color.
    Value<ColorSpec> text_color() const;
    
private:
    std::string text_;
    Coordinate position_;
    Value<Alignment> text_align_;
    Value<double> rotation_;
    bool front_;
    Value<std::string> font_name_;
    Value<unsigned int> font_size_;
    Value<ColorSpec> font_color_;
    //TODO: point type
    //TODO: offset
    
    friend std::ostream& operator<<(std::ostream& stream, const Label& label);
};
    
class Arrow{
public:
	enum HeadType{nohead=0b00, head=0b01, backhead=0b10, heads=0b11};
	enum HeadFill{filled, empty, nofilled};
	
	Arrow():ht_(head),hf_(filled),front_(true){}
    
	Arrow(Coordinate from, Coordinate to, HeadType ht=head, HeadFill hf=filled):
	from_(from),to_(to),ht_(ht),hf_(hf),front_(true){}
	
	Arrow& from(Coordinate f){ from_=f; return(*this); }
	Coordinate& from(){ return(from_); }
	Arrow& to(Coordinate t){ to_=t; return(*this); }
	Coordinate& to(){ return(to_); }
	Arrow& width(double w){ width_=w; return(*this); }
	Value<double>& width(){ return(width_); }
	Arrow& fill(HeadFill f){ hf_=f; return(*this); }
	HeadFill& fill(){ return(hf_); }
	Arrow& front(bool f){ front_=f; return(*this); }
	bool& front(){ return(front_); }
	Arrow& head_length(double hl){ head_length_=hl; return(*this); }
	Value<double>& head_length(){ return(head_length_); }
	Arrow& head_angle(double ha){ head_angle_=ha; return(*this); }
	Value<double>& head_angle(){ return(head_angle_); }
	Arrow& head_back_angle(double hba){ head_back_angle_=hba; return(*this); }
	Value<double>& head_back_angle(){ return(head_back_angle_); }
	
private:
	Coordinate from_, to_;
	HeadType ht_;
	HeadFill hf_;
	Value<double> width_;
	bool front_;
	//TODO: size: length, angle, backangle,
	Value<double> head_length_, head_angle_, head_back_angle_;
	
	friend std::ostream& operator<<(std::ostream&, const Arrow&);
};
	
class Key{
public:
	Key();
	
	Key& on();
	Key& off();
	Key& show(bool s=true);
    
    Key& left();
    Key& right();
    Key& center_h();
    Key& top();
    Key& bottom();
    Key& center_v();
    Key& lmargin();
    Key& rmargin();
    Key& tmargin();
    Key& bmargin();
    Key& above();
    Key& over();
    Key& below();
    Key& under();
    Key& at(Coordinate pos);
    //inside/outside?
    
    
    Key& vertical();
    Key& horizontal();
    
    Key& left_justify();
    Key& right_justify();
    
    Key& opaque(bool o=true);
    Key& reverse(bool r=true);
    Key& invert(bool i=true);
    
    Key& sample_length(double len);
    Key& spacing(double space);
    Key& width(double w);
    Key& height(double h);
    
    //TODO: box settings
    
    Key& maxcols(unsigned int max);
    Key& maxrows(unsigned int max);
    
private:
    
	bool show_;
	bool inside_;
};
    
template<typename BaseType, typename MaybeDerivedType>
using IfIsBaseOf=std::enable_if<std::is_base_of<BaseType,typename std::remove_reference<MaybeDerivedType>::type>::value,void>;

namespace keywords{
    
    struct output_file{
        std::string filename;
        output_file(const std::string& filename):filename(filename){}
    };
    struct grid{
        bool show;
        grid(bool show):show(show){}
    };
    struct table{
        bool use_table;
        table(bool use_table):use_table(use_table){}
    };
    struct surface{
        bool show;
        surface(bool show):show(show){}
    };
    struct line_style{
        unsigned int index;
        LineStyle style;
        line_style(unsigned int index, LineStyle style):index(index),style(style){
            if(index<1)
                throw std::runtime_error("Line styel indices are one based");
        }
    };
    struct bar_size{
        double size;
        bar_size(double size):size(size){}
    };
    struct bars_front{
        bool front;
        bars_front(bool front):front(front){}
    };
    struct contour{
        bool base, surface;
        contour(bool base, bool surface):base(base),surface(surface){}
    };
    struct contour_levels{
        std::vector<double> levels;
        template<typename DoubleContainer>
        contour_levels(const DoubleContainer& levels){
            std::copy(std::begin(levels),std::end(levels),std::back_inserter(this->levels));
        }
        template<typename DoubleIterator>
        contour_levels(DoubleIterator levels_begin, DoubleIterator levels_end){
            std::copy(levels_begin,levels_end,std::back_inserter(levels));
        }
    };
    struct manual_settings{
        std::string settings;
        manual_settings(const std::string& settings):settings(settings){}
    };
    
}

/// @brief Plot things with Gnuplot.
///
/// This class provides an interface to the external Gnuplot binary
/// (assumed by default to be called "gnuplot").  It has accessors for
/// adjusting global settings and other methods for plotting Plottable
/// and Plottable3D objects.
class Gnuplot {
public:
    /// @brief Construct a Gnuplot object.
    /// @param path_to_gnuplot The name of or path to the gnuplot binary
    ///
    /// This constructor sets up a Gnuplot object, opening a Pipe to
    /// and from the specified gnuplot binary (or just "gnuplot" by
    /// default).
    Gnuplot (const /*boost::filesystem::path*/std::string& path_to_gnuplot = "gnuplot");
    ~Gnuplot ();

    /// @brief Get the plot title.
    /// @return The title.
    std::string     title ();

    /// @brief Set the plot title.
    /// @param title The new title.
    void            title (const std::string& title);

    /// @brief Get the x Axis object.
    /// @return A reference to the x Axis object.
    ///
    /// This method is intended to be combined with Axis's method
    /// chaining.  For example:
    /// @code
    /// Gnuplot g;
    /// g.x_axis ().min (0).label ("t (ns)");
    /// @endcode
    Axis&           x_axis ();

    /// @brief Like Gnuplot::x_axis(), but for the y Axis.
    /// @return A reference to the y Axis object.
    Axis&           y_axis ();

    /// @brief Like Gnuplot::x_axis(), but for the z Axis.
    /// @return A reference to the z Axis object.
    Axis&           z_axis ();
    
    /// @brief Like Gnuplot::x_axis(), but for the color Axis.
    /// @return A reference to the color Axis object.
    Axis&           cb_axis ();
    
    /// @brief Get a line style by index, allowing it to be altered.
    /// @param num The index of the line style to fetch, which must be >0.
    /// @return The requested global line style object.
    LineStyle&      line_style(unsigned int num);
    
    /// @brief Immediately unset a previously set line style.
    /// @param num The index of the line style to unset, which must be >0.
    void            unset_line_style(unsigned int num);
    
    /// @brief Get or create a label by index.
    /// @param num The index of the label, which must be >0.
    /// @return The requested label.
    Label&          label(unsigned int num);
    
    /// @brief Remove a label.
    /// @param num The index of the label to remove, which must be >0.
    void            unset_label(unsigned int num);
    
    /// @brief Add a label.
    /// @param l The label to add.
    /// @return The index at which the newly added label was stored.
    unsigned int    add_label(const Label& l);
    
    /// @brief Add a label.
    /// @param l The label to add.
    /// @return The index at which the newly added label was stored.
    unsigned int    add_label(Label&& l);
    
    /// @brief Get or create an arrow by index.
    /// @param num The index of the arrow, which must be >0.
    /// @return The requested arrow.
    Arrow&          arrow(unsigned int num);
    
    /// @brief Remove an arrow.
    /// @param num The index of the arrow to remove, which must be >0.
    void            unset_arrow(unsigned int num);
    
    /// @brief Add an arrow.
    /// @param a The arrow to add.
    /// @return The index at which the newly added arrow was stored.
    unsigned int    add_arrow(const Arrow& l);
    
    /// @brief Add an arrow.
    /// @param a The arrow to add.
    /// @return The index at which the newly added arow was stored.
    unsigned int    add_arrow(Arrow&& l);
    
    /// @brief Set the size of tics on the ends of error bars.
    /// @param size The size to set. If negative this will be translated to
    ///             "fullwidth" which makes the widths of error bar tics the
    ///             same as the width of the item with which the error bar is
    ///             associated.
    void            set_bar_size(double size);
    
    /// @brief Set whether error bars are drawn in front.
    /// @param front Whether error bars should be drawn in front (true) or
    ///              in back (false).
    void            set_bars_front(bool front);
    
    /// @brief Set whether the 3D surfaces are shown.
    /// @param show Whether surfaces should be shown.
    void show_surface(bool show=true);
    
    /// @brief Set whether the coordinate grid overlay is shown.
    /// @param show Whether the grid should be shown.
    /// @todo Replace this with a system which allows setting all grid properties
    void show_grid(bool show=true);
    
    /// @brief Activate 'table' mode.
    ///
    /// In 'table' mode which plot data is written as a table to the specified output
    /// location instead of being drawn normally.
    void set_table(bool table_output=true){ table_output_=table_output; }
    
    /// @brief Deactivate 'table' mode.
    void unset_table(){ table_output_=false; }
    
    /// @brief Control drawing of contours.
    ///
    /// @param base Whether to draw contours at the base of the plot.
    /// @param surface Whether to draw contours on the plotted surfaces.
    void set_contour(bool base, bool surface);
    
    template<typename DoubleContainer>
    void set_contour_levels(DoubleContainer levels){
        contour_levels_.clear();
        std::copy(std::begin(levels),std::end(levels),std::back_inserter(contour_levels_));
    }
    
    void manual_settings(const std::string& settings);
    
    //EXPERIMENTAL
public:
    void set_output(const std::string& file){
        output_ = file;
    }
    
    /*template <typename TermType=Terminal, typename... ArgTypes>
    void set_terminal(ArgTypes&&... args){
        terminal_ = std::unique_ptr<Terminal>(new TermType(std::forward<ArgTypes>(args)...));
        output_.clear();
    }*/
    
    void set_terminal(std::unique_ptr<Terminal>&& terminal){
        terminal_ = std::move(terminal);
        output_.clear();
    }
    
    template <typename TermType>
    void set_terminal(TermType&& terminal){
        terminal_.reset(new typename std::remove_reference<TermType>::type(std::forward<TermType>(terminal)));
        output_.clear();
    }
    
    template <typename TermType>
    void set_terminal(const TermType& terminal){
        terminal_.reset(new typename std::remove_reference<TermType>::type(terminal));
        output_.clear();
    }
    
    
private:
    template<unsigned int Index, typename PlotType, typename... Others>
    void checkPlotArguments(const PlotType& plot, const Others&... others){
        static_assert(std::is_base_of<Plottable,typename std::remove_reference<PlotType>::type>::value
                      || std::is_base_of<Plottable3D,typename std::remove_reference<PlotType>::type>::value,
                      "All arguments to Gnuplot::plot() must be Plottable objects, "
                      "except the last, which may be a string");
        checkPlotArguments<Index+1>(others...);
    }
    
    template<unsigned int Index>
    void checkPlotArguments(){
        static_assert(Index>0,"Gnuplot::plot() requires at least one argument");
    }
    
    template<unsigned int Index>
    void checkPlotArguments(const std::string& output){
        static_assert(Index>0,"Gnuplot::plot() requires at least Plottable argument");
    }
    
    template<unsigned int Index>
    void checkPlotArguments(const char* output){
        static_assert(Index>0,"Gnuplot::plot() requires at least Plottable argument");
    }
    
    
    template<unsigned int Index, typename PlotType, typename... Others>
    void checkPlot3DArguments(const PlotType& plot, const Others&... others){
        static_assert(std::is_base_of<Plottable3D,typename std::remove_reference<PlotType>::type>::value,
                      "All arguments to Gnuplot::plot3d() must be Plottable3D objects, "
                      "except the last, which may be a string");
        checkPlot3DArguments<Index+1>(others...);
    }
    
    template<unsigned int Index>
    void checkPlot3DArguments(){
        static_assert(Index>0,"Gnuplot::plot() requires at least one argument");
    }
    
    template<unsigned int Index>
    void checkPlot3DArguments(const std::string& output){
        static_assert(Index>0,"Gnuplot::plot() requires at least Plottable argument");
    }
    
    template<unsigned int Index>
    void checkPlot3DArguments(const char* output){
        static_assert(Index>0,"Gnuplot::plot() requires at least Plottable argument");
    }
    
    
    
    //Base case with a file specified
    void maybe_set_output(const std::string& file){
        set_output(file);
    }
    void maybe_set_output(const char* file){
        set_output(file);
    }
    
    //Base case with no file specified
    void maybe_set_output(){ /* Do Nothing */}
    
    template<typename First, typename... Others>
    void maybe_set_output(const First&, const Others&... others){
        maybe_set_output(others...);
    }
	
    void do_stuff(){}
    
    template<typename PlotType, typename... OtherArgs>
	typename IfIsBaseOf<Plottable,PlotType>::type
    do_stuff(const PlotType& plot, OtherArgs&&... others){
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
	template<typename... OtherArgs>
	void do_stuff(const Axis& axis, OtherArgs&&... others){
		if(axis.name()=="x")
            x_axis_=axis;
        else if(axis.name()=="y")
            y_axis_=axis;
        else if(axis.name()=="z")
            z_axis_=axis;
        else if(axis.name()=="cb")
            cb_axis_=axis;
        else
            throw std::runtime_error("Unrecognized axis name: "+axis.name());
        do_stuff(std::forward<OtherArgs>(others)...);
	}
    
    template<typename TermType, typename... OtherArgs>
	typename IfIsBaseOf<Terminal,TermType>::type
    do_stuff(const TermType& term, OtherArgs&&... others){
        set_terminal(term);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename TermType, typename... OtherArgs>
    typename IfIsBaseOf<Terminal,TermType>::type
    do_stuff(TermType&& term, OtherArgs&&... others){
        set_terminal(std::move(term));
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(const Label& label, OtherArgs&&... others){
        add_label(label);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(Label&& label, OtherArgs&&... others){
        add_label(std::move(label));
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(const std::string& outputFileName, OtherArgs&&... others){
        set_output(outputFileName);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(const char* outputFileName, OtherArgs&&... others){
        set_output(outputFileName);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::output_file outputFileName, OtherArgs&&... others){
        set_output(outputFileName.filename);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::grid show, OtherArgs&&... others){
        show_grid(show.show);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::surface show, OtherArgs&&... others){
        show_surface(show.show);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::table table, OtherArgs&&... others){
        set_table(table.use_table);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::line_style style, OtherArgs&&... others){
        line_style(style.index)=style.style;
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::bar_size size, OtherArgs&&... others){
        set_bar_size(size.size);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::bars_front front, OtherArgs&&... others){
        set_bars_front(front.front);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::contour c, OtherArgs&&... others){
        set_contour(c.base,c.surface);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::contour_levels levels, OtherArgs&&... others){
        set_contour_levels(levels.levels);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    template<typename... OtherArgs>
	void do_stuff(keywords::manual_settings settings, OtherArgs&&... others){
        manual_settings(settings.settings);
        do_stuff(std::forward<OtherArgs>(others)...);
    }
    
    
    
    template<unsigned int Index, typename OpType, typename PlotType, typename... Others,
             typename enable=typename IfIsBaseOf<Plottable,PlotType>::type>
    void accumulate_to_stream(std::ostream& os, const std::string& separator,
                              OpType op, PlotType&& plot, Others&&... others){
        if(Index)
            os << separator;
        os << op(plot);
        accumulate_to_stream<Index+1>(os,separator,op,std::forward<Others>(others)...);
    }
    
//    template<unsigned int Index, typename OpType>
//    void accumulate_to_stream(std::ostream& os, const std::string& separator, OpType op, std::string last){
//        /* Do Nothing */
//    }
//    
//    template<unsigned int Index, typename OpType>
//    void accumulate_to_stream(std::ostream& os, const std::string& separator, OpType op, const char* last){
//        /* Do Nothing */
//    }
//    
//    template<unsigned int Index, typename OpType>
//    void accumulate_to_stream(std::ostream& os, const std::string& separator, OpType op){
//        /* Do Nothing */
//    }
    
    template<unsigned int Index, typename OpType, typename... Stuff>
    void accumulate_to_stream(std::ostream& os, const std::string& separator,
                              OpType op, const Stuff&... stuff){
        //stop
    }
    
public:
    ///\brief Plot an arbitrary number of Plottables using the current terminal
    ///The final argumenat may optionally be a string, in which case it will be
    ///interpreted as the output file to which to write the plot
    template<typename... ArgTypes>
    void plot(ArgTypes&&... args){
        //checkPlotArguments<0>(std::forward<ArgTypes>(args)...);
        do_stuff(args...);
        std::ostringstream gp_code_plot_line;
        gp_code_plot_line << "plot ";
        accumulate_to_stream<0>(gp_code_plot_line,", ",
                                [](const Plottable& plot){ return(plot.get_plot_code()); },
                                std::forward<ArgTypes>(args)...);
        gp_code_plot_line << std::endl;
        std::ostringstream gp_code_pre_norm;
        gp_code_pre_norm << "reset" << std::endl
        << get_settings_code() << gp_code_plot_line.str();
        accumulate_to_stream<0>(gp_code_pre_norm,"",
                                [](const Plottable& plot){ return(plot.get_data_code()); },
                                std::forward<ArgTypes>(args)...);
        // deal with normalization
        NormalizeResult norm = normalize_once (gp_code_pre_norm.str ());
        // do real plotting
        std::ostringstream gp_code_post_norm;
        gp_code_post_norm << "reset" << std::endl
        << get_settings_code() << norm << gp_code_plot_line.str();
        accumulate_to_stream<0>(gp_code_post_norm,"",
                                [&norm](const Plottable& plot){ return(plot.get_data_code(norm.x_min, norm.x_max, norm.y_min, norm.y_max)); },
                                std::forward<ArgTypes>(args)...);
        //std::cout << "Settings: " << gp_code_post_norm.str() << std::endl;
        commands << *terminal_;
        //maybe_set_output(args...);
        //do_stuff(args...);
        if(!output_.empty())
            commands << "set output '" << output_ << '\'' << std::endl;
        else
            commands << "set output" << std::endl;
        if(table_output_)
            commands << "set table" << std::endl;
        commands << gp_code_post_norm.str();
        if(table_output_)
            commands << "unset table" << std::endl;
    }
    
    template<typename IteratorType, typename = typename std::enable_if<!std::is_base_of<Plottable,IteratorType>::value>::type>
    void plotI(IteratorType begin, IteratorType end, std::string output=""){
        //TODO: type checking to ensure that the function is only called with pointer-like objects to Plottable
        std::ostringstream gp_code_plot_line;
        gp_code_plot_line << "plot ";
        for(IteratorType it=begin; it!=end; it++){
            if(it!=begin)
                gp_code_plot_line << ", ";
            gp_code_plot_line << (*it)->get_plot_code();
        }
        gp_code_plot_line << std::endl;
        std::ostringstream gp_code_pre_norm;
        gp_code_pre_norm << "reset" << std::endl
        << get_settings_code() << gp_code_plot_line.str();
        for(IteratorType it=begin; it!=end; it++)
            gp_code_pre_norm << (*it)->get_data_code();
		//std::cout << "Pre-norm Settings: " << gp_code_pre_norm.str() << std::endl;
        // deal with normalization
        NormalizeResult norm = normalize_once (gp_code_pre_norm.str ());
        // do real plotting
        std::ostringstream gp_code_post_norm;
        gp_code_post_norm << "reset" << std::endl
        << get_settings_code() << norm << gp_code_plot_line.str();
        for(IteratorType it=begin; it!=end; it++)
            gp_code_post_norm << (*it)->get_data_code(norm.x_min, norm.x_max, norm.y_min, norm.y_max);
        //std::cout << "Settings: " << gp_code_post_norm.str() << std::endl;
        commands << *terminal_;
        if(!output.empty())
            output_=output;
        if(!output_.empty())
            commands << "set output '" << output_ << '\'' << std::endl;
        if(table_output_)
            commands << "set table" << std::endl;
        commands << gp_code_post_norm.str();
        if(table_output_)
            commands << "unset table" << std::endl;
    }
    
    template<typename IterableType, typename = typename std::enable_if<!std::is_base_of<Plottable,IterableType>::value>::type>
    void plotI(const IterableType& container, std::string output=""){
        plotI(std::begin(container),std::end(container),output);
    }
    
    template<typename... ArgTypes>
    void plot3d(ArgTypes&&... args){
        checkPlot3DArguments<0>(std::forward<ArgTypes>(args)...);
        std::ostringstream gp_code_plot_line;
        gp_code_plot_line << "splot ";
        accumulate_to_stream<0>(gp_code_plot_line,", ",
                                [](const Plottable3D& plot){ return(plot.get_plot_code()); },
                                std::forward<ArgTypes>(args)...);
        gp_code_plot_line << std::endl;
        std::ostringstream gp_code;
        gp_code << "reset" << std::endl
            << get_settings_code() << gp_code_plot_line.str();
        accumulate_to_stream<0>(gp_code,"",
                                [](const Plottable3D& plot){ return(plot.get_data_code()); },
                                std::forward<ArgTypes>(args)...);
        commands << *terminal_;
        maybe_set_output(args...);
        if(!output_.empty())
            commands << "set output '" << output_ << '\'' << std::endl;
        if(table_output_)
            commands << "set table" << std::endl;
        commands << gp_code.str();
        if(table_output_)
            commands << "unset table" << std::endl;
    }
    
    //TODO: this won't play nicely with plot(), since it issues multiple 'plot' commands
    void begin_multiplot(unsigned int rows=0, unsigned int columns=0);
    void end_multiplot();

    /// \brief Allow the user to interact directly with the child gnuplot process.
    ///
    /// Sends user input from cin to gnuplot's stdin, and output from gnuplot's
    /// stderr to cout. The interactive subsession can be ended by typing "exit"
    /// or ".q"; in either case the child gnuplot process will continue running.
    /// This functionality depends on using predetermined marker strings written
    /// through gnuplot's `print` command, so it is important that the `set print`
    /// command NOT be used to redirect the output of `print` to anywwhere other
    /// then stderr.
    ///
    /// This faeture is still a work in progress. known bugs include:
    ///     *gnuplot output to stdout is not shown
    ///     *signals are not forwarded to gnuplot
    ///     *paging in gnuplot help is broken (it triggers an infinite loop)
    void interact();
    
    /// \brief Fetches any text written by gnuplot to stdout.
    ///
    /// This function is intended for use with set_output("") and set_table().
    std::string get_raw_output();
    
    std::string get_raw_errors();
    
    //END EXPERIMENTAL

private:

    std::string     get_settings_code () const;

    struct NormalizeResult {
        double x_min, x_max, y_min, y_max;
    };

    friend std::ostream& operator<< (std::ostream& stream,
                                     const NormalizeResult& normalization);

    NormalizeResult normalize_once (const std::string& gnuplot_code) const;

    Axis                    x_axis_;
    Axis                    y_axis_;
    Axis                    z_axis_;
    Axis                    cb_axis_;
    std::string             title_;
    std::map<unsigned int, LineStyle> line_styles_;
    std::map<unsigned int, Label> labels_;
    std::map<unsigned int, Arrow> arrows_;
    phys_tools::Value<std::string> error_bar_size;
    phys_tools::Value<std::string> error_bar_front;
    phys_tools::Value<bool> show_surface_;
    phys_tools::Value<bool> show_grid_;
    bool table_output_;
    bool base_contours_, surface_contours_;
    std::vector<double> contour_levels_;
    std::string manual_settings_;

    mutable std::ofstream   commands;
    mutable std::ifstream   outputs;
    mutable std::ifstream   errors;
    Pipe                    pipe;
    
    std::unique_ptr<Terminal> terminal_;
    std::string output_;

};

/// @}

} // namespace phys_tools::gnuplot
} // namespace phys_tools


#include "gnuplot.tcpp"
#endif  /* GNUPLOT_HPP */
