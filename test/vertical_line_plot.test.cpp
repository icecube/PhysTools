#include <PhysTools/gnuplot.h>

int main(){
    using namespace phys_tools::gnuplot;
    
    Gnuplot g;
    g.set_terminal(Terminal("dumb size 60,15"));
    g.x_axis().min(-10).max(10);
    g.y_axis().min(0).auto_max();
    g.plot(FunctionPlot("x**2","quadratic"),VerticalLinePlot(-5));
    std::cout << g.get_raw_output() << std::endl;
}