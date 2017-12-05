#include <PhysTools/gnuplot.h>

int main(){
    using namespace phys_tools::gnuplot;
    
    Gnuplot g;
    g.set_terminal(Terminal("dumb size 60,15"));
    g.x_axis().min(-10).max(10);
    g.y_axis().min(0).max(100);
    g.set_table();
    g.plot(FunctionPlot("x**2","quadratic"),"");
    std::cout << g.get_raw_output() << std::endl;
    g.unset_table();
    g.plot(FunctionPlot("x**2","quadratic"));
    std::cout << g.get_raw_output() << std::endl;
}