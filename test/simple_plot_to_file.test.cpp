#include <PhysTools/gnuplot.h>
#include <unistd.h>

int main(){
    using namespace phys_tools::gnuplot;
	
    const std::string filename="simple_plot_to_file.dat";
	{
		Gnuplot g;
		g.set_terminal(Terminal("dumb size 60,15"));
		g.x_axis().min(-10).max(10);
		g.y_axis().min(0).max(100);
		g.plot(FunctionPlot("x**2","quadratic"),filename);
	} //important: need to force gnuplot process to end
	// causing open files to be flushed to disk and closed!
	
	{
		std::ifstream plotfile(filename);
		if(!plotfile.good()){
			std::cout << "Failed to open plot data. Not written?" << std::endl;
			return(0);
		}
		std::string line;
		while(getline(plotfile,line))
			std::cout << line << '\n';
		unlink(filename.c_str());
	}
}