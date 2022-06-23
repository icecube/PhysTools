#include <PhysTools/likelihood/likelihood.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <random>


bool double_compare(double x, double y){
    double diff = std::abs(x-y);
    double eps=(1e-6);
    return (diff<eps);
}

int main(int argc, char** argv){
    /*
    Checking the 2D case against the known working Gaussian 2D prior
    */


    using namespace phys_tools::likelihood;

    // these are chosen to mimic the Snowstorm effective gradients used in MEOWS
    double alpha = 0.05091035738186185100;
    std::vector<std::vector<double>> corr{{1, alpha},
                                          {alpha, 1}};
    std::vector<double> means{0.0, 0.0};
    std::vector<double> sigmas{1.0, 1.0};

    const size_t n_vals = static_cast<size_t>(2);

    GaussianNDPrior<2> testprio(means, sigmas, corr);
    Gaussian2DPrior otherprio(means[0], means[1], sigmas[0], sigmas[1], alpha);

    std::vector<double> test_values = {0,0.5,1.0,1.5,2.0,2.5};

    for (unsigned int i=0; i<test_values.size(); i++){
        for (unsigned int j=0; j<test_values.size(); j++){
        
        std::vector<double> this_test = {test_values[i], test_values[j]};
            auto result = testprio(this_test);
            auto other = otherprio(test_values[i], test_values[j]);
            bool are_equal = double_compare(result, other);
            if (!are_equal){
                throw std::runtime_error("Test failed!");
            }

        }
    }
    return 0;

}