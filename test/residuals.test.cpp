#include <iostream>
#include "../PhysTools/residuals.h"

int main() {
    std::vector<double> z = {0};
    std::vector<unsigned int> n {5};
    std::vector<double> s = {0.0001, 0.0001001, 0.0001002, 0.0001003};
    std::vector<unsigned int> m = {1, 1, 1, 1};
    std::vector<double> c = phys_tools::residuals::residual_computer<double>()(z, n, s, m);
    double c_sum = 0;
    for(unsigned int i=0; i<c.size(); ++i) {
        std::cout << c[i] << std::endl;
        c_sum += c[i];
    }
    std::cout << c_sum << std::endl;
    return 0;
}
