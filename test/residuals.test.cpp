#include <iostream>
#include "../PhysTools/residuals.h"

int main() {
    std::vector<double> z = {-3., 1., 2., 3.};
    std::vector<unsigned int> n {1, 3, 3, 1};
    std::vector<double> s = {0., -1., -2.};
    std::vector<unsigned int> m = {1, 4, 1};
    std::vector<double> c = phys_tools::residuals::residual_computer<double>()(z, n, s, m);
    for(unsigned int i=0; i<c.size(); ++i) {
        std::cout << c[i] << std::endl;
    }
    return 0;
}
