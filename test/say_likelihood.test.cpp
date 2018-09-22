#include <iostream>
#include <cmath>
#include <random>
#include <numeric>
#include <cassert>
#include <algorithm>

#include <PhysTools/likelihood/likelihood.h>

int main() {
    phys_tools::likelihood::SAYLikelihood cat;
    phys_tools::likelihood::poissonLikelihood p;
    int nEvents = 1e4;
    std::vector<double> weights(nEvents);
    std::vector<double> weightsSq(nEvents);

    double xMin = 1e-12, xMax=1e0;
	std::mt19937 rng(52768);

    std::uniform_real_distribution<double> dist(xMin,xMax);
    for(size_t i=0; i<nEvents; i++) {
        weights[i] = dist(rng);
        weightsSq[i] = dist(rng);
    }

    std::vector<double> w = {std::accumulate(weights.begin(), weights.end(), double(0), std::plus<double>())};
    std::vector<double> w2 = {std::accumulate(weightsSq.begin(), weightsSq.end(), double(0), std::plus<double>())};
    std::vector<double> results;
    double poisson_result = p(0,w,w2);


    for(double i = -7; i<10; i+=0.5) {
        for(double k=0; k<20; ++k) {
            assert(cat(k,w,w2) > p(k,w,w2));
            //std::cout << cat(k,w,w2) << " " << p(k,w,w2) << std::endl;
        }
        double cat_res = cat(4000,w,w2);
        std::cout << cat_res << " " << poisson_result << std::endl;
        results.push_back(cat_res);
        w2[0]*=pow(10,0.5);
    }
    w2[0] = w2[0] / std::pow(10, 10*2);
    std::cout << cat(0,w,w2) << std::endl;
    std::cout << p(0,w,w2) << std::endl;
    std::cout << cat(0,w,w2) - p(0,w,w2) << std::endl;
    std::cout << (cat(0,w,w2) - p(0,w,w2))/p(0,w,w2) << std::endl;
    std::cout << std::abs((cat(0,w,w2) - p(0,w,w2))/p(0,w,w2)) << std::endl;
    std::cout << w[0] << std::endl;
    assert(std::is_sorted(results.begin(), results.end()));
}
