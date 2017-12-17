#ifndef __RESIDUALS_H__
#define __RESIDUALS_H__

#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <deque>

#include <boost/math/constants/constants.hpp>

#include <boost/mpl/if.hpp>

#include "./autodiff.h"

namespace phys_tools {
namespace residuals {

template<typename T>
struct residual_computer {
    std::vector<T> operator()(std::vector<T>const& zj, std::vector<T>const& nj, std::vector<T>const& sj, std::vector<T>const& mj) {
        unsigned int N = nj.size();
        unsigned int M = mj.size();
        std::vector<T> ckmk();
        for(unsigned int k=0; k<M; ++k) {
            for(unsigned int j=0; j<N; ++j) {
            
            }
        }
    }
    
}

} // namespace residuals
} // namespace phys_tools

#endif
