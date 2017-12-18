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

namespace phys_tools {
namespace residuals {

template<typename T>
struct residual_computer {
    std::vector<T> operator()(std::vector<T>const& z, std::vector<unsigned int>const& n, std::vector<T>const& s, std::vector<unsigned int>const& m) {
        // Useful things to precompute
        const unsigned int N = n.size();
        const unsigned int M = m.size();
        const unsigned int max_n = *std::max_element(n.begin(), n.end());
        const unsigned int max_m = *std::max_element(m.begin(), m.end());

        std::vector<T> ss_diff(M*M);
        std::vector<T> zs_diff(N*M);
        std::vector<T> c(M*(max_m), T(0));
        #define ss_diff(r, c) (ss_diff[(r)*M + (c)])
        #define zs_diff(r, c) (zs_diff[(r)*M + (c)])
        #define c(i, j) (c[(i)*max_m + (j)])

		for(unsigned int j=0; j<M; ++j) {
	    	for(unsigned int i=0; i<M; ++i) {
                ss_diff(i, j) = s[i] - s[j];
            }
	    	for(unsigned int i=0; i<N; ++i) {
                zs_diff(i, j) = z[i] - s[j];
            }
		}

        // The algorithm
        for(unsigned int k=0; k<M; ++k) {
            T numerator = T(1);
            T denominator = T(1);
            for(unsigned int j=0; j<N; ++j) {
                //numerator *= pow(s[k]-z[j], n[j]);
                numerator *= pow(-zs_diff(j, k), n[j]);
            }
            for(unsigned int j=0; j<M; ++j) {
                if(j != k)
                    //denominator *= pow(s[k]-s[j], m[j]);
                    denominator *= pow(ss_diff(k, j), m[j]);
            }
            c(k, m[k]-1) = numerator/denominator;
        }
        for(unsigned int k=0; k<M; ++k) {
            std::vector<T> lambda(m[k]-1);
            for(unsigned int i=0; i<m[k]-1; ++i) {
                lambda[i] = T(0);
                for(unsigned int j=0; j<M; ++j) {
                    if(j != k)
                        lambda[i] += m[j]/pow(ss_diff(j, k), i+1);
                }
                for(unsigned int j=0; j<N; ++j) {
                    lambda[i] -= n[j]/pow(zs_diff(j, k), i+1);
                }
            }
            for(unsigned int L=m[k]-1; L>=1; --L) {
                //c(k, L) = 0;
                for(unsigned int i=1; i<=m[k] - L; ++i) {
                    c(k, L-1) += c(k, L+i-1) * lambda[i-1];
                }
                c(k, L-1) = c(k, L-1) / (m[k] - L);
            }
        }

        // Clean up the nasty stuff we defined
        #undef ss_diff
        #undef zs_diff
		#undef c
        return c;
    }
    
};

template<typename T>
struct contour_integral {
    T operator()(std::vector<T>const& z, std::vector<unsigned int>const& n, std::vector<T>const& s, std::vector<unsigned int>const& m) {
        const unsigned int M = m.size();
        const unsigned int max_m = *std::max_element(m.begin(), m.end());
        std::vector<T> c = residual_computer<T>()(z, n, s, m);
        #define c(i, j) (c[(i)*max_m + (j)])

        T residual_sum(0);

        std::cout << "residuals:" << std::endl;
        for(unsigned int k=0; k<M; ++k) {
            std::cout << "c(" << k << ",0) = " << c(k,0) << std::endl;
            residual_sum += c(k,0);
        }
        #undef c
        return residual_sum;
    }
};

} // namespace residuals
} // namespace phys_tools

#endif
