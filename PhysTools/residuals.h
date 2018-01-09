#ifndef __RESIDUALS_H__
#define __RESIDUALS_H__

#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <deque>
#include <cfenv>
#include <ttmath/ttmath.h>

#include <boost/math/constants/constants.hpp>
#include <boost/range/adaptor/strided.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <boost/mpl/if.hpp>

#include "autodiff.h"

namespace phys_tools {
namespace residuals {
template <class InIt>
typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
    typedef typename std::iterator_traits<InIt>::value_type real;
    real sum = real();
    real running_error = real();
    real temp;
    real difference;

    for (; begin != end; ++begin) {
        difference = *begin;
        difference -= running_error;
        temp = sum;
        temp += difference;
        running_error = temp;
        running_error -= sum;
        running_error -= difference;
        sum = std::move(temp);
    }
    return sum;
}

template<typename T>
struct almost_equal {
    typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
        operator()(T x, T y, int ulp)
    {
        // the machine epsilon has to be scaled to the magnitude of the values used
        // and multiplied by the desired precision in ULPs (units in the last place)
        return fabs(x-y) <= std::numeric_limits<T>::epsilon() * fabs(x+y) * ulp
        // unless the result is subnormal
               || (x-y) < std::numeric_limits<T>::min();
    }
};

template<typename T, unsigned int Dim>
struct almost_equal<phys_tools::autodiff::FD<Dim, T> > {
    using result_type = phys_tools::autodiff::FD<Dim, T>;
    typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
        operator()(result_type x, result_type y, int ulp)
    {
        std::function<bool(T, T)> comp = [&](T x, T y) {
            // the machine epsilon has to be scaled to the magnitude of the values used
            // and multiplied by the desired precision in ULPs (units in the last place)
            return fabs(x-y) <= std::numeric_limits<T>::epsilon() * fabs(x+y) * ulp
            // unless the result is subnormal
                   || fabs(x-y) < std::numeric_limits<T>::min();
        };
        unsigned int n = phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,T>::nVars(x);
        unsigned int m = phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,T>::nVars(y);

        if(n != m)
            return false;

        bool res = comp(x.value(), y.value());
        for(unsigned int i=0; res && i<n; ++i) {
            res &= comp(x.derivative(i), y.derivative(i));
        }

        return res;
    }
};

template<typename T>
struct residual_computer {
    void operator()(std::vector<T>const& z, std::vector<unsigned int>const& n, std::vector<T>const& s, std::vector<unsigned int>const& m, std::vector<T>& res) {
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        //std::cout << std::setprecision(16);
        std::cout << "residual_computer" << std::endl;
        // Useful things to precompute
        const unsigned int N = n.size();
        const unsigned int M = m.size();
        const unsigned int max_n = *std::max_element(n.begin(), n.end());
        const unsigned int max_m = *std::max_element(m.begin(), m.end());
        
        std::vector<T> & c = res;
        c.resize(M*max_m, T(0));
        //std::vector<T> c(M*max_m, T(0));
        #define c(i, j) (c[(i)*max_m + (j)])

        // The algorithm
        for(unsigned int k=0; k<M; ++k) {
            T zero = T(0);
            T one = T(1);
            T denominator = T(0);
            T this_c = T(0);
            bool set_n=false;
            bool set_d=false;
            std::vector<std::pair<T, int> > prod_tokens;
            std::vector<std::pair<T, int> > quot_tokens;
            unsigned int token_pos = 0;
            std::function<void(unsigned int)> acc_num = [&](unsigned int j) {
                T diff = s[k]-z[j];
                int p = n[j];
                T val = pow(diff, p);
                /*
                if(fabs(val) > 1) {
                    prod_tokens.push_back(std::pair<T, int>(diff, p));
                }
                else {
                    quot_tokens.push_back(std::pair<T, int>(diff, p));
                }
                */
                //std::cout << "c *= " << val << std::endl;
                if(set_n) {
                    this_c *= val;
                }
                else {
                    this_c = val;
                    set_n = true;
                }
            };
            std::function<void(unsigned int)> acc_den = [&](unsigned int j) {
                if(j != k) {
                    T diff = s[k]-s[j];
                    //std::cout << "diff: " << diff << std::endl;
                    int p = -int(m[j]);
                    T val(pow(diff, -p));
                    /*
                    if(fabs(val) > 1) {
                        quot_tokens.push_back(std::pair<T, int>(diff, p));
                    }
                    else{
                        prod_tokens.push_back(std::pair<T, int>(diff, p));
                    }
                    */
                    assert(!(val == zero));
                    //std::cout << "c /= " << val << std::endl;
                    if(set_d) {
                        this_c /= val;
                    }
                    else {
                        if(set_d) {
                            denominator *= val;
                        }
                        else {
                            denominator = one / val;
                            set_d = true;
                        }
                    }
                }
            };

            // Naive way
            unsigned int max_NM = std::max(N,M);
            unsigned int min_NM = std::min(N,M);
            for(unsigned int j=0; j<min_NM; ++j) {
                acc_num(j);
                acc_den(j);
            }
            if(max_NM == N) {
                for(unsigned int j=min_NM; j<max_NM; ++j) {
                    acc_num(j);
                }
            }
            else {
                for(unsigned int j=min_NM; j<max_NM; ++j) {
                    acc_den(j);
                }
            }
            assert(set_n);
            assert(set_d);
            if(set_d)
                this_c /= denominator;
            
            //std::cout << "denominator == " << denominator <<std::endl;
            //std::cout << "c[" << k << ", " << m[k]-1 << "] = " << this_c << std::endl;
            
            /*
            // Smart way
            auto prod_it = prod_tokens.begin();
            auto prod_end = prod_tokens.end();
            unsigned int prod_pos = 0;
            auto quot_it = quot_tokens.begin();
            auto quot_end = quot_tokens.end();
            unsigned int quot_pos = 0;

            T new_c(1);
    
            while(prod_it != prod_end || quot_it != quot_end) {
                if((fabs(new_c) < 1 || quot_it == quot_end) && prod_it != prod_end) {
                    // Make it bigger
                    T & val = prod_it->first;
                    int p = prod_it->second;
                    bool pos = p>0;
                    unsigned int index = abs(p);
                    if(prod_pos >= index) {
                        prod_pos = 0;
                        ++prod_it;
                        continue;
                    }
                    else {
                        ++prod_pos;
                    }
                    std::cout << "c = " << new_c << std::endl;
                    std::cout << "Make it bigger" << std::endl;
                    if(pos) {
                        new_c *= val;
                        std::cout <<  "c *="  << val << std::endl;
                    }
                    else {
                        new_c /= val;
                        std::cout <<  "c /="  << val << std::endl;
                    }
                }
                else if(quot_it != quot_end){
                    // Make it smaller
                    T & val = quot_it->first;
                    int p = quot_it->second;
                    bool pos = p>0;
                    unsigned int index = abs(p);
                    if(quot_pos >= index) {
                        quot_pos = 0;
                        ++quot_it;
                        continue;
                    }
                    else {
                        ++quot_pos;
                    }
                    std::cout << "c = " << new_c << std::endl;
                    std::cout << "Make it smaller" << std::endl;
                    if(pos) {
                        new_c *= val;
                        std::cout <<  "c *="  << val << std::endl;
                    }
                    else {
                        new_c /= val;
                        std::cout <<  "c /="  << val << std::endl;
                    }
                }
            }
            */

            //std::cout << "c[" << k << ", " << m[k]-1 << "] = " << new_c << std::endl;
            //std::cout << "new_c[" << k << ", " << m[k]-1 << "] = " << new_c << std::endl;
            //std::cout << "old_c[" << k << ", " << m[k]-1 << "] = " << this_c << std::endl;
            //assert(almost_equal<T>()(new_c, this_c, 2));
            c(k, m[k]-1) = this_c;
            //c(k, m[k]-1) = new_c;
        }
        for(unsigned int k=0; k<M; ++k) {
            std::vector<T> lambda(m[k]-1);
            for(unsigned int i=0; i<m[k]-1; ++i) {
                lambda[i] = T(0);
                for(unsigned int j=0; j<M; ++j) {
                    if(j != k)
                        lambda[i] += m[j]/pow(s[j]-s[k], i+1);
                }
                for(unsigned int j=0; j<N; ++j) {
                    lambda[i] -= n[j]/pow(z[j]-s[k], i+1);
                }
            }
            for(unsigned int L=m[k]-1; L>=1; --L) {
                c(k, L-1) = T(0);
                for(unsigned int i=1; i<=m[k] - L; ++i) {
                    c(k, L-1) += c(k, L+i-1) * lambda[i-1];
                }
                c(k, L-1) = c(k, L-1) / (m[k] - L);
            }
        }

        // Clean up the nasty stuff we defined
		#undef c
        fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        //return c;
    }
};

template<unsigned int digits10, typename T>
struct residual_computer_bignum {
    using big_type=boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<digits10, boost::multiprecision::allocate_stack> >;
    void operator()(std::vector<T>const& z, std::vector<unsigned int>const& n, std::vector<T>const& s, std::vector<unsigned int>const& m, std::vector<big_type>& big_res) {
        std::vector<big_type> big_z(z.begin(), z.end());
        std::vector<big_type> big_s(s.begin(), s.end());
        residual_computer<big_type>()(big_z, n, big_s, m, big_res);
        //std::vector<T> res(big_res.begin(), big_res.end());
        //return big_res;
    }
};

template<unsigned int digits10, typename T>
struct contour_integral_bignum {
    using big_type=boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<digits10, boost::multiprecision::allocate_stack> >;
    big_type operator()(std::vector<T>const& z, std::vector<unsigned int>const& n, std::vector<T>const& s, std::vector<unsigned int>const& m) {
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        const unsigned int M = m.size();
        const unsigned int max_m = *std::max_element(m.begin(), m.end());
        std::vector<big_type> big_z(z.begin(), z.end());
        std::vector<big_type> big_s(s.begin(), s.end());
        std::vector<big_type> c;
        residual_computer<big_type>()(big_z, n, big_s, m, c);
        auto c_range = boost::adaptors::stride(c, max_m);
        std::vector<big_type> c_nonzero(c_range.begin(), c_range.end());
        std::sort(c_nonzero.begin(), c_nonzero.end(), [](big_type const & a, big_type const & b)->bool{return fabs(a) < fabs(b);});
        //auto c_range = iter::slice(c,0,c.size(),max_m);
        //#define c(i, j) (c[(i)*max_m + (j)])

        big_type residual_sum = accumulate(c_nonzero.begin(), c_nonzero.end());

        //for(unsigned int k=0; k<M; ++k) {
        //    residual_sum += c(k,0);
        //}
        //#undef c
        fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        return residual_sum;
    }
};


template<unsigned int digits10, unsigned int Dim, typename T>
struct residual_computer_bignum<digits10, phys_tools::autodiff::FD<Dim, T> > {
    using result_type=phys_tools::autodiff::FD<Dim, T>;
    using big_type=phys_tools::autodiff::FD<Dim, boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<digits10, boost::multiprecision::allocate_stack> > >;
    void operator()(std::vector<result_type>const& z, std::vector<unsigned int>const& n, std::vector<result_type>const& s, std::vector<unsigned int>const& m, std::vector<big_type>& c) {
        std::vector<big_type> big_z(z.begin(), z.end());
        std::vector<big_type> big_s(s.begin(), s.end());
        residual_computer<big_type>()(big_z, n, big_s, m, c);
    }
};

template<unsigned int digits10, unsigned int Dim, typename T>
struct contour_integral_bignum<digits10, phys_tools::autodiff::FD<Dim, T> > {
    using result_type=phys_tools::autodiff::FD<Dim, T>;
    using big_type=phys_tools::autodiff::FD<Dim, boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<digits10, boost::multiprecision::allocate_stack> > >;
    big_type operator()(std::vector<result_type>const& z, std::vector<unsigned int>const& n, std::vector<result_type>const& s, std::vector<unsigned int>const& m) {
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        const unsigned int M = m.size();
        const unsigned int max_m = *std::max_element(m.begin(), m.end());
        std::vector<big_type> big_z(z.begin(), z.end());
        std::vector<big_type> big_s(s.begin(), s.end());
        std::vector<big_type> c;
        residual_computer<big_type>()(big_z, n, big_s, m, c);
        auto c_range = boost::adaptors::stride(c, max_m);
        std::vector<big_type> c_nonzero(c_range.begin(), c_range.end());
        std::sort(c_nonzero.begin(), c_nonzero.end(), [](big_type const & a, big_type const & b)->bool{return fabs(a) < fabs(b);});
        //auto c_range = iter::slice(c,0,c.size(),max_m);
        //#define c(i, j) (c[(i)*max_m + (j)])

        big_type residual_sum = accumulate(c_nonzero.begin(), c_nonzero.end());

        //for(unsigned int k=0; k<M; ++k) {
        //    residual_sum += c(k,0);
        //}
        //#undef c
        fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        return residual_sum;
    }
};

template<typename T>
struct contour_integral_fast {
    T operator()(std::vector<T>const& z, std::vector<unsigned int>const& n, std::vector<T>const& s, std::vector<unsigned int>const& m) {
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        const unsigned int M = m.size();
        const unsigned int max_m = *std::max_element(m.begin(), m.end());
        std::vector<big_type> big_s(s.begin(), s.end());
        std::vector<big_type> c;
    
        fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    }
};

template<typename T>
struct contour_integral {
    T operator()(std::vector<T>const& z, std::vector<unsigned int>const& n, std::vector<T>const& s, std::vector<unsigned int>const& m) {
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        const unsigned int M = m.size();
        const unsigned int max_m = *std::max_element(m.begin(), m.end());
        std::vector<T> c;
        residual_computer<T>()(z, n, s, m, c);
        #define c(i, j) (c[(i)*max_m + (j)])

        T residual_sum(0);

        for(unsigned int k=0; k<M; ++k) {
            residual_sum += c(k,0);
        }
        #undef c
        fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        return residual_sum;
    }
};

} // namespace residuals
} // namespace phys_tools

#endif
