#ifndef __BRENT_H__
#define __BRENT_H__

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <functional>

using namespace std;

#include <vector>
namespace brent {

template<typename T>
class func_base{
    public:
        virtual T operator() (T) = 0;
};

template<typename T>
class monicPoly : public func_base<T> {
    public:
        std::vector<T> coeff;
        virtual T operator() (T x);
        // constructors:
        monicPoly(const size_t degree)
            : coeff(degree) {}
        monicPoly(const std::vector<T>& v)
            : coeff(v) {}
        monicPoly(const T* c, size_t degree)
            : coeff(std::vector<T>(c, c+degree)) {}
};

template<typename T>
class Poly : public func_base<T> {
    public:
        std::vector<T> coeff;    // a vector of size nterms i.e. 1+degree
        virtual T operator() (T x);
        // constructors:
        Poly(const size_t degree)
            : coeff(1+degree) {}
        Poly(const std::vector<T>& v)
            : coeff(v) {}
        Poly(const T* c, size_t degree)
            : coeff(std::vector<T>(c, 1+c+degree)) {}
};

template<typename T>
T glomin ( T a, T b, T c, T m, T e, T t,
        func_base<T>& f, T &x );
template<typename T>
T local_min ( T a, T b, T t, func_base<T>& f,
        T &x );
template<typename T>
T local_min_rc ( T &a, T &b, int &status, T value );
template<typename T>
T r8_epsilon ( );
template<typename T>
T r8_max ( T x, T y );
template<typename T>
T r8_sign ( T x );
template<typename T>
void timestamp ( );
template<typename T>
T zero ( T a, T b, T t, func_base<T>& f );
template<typename T>
void zero_rc ( T a, T b, T t, T &arg, int &status,
        T value );

// === simple wrapper functions
// === for convenience and/or compatibility
template<typename T>
T glomin ( T a, T b, T c, T m, T e, T t,
        T f ( T x ), T &x );
template<typename T>
T local_min ( T a, T b, T t, T f ( T x ),
        T &x );
template<typename T>
T zero ( T a, T b, T t, T f ( T x ) );

//****************************************************************************80

template<typename T> 
T glomin ( T a, T b, T c, T m, T e, T t,
        func_base<T>& f, T &x )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
    //
    //  Discussion:
    //
    //    This function assumes that F(X) is twice continuously differentiable
    //    over [A,B] and that F''(X) <= M for all X in [A,B].
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    17 July 2011
    //
    //  Author:
    //
    //    Original FORTRAN77 version by Richard Brent.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //
    //    Richard Brent,
    //    Algorithms for Minimization Without Derivatives,
    //    Dover, 2002,
    //    ISBN: 0-486-41998-3,
    //    LC: QA402.5.B74.
    //
    //  Parameters:
    //
    //    Input, T A, B, the endpoints of the interval.
    //    It must be the case that A < B.
    //
    //    Input, T C, an initial guess for the global
    //    minimizer.  If no good guess is known, C = A or B is acceptable.
    //
    //    Input, T M, the bound on the second derivative.
    //
    //    Input, T E, a positive tolerance, a bound for the
    //    absolute error in the evaluation of F(X) for any X in [A,B].
    //
    //    Input, T T, a positive error tolerance.
    //
    //    Input, func_base& F, a user-supplied c++ functor whose
    //    global minimum is being sought.  The input and output
    //    of F() are of type T.
    //
    //    Output, T &X, the estimated value of the abscissa
    //    for which F attains its global minimum value in [A,B].
    //
    //    Output, T GLOMIN, the value F(X).
    //
{
    T a0;
    T a2;
    T a3;
    T d0;
    T d1;
    T d2;
    T h;
    int k;
    T m2;
    T macheps;
    T p;
    T q;
    T qs;
    T r;
    T s;
    T sc;
    T y;
    T y0;
    T y1;
    T y2;
    T y3;
    T yb;
    T z0;
    T z1;
    T z2;

    a0 = b;
    x = a0;
    a2 = a;
    y0 = f ( b );
    yb = y0;
    y2 = f ( a );
    y = y2;

    if ( y0 < y )
    {
        y = y0;
    }
    else
    {
        x = a;
    }

    if ( m <= 0.0 || b <= a )
    {
        return y;
    }

    macheps = r8_epsilon<T> ( );

    m2 = 0.5 * ( 1.0 + 16.0 * macheps ) * m;

    if ( c <= a || b <= c )
    {
        sc = 0.5 * ( a + b );
    }
    else
    {
        sc = c;
    }

    y1 = f ( sc );
    k = 3;
    d0 = a2 - sc;
    h = 9.0 / 11.0;

    if ( y1 < y )
    {
        x = sc;
        y = y1;
    }
    //
    //  Loop.
    //
    for ( ; ; )
    {
        d1 = a2 - a0;
        d2 = sc - a0;
        z2 = b - a2;
        z0 = y2 - y1;
        z1 = y2 - y0;
        r = d1 * d1 * z0 - d0 * d0 * z1;
        p = r;
        qs = 2.0 * ( d0 * z1 - d1 * z0 );
        q = qs;

        if ( k < 1000000 || y2 <= y )
        {
            for ( ; ; )
            {
                if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) <
                        z2 * m2 * r * ( z2 * q - r ) )
                {
                    a3 = a2 + r / q;
                    y3 = f ( a3 );

                    if ( y3 < y )
                    {
                        x = a3;
                        y = y3;
                    }
                }
                k = ( ( 1611 * k ) % 1048576 );
                q = 1.0;
                r = ( b - a ) * 0.00001 * ( T ) ( k );

                if ( z2 <= r )
                {
                    break;
                }
            }
        }
        else
        {
            k = ( ( 1611 * k ) % 1048576 );
            q = 1.0;
            r = ( b - a ) * 0.00001 * ( T ) ( k );

            while ( r < z2 )
            {
                if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) <
                        z2 * m2 * r * ( z2 * q - r ) )
                {
                    a3 = a2 + r / q;
                    y3 = f ( a3 );

                    if ( y3 < y )
                    {
                        x = a3;
                        y = y3;
                    }
                }
                k = ( ( 1611 * k ) % 1048576 );
                q = 1.0;
                r = ( b - a ) * 0.00001 * ( T ) ( k );
            }
        }

        r = m2 * d0 * d1 * d2;
        s = sqrt ( ( ( y2 - y ) + t ) / m2 );
        h = 0.5 * ( 1.0 + h );
        p = h * ( p + 2.0 * r * s );
        q = q + 0.5 * qs;
        r = - 0.5 * ( d0 + ( z0 + 2.01 * e ) / ( d0 * m2 ) );

        if ( r < s || d0 < 0.0 )
        {
            r = a2 + s;
        }
        else
        {
            r = a2 + r;
        }

        if ( 0.0 < p * q )
        {
            a3 = a2 + p / q;
        }
        else
        {
            a3 = r;
        }

        for ( ; ; )
        {
            a3 = r8_max ( a3, r );

            if ( b <= a3 )
            {
                a3 = b;
                y3 = yb;
            }
            else
            {
                y3 = f ( a3 );
            }

            if ( y3 < y )
            {
                x = a3;
                y = y3;
            }

            d0 = a3 - a2;

            if ( a3 <= r )
            {
                break;
            }

            p = 2.0 * ( y2 - y3 ) / ( m * d0 );

            if ( ( 1.0 + 9.0 * macheps ) * d0 <= fabs ( p ) )
            {
                break;
            }

            if ( 0.5 * m2 * ( d0 * d0 + p * p ) <= ( y2 - y ) + ( y3 - y ) + 2.0 * t )
            {
                break;
            }
            a3 = 0.5 * ( a2 + a3 );
            h = 0.9 * h;
        }

        if ( b <= a3 )
        {
            break;
        }

        a0 = sc;
        sc = a2;
        a2 = a3;
        y0 = y1;
        y1 = y2;
        y2 = y3;
    }

    return y;
}
//****************************************************************************80

template<typename T> 
T local_min ( T a, T b, T t, func_base<T>& f,
        T &x )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
    //
    //  Discussion:
    //
    //    The method used is a combination of golden section search and
    //    successive parabolic interpolation.  Convergence is never much slower
    //    than that for a Fibonacci search.  If F has a continuous second
    //    derivative which is positive at the minimum (which is not at A or
    //    B), then convergence is superlinear, and usually of the order of
    //    about 1.324....
    //
    //    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
    //    F is never evaluated at two points closer than TOL.
    //
    //    If F is a unimodal function and the computed values of F are always
    //    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
    //    LOCAL_MIN approximates the abscissa of the global minimum of F on the
    //    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
    //
    //    If F is not unimodal, then LOCAL_MIN may approximate a local, but
    //    perhaps non-global, minimum to the same accuracy.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    17 July 2011
    //
    //  Author:
    //
    //    Original FORTRAN77 version by Richard Brent.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //
    //    Richard Brent,
    //    Algorithms for Minimization Without Derivatives,
    //    Dover, 2002,
    //    ISBN: 0-486-41998-3,
    //    LC: QA402.5.B74.
    //
    //  Parameters:
    //
    //    Input, T A, B, the endpoints of the interval.
    //
    //    Input, T T, a positive absolute error tolerance.
    //
    //    Input, func_base& F, a user-supplied c++ functor whose
    //    local minimum is being sought.  The input and output
    //    of F() are of type T.
    //
    //    Output, T &X, the estimated value of an abscissa
    //    for which F attains a local minimum value in [A,B].
    //
    //    Output, T LOCAL_MIN, the value F(X).
    //
{
    T c;
    T d;
    T e;
    T eps;
    T fu;
    T fv;
    T fw;
    T fx;
    T m;
    T p;
    T q;
    T r;
    T sa;
    T sb;
    T t2;
    T tol;
    T u;
    T v;
    T w;
    //
    //  C is the square of the inverse of the golden ratio.
    //
    c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

    eps = sqrt ( r8_epsilon<T> ( ) );

    sa = a;
    sb = b;
    x = sa + c * ( b - a );
    w = x;
    v = w;
    e = 0.0;
    fx = f ( x );
    fw = fx;
    fv = fw;

    for ( ; ; )
    {
        m = 0.5 * ( sa + sb ) ;
        tol = eps * fabs ( x ) + t;
        t2 = 2.0 * tol;
        //
        //  Check the stopping criterion.
        //
        if ( fabs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
        {
            break;
        }
        //
        //  Fit a parabola.
        //
        r = 0.0;
        q = r;
        p = q;

        if ( tol < fabs ( e ) )
        {
            r = ( x - w ) * ( fx - fv );
            q = ( x - v ) * ( fx - fw );
            p = ( x - v ) * q - ( x - w ) * r;
            q = 2.0 * ( q - r );
            if ( 0.0 < q )
            {
                p = - p;
            }
            q = fabs ( q );
            r = e;
            e = d;
        }

        if ( fabs ( p ) < fabs ( 0.5 * q * r ) &&
                q * ( sa - x ) < p &&
                p < q * ( sb - x ) )
        {
            //
            //  Take the parabolic interpolation step.
            //
            d = p / q;
            u = x + d;
            //
            //  F must not be evaluated too close to A or B.
            //
            if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
            {
                if ( x < m )
                {
                    d = tol;
                }
                else
                {
                    d = - tol;
                }
            }
        }
        //
        //  A golden-section step.
        //
        else
        {
            if ( x < m )
            {
                e = sb - x;
            }
            else
            {
                e = sa - x;
            }
            d = c * e;
        }
        //
        //  F must not be evaluated too close to X.
        //
        if ( tol <= fabs ( d ) )
        {
            u = x + d;
        }
        else if ( 0.0 < d )
        {
            u = x + tol;
        }
        else
        {
            u = x - tol;
        }

        fu = f ( u );
        //
        //  Update A, B, V, W, and X.
        //
        if ( fu <= fx )
        {
            if ( u < x )
            {
                sb = x;
            }
            else
            {
                sa = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if ( u < x )
            {
                sa = u;
            }
            else
            {
                sb = u;
            }

            if ( fu <= fw || w == x )
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if ( fu <= fv || v == x || v == w )
            {
                v = u;
                fv = fu;
            }
        }
    }
    return fx;
}
//****************************************************************************80

template<typename T> 
T local_min_rc ( T &a, T &b, int &status, T value )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
    //
    //  Discussion:
    //
    //    This routine seeks an approximation to the point where a function
    //    F attains a minimum on the interval (A,B).
    //
    //    The method used is a combination of golden section search and
    //    successive parabolic interpolation.  Convergence is never much
    //    slower than that for a Fibonacci search.  If F has a continuous
    //    second derivative which is positive at the minimum (which is not
    //    at A or B), then convergence is superlinear, and usually of the
    //    order of about 1.324...
    //
    //    The routine is a revised version of the Brent local minimization
    //    algorithm, using reverse communication.
    //
    //    It is worth stating explicitly that this routine will NOT be
    //    able to detect a minimizer that occurs at either initial endpoint
    //    A or B.  If this is a concern to the user, then the user must
    //    either ensure that the initial interval is larger, or to check
    //    the function value at the returned minimizer against the values
    //    at either endpoint.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    17 July 2011
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    Richard Brent,
    //    Algorithms for Minimization Without Derivatives,
    //    Dover, 2002,
    //    ISBN: 0-486-41998-3,
    //    LC: QA402.5.B74.
    //
    //    David Kahaner, Cleve Moler, Steven Nash,
    //    Numerical Methods and Software,
    //    Prentice Hall, 1989,
    //    ISBN: 0-13-627258-4,
    //    LC: TA345.K34.
    //
    //  Parameters
    //
    //    Input/output, T &A, &B.  On input, the left and right
    //    endpoints of the initial interval.  On output, the lower and upper
    //    bounds for an interval containing the minimizer.  It is required
    //    that A < B.
    //
    //    Input/output, int &STATUS, used to communicate between
    //    the user and the routine.  The user only sets STATUS to zero on the first
    //    call, to indicate that this is a startup call.  The routine returns STATUS
    //    positive to request that the function be evaluated at ARG, or returns
    //    STATUS as 0, to indicate that the iteration is complete and that
    //    ARG is the estimated minimizer.
    //
    //    Input, T VALUE, the function value at ARG, as requested
    //    by the routine on the previous call.
    //
    //    Output, T LOCAL_MIN_RC, the currently considered point.
    //    On return with STATUS positive, the user is requested to evaluate the
    //    function at this point, and return the value in VALUE.  On return with
    //    STATUS zero, this is the routine's estimate for the function minimizer.
    //
    //  Local parameters:
    //
    //    C is the squared inverse of the golden ratio.
    //
    //    EPS is the square root of the relative machine precision.
    //
{
    static T arg;
    static T c;
    static T d;
    static T e;
    static T eps;
    static T fu;
    static T fv;
    static T fw;
    static T fx;
    static T midpoint;
    static T p;
    static T q;
    static T r;
    static T tol;
    static T tol1;
    static T tol2;
    static T u;
    static T v;
    static T w;
    static T x;
    //
    //  STATUS (INPUT) = 0, startup.
    //
    if ( status == 0 )
    {
        if ( b <= a )
        {
            cout << "\n";
            cout << "LOCAL_MIN_RC - Fatal error!\n";
            cout << "  A < B is required, but\n";
            cout << "  A = " << a << "\n";
            cout << "  B = " << b << "\n";
            status = -1;
            exit ( 1 );
        }
        c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

        eps = sqrt ( r8_epsilon<T> ( ) );
        tol = r8_epsilon<T> ( );

        v = a + c * ( b - a );
        w = v;
        x = v;
        e = 0.0;

        status = 1;
        arg = x;

        return arg;
    }
    //
    //  STATUS (INPUT) = 1, return with initial function value of FX.
    //
    else if ( status == 1 )
    {
        fx = value;
        fv = fx;
        fw = fx;
    }
    //
    //  STATUS (INPUT) = 2 or more, update the data.
    //
    else if ( 2 <= status )
    {
        fu = value;

        if ( fu <= fx )
        {
            if ( x <= u )
            {
                a = x;
            }
            else
            {
                b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if ( u < x )
            {
                a = u;
            }
            else
            {
                b = u;
            }

            if ( fu <= fw || w == x )
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if ( fu <= fv || v == x || v == w )
            {
                v = u;
                fv = fu;
            }
        }
    }
    //
    //  Take the next step.
    //
    midpoint = 0.5 * ( a + b );
    tol1 = eps * fabs ( x ) + tol / 3.0;
    tol2 = 2.0 * tol1;
    //
    //  If the stopping criterion is satisfied, we can exit.
    //
    if ( fabs ( x - midpoint ) <= ( tol2 - 0.5 * ( b - a ) ) )
    {
        status = 0;
        return arg;
    }
    //
    //  Is golden-section necessary?
    //
    if ( fabs ( e ) <= tol1 )
    {
        if ( midpoint <= x )
        {
            e = a - x;
        }
        else
        {
            e = b - x;
        }
        d = c * e;
    }
    //
    //  Consider fitting a parabola.
    //
    else
    {
        r = ( x - w ) * ( fx - fv );
        q = ( x - v ) * ( fx - fw );
        p = ( x - v ) * q - ( x - w ) * r;
        q = 2.0 * ( q - r );
        if ( 0.0 < q )
        {
            p = - p;
        }
        q = fabs ( q );
        r = e;
        e = d;
        //
        //  Choose a golden-section step if the parabola is not advised.
        //
        if (
                ( fabs ( 0.5 * q * r ) <= fabs ( p ) ) ||
                ( p <= q * ( a - x ) ) ||
                ( q * ( b - x ) <= p ) )
        {
            if ( midpoint <= x )
            {
                e = a - x;
            }
            else
            {
                e = b - x;
            }
            d = c * e;
        }
        //
        //  Choose a parabolic interpolation step.
        //
        else
        {
            d = p / q;
            u = x + d;

            if ( ( u - a ) < tol2 )
            {
                d = tol1 * r8_sign ( midpoint - x );
            }

            if ( ( b - u ) < tol2 )
            {
                d = tol1 * r8_sign ( midpoint - x );
            }
        }
    }
    //
    //  F must not be evaluated too close to X.
    //
    if ( tol1 <= fabs ( d ) )
    {
        u = x + d;
    }
    if ( fabs ( d ) < tol1 )
    {
        u = x + tol1 * r8_sign ( d );
    }
    //
    //  Request value of F(U).
    //
    arg = u;
    status = status + 1;

    return arg;
}
//****************************************************************************80

template<typename T> 
T r8_epsilon ( )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8_EPSILON returns the R8 roundoff unit.
    //
    //  Discussion:
    //
    //    The roundoff unit is a number R which is a power of 2 with the
    //    property that, to the precision of the computer's arithmetic,
    //      1 < 1 + R
    //    but
    //      1 = ( 1 + R / 2 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    01 September 2012
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Output, T R8_EPSILON, the R8 round-off unit.
    //
{
    const T value = 2.220446049250313E-016;

    return value;
}
//****************************************************************************80

template<typename T> 
T r8_max ( T x, T y )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8_MAX returns the maximum of two R8's.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    18 August 2004
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, T X, Y, the quantities to compare.
    //
    //    Output, T R8_MAX, the maximum of X and Y.
    //
{
    T value;

    if ( y < x )
    {
        value = x;
    }
    else
    {
        value = y;
    }
    return value;
}
//****************************************************************************80

template<typename T> 
T r8_sign ( T x )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8_SIGN returns the sign of an R8.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    18 October 2004
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, T X, the number whose sign is desired.
    //
    //    Output, T R8_SIGN, the sign of X.
    //
{
    T value;

    if ( x < 0.0 )
    {
        value = -1.0;
    }
    else
    {
        value = 1.0;
    }
    return value;
}
//****************************************************************************80

template<typename T>
void timestamp ( )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    TIMESTAMP prints the current YMDHMS date as a time stamp.
    //
    //  Example:
    //
    //    31 May 2001 09:45:54 AM
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    24 September 2003
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    None
    //
{
    const int TIME_SIZE(40);

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    cout << time_buffer << "\n";

    return;
}
//****************************************************************************80

template<typename T> 
T zero ( T a, T b, T t, func_base<T>& f )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    ZERO seeks the root of a function F(X) in an interval [A,B].
    //
    //  Discussion:
    //
    //    The interval [A,B] must be a change of sign interval for F.
    //    That is, F(A) and F(B) must be of opposite signs.  Then
    //    assuming that F is continuous implies the existence of at least
    //    one value C between A and B for which F(C) = 0.
    //
    //    The location of the zero is determined to within an accuracy
    //    of 6 * MACHEPS * abs ( C ) + 2 * T.
    //
    //    Thanks to Thomas Secretin for pointing out a transcription error in the
    //    setting of the value of P, 11 February 2013.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    11 February 2013
    //
    //  Author:
    //
    //    Original FORTRAN77 version by Richard Brent.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //
    //    Richard Brent,
    //    Algorithms for Minimization Without Derivatives,
    //    Dover, 2002,
    //    ISBN: 0-486-41998-3,
    //    LC: QA402.5.B74.
    //
    //  Parameters:
    //
    //    Input, T A, B, the endpoints of the change of sign interval.
    //
    //    Input, T T, a positive error tolerance.
    //
    //    Input, func_base& F, the name of a user-supplied c++ functor
    //    whose zero is being sought.  The input and output
    //    of F() are of type T.
    //
    //    Output, T ZERO, the estimated value of a zero of
    //    the function F.
    //
{
    T c;
    T d;
    T e;
    T fa;
    T fb;
    T fc;
    T m;
    T macheps;
    T p;
    T q;
    T r;
    T s;
    T sa;
    T sb;
    T tol;
    //
    //  Make local copies of A and B.
    //
    sa = a;
    sb = b;
    fa = f ( sa );
    fb = f ( sb );

    c = sa;
    fc = fa;
    e = sb - sa;
    d = e;

    macheps = r8_epsilon<T> ( );

    for ( ; ; )
    {
        if ( fabs ( fc ) < fabs ( fb ) )
        {
            sa = sb;
            sb = c;
            c = sa;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol = 2.0 * macheps * fabs ( sb ) + t;
        m = 0.5 * ( c - sb );

        if ( fabs ( m ) <= tol || fb == 0.0 )
        {
            break;
        }

        if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
        {
            e = m;
            d = e;
        }
        else
        {
            s = fb / fa;

            if ( sa == c )
            {
                p = 2.0 * m * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
                q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
            }

            if ( 0.0 < p )
            {
                q = - q;
            }
            else
            {
                p = - p;
            }

            s = e;
            e = d;

            if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
                    p < fabs ( 0.5 * s * q ) )
            {
                d = p / q;
            }
            else
            {
                e = m;
                d = e;
            }
        }
        sa = sb;
        fa = fb;

        if ( tol < fabs ( d ) )
        {
            sb = sb + d;
        }
        else if ( 0.0 < m )
        {
            sb = sb + tol;
        }
        else
        {
            sb = sb - tol;
        }

        fb = f ( sb );

        if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
        {
            c = sa;
            fc = fa;
            e = sb - sa;
            d = e;
        }
    }
    return sb;
}
//****************************************************************************80

template<typename T> 
void zero_rc ( T a, T b, T t, T &arg, int &status,
        T value )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    ZERO_RC seeks the root of a function F(X) using reverse communication.
    //
    //  Discussion:
    //
    //    The interval [A,B] must be a change of sign interval for F.
    //    That is, F(A) and F(B) must be of opposite signs.  Then
    //    assuming that F is continuous implies the existence of at least
    //    one value C between A and B for which F(C) = 0.
    //
    //    The location of the zero is determined to within an accuracy
    //    of 6 * MACHEPS * abs ( C ) + 2 * T.
    //
    //    The routine is a revised version of the Brent zero finder
    //    algorithm, using reverse communication.
    //
    //    Thanks to Thomas Secretin for pointing out a transcription error in the
    //    setting of the value of P, 11 February 2013.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    11 February 2013
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    Richard Brent,
    //    Algorithms for Minimization Without Derivatives,
    //    Dover, 2002,
    //    ISBN: 0-486-41998-3,
    //    LC: QA402.5.B74.
    //
    //  Parameters:
    //
    //    Input, T A, B, the endpoints of the change of sign interval.
    //
    //    Input, T T, a positive error tolerance.
    //
    //    Output, T &ARG, the currently considered point.  The user
    //    does not need to initialize this value.  On return with STATUS positive,
    //    the user is requested to evaluate the function at ARG, and return
    //    the value in VALUE.  On return with STATUS zero, ARG is the routine's
    //    estimate for the function's zero.
    //
    //    Input/output, int &STATUS, used to communicate between
    //    the user and the routine.  The user only sets STATUS to zero on the first
    //    call, to indicate that this is a startup call.  The routine returns STATUS
    //    positive to request that the function be evaluated at ARG, or returns
    //    STATUS as 0, to indicate that the iteration is complete and that
    //    ARG is the estimated zero
    //
    //    Input, T VALUE, the function value at ARG, as requested
    //    by the routine on the previous call.
    //
{
    static T c;
    static T d;
    static T e;
    static T fa;
    static T fb;
    static T fc;
    T m;
    static T macheps;
    T p;
    T q;
    T r;
    T s;
    static T sa;
    static T sb;
    T tol;
    //
    //  Input STATUS = 0.
    //  Initialize, request F(A).
    //
    if ( status == 0 )
    {
        macheps = r8_epsilon<T> ( );

        sa = a;
        sb = b;
        e = sb - sa;
        d = e;

        status = 1;
        arg = a;
        return;
    }
    //
    //  Input STATUS = 1.
    //  Receive F(A), request F(B).
    //
    else if ( status == 1 )
    {
        fa = value;
        status = 2;
        arg = sb;
        return;
    }
    //
    //  Input STATUS = 2
    //  Receive F(B).
    //
    else if ( status == 2 )
    {
        fb = value;

        if ( 0.0 < fa * fb )
        {
            status = -1;
            return;
        }
        c = sa;
        fc = fa;
    }
    else
    {
        fb = value;

        if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
        {
            c = sa;
            fc = fa;
            e = sb - sa;
            d = e;
        }
    }
    //
    //  Compute the next point at which a function value is requested.
    //
    if ( fabs ( fc ) < fabs ( fb ) )
    {
        sa = sb;
        sb = c;
        c = sa;
        fa = fb;
        fb = fc;
        fc = fa;
    }

    tol = 2.0 * macheps * fabs ( sb ) + t;
    m = 0.5 * ( c - sb );

    if ( fabs ( m ) <= tol || fb == 0.0 )
    {
        status = 0;
        arg = sb;
        return;
    }

    if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
    {
        e = m;
        d = e;
    }
    else
    {
        s = fb / fa;

        if ( sa == c )
        {
            p = 2.0 * m * s;
            q = 1.0 - s;
        }
        else
        {
            q = fa / fc;
            r = fb / fc;
            p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
            q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
        }

        if ( 0.0 < p )
        {
            q = - q;
        }
        else
        {
            p = - p;
        }
        s = e;
        e = d;

        if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
                p < fabs ( 0.5 * s * q ) )
        {
            d = p / q;
        }
        else
        {
            e = m;
            d = e;
        }
    }

    sa = sb;
    fa = fb;

    if ( tol < fabs ( d ) )
    {
        sb = sb + d;
    }
    else if ( 0.0 < m )
    {
        sb = sb + tol;
    }
    else
    {
        sb = sb - tol;
    }

    arg = sb;
    status = status + 1;

    return;
}

// ======================================================================
// === Simple wrapper functions
// === for convenience and/or compatibility.
//
// === The three functions are the same as above,
// === except that they take a plain function F
// === instead of a c++ functor.  In all cases, the
// === input and output of F() are of type T.


template<typename T> 
class std_func_wrapper : public func_base<T> {
    std::function<T(T)> func;
    public:
    std_func_wrapper(std::function<T(T)> f):func(f) {}
    T operator() (T x){
        return func(x);
    }
};

template<typename T> 
class func_wrapper : public func_base<T> {
    typedef T TOfT (T);
    TOfT* func;
    public:
    func_wrapper(TOfT* f) {
        func = f;
    }
    virtual T operator() (T x){
        return func(x);
    }
};

//****************************************************************************80

template<typename T> 
T glomin ( T a, T b, T c, T m, T e,
        T t, T f ( T x ), T &x ){
    func_wrapper<T> foo(f);
    return glomin<T>(a, b, c, m, e, t, foo, x);
}

//****************************************************************************80

template<typename T> 
T local_min ( T a, T b, T t, T f ( T x ),
        T &x ){
    func_wrapper<T> foo(f);
    return local_min<T>(a, b, t, foo, x);
}

//****************************************************************************80

template<typename T> 
T zero ( T a, T b, T t, std::function<T(T)> f){
    std_func_wrapper<T> foo(f);
    return zero<T>(a, b, t, foo);
}

template<typename T> 
T zero ( T a, T b, T t, T f ( T x ) ){
    func_wrapper<T> foo(f);
    return zero<T>(a, b, t, foo);
}

// ======================================================================
// Generally useful functor to evaluate a monic polynomial.
// For details, see class definition in brent.hpp

template<typename T> 
T monicPoly<T>::operator()(T x){
    T rslt(1);
    for (int ii = coeff.size()-1; ii >= 0; ii--){
        rslt *= x;
        rslt += coeff[ii];
    }
    return rslt;
}

// Similarly, evaluate a general polynomial (not necessarily monic):
template<typename T> 
T Poly<T>::operator()(T x){
    T rslt(0);
    for (int ii = coeff.size()-1; ii >= 0; ii--){
        rslt *= x;
        rslt += coeff[ii];
    }
    return rslt;
}

} // end namespace brent
#endif
