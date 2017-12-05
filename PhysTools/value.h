/// @file value.hpp
/// @brief Contain a value or the absence of a value.


#ifndef VALUE_HPP
#define VALUE_HPP

//#warning TODO: Replace uses of value with boost::optional

#include <iostream>
#include <stdexcept>


/// @mainpage
/// This is the documentation for the Physics Tools.  That name could
/// perhaps use some work, so it's not set in stone.  At the moment,
/// Physics Tools is by no means a feature complete, solve-all-problems
/// package; rather, it is a generalization of some useful home-grown
/// code for solving particular problems easily.  At present, this
/// includes
/// @li Histogramming
/// @li Plotting
///
/// To get started with histogramming, see @ref Histogramming.  To get
/// started with plotting, see phys_tools::gnuplot.
///
/// At present, the simplest way to get the library source is
/// @code
/// svn checkout http://code.icecube.wisc.edu/svn/sandbox/cweaver/PhysTools/ PhysTools
/// @endcode
/// Of course, you will need a valid IceCube subversion account.
///
/// The library depends on some functionality provided by Boost.
/// Specifically, we use Boost's shared_ptr, tuple, filesystem,
/// iostreams, program options and regex libraries.


/// @brief Namespace for the Physics Tools.
///
/// This namespace contains all of the code making up the Physics Tools.
namespace phys_tools {

/// @brief Store a value or the absence of a value.
///
/// Many situations call for either a value or the explicit absence of
/// a value.  This class is a nearly transparent wrapper for a value,
/// except that a MissingValueError exception is thrown if an attempt
/// is made to use a Value that contains no value.
///
/// A Value can implicitly converted to the underlying value type, or
/// its stored value can be accessed implicitly with the
/// value() method.  It can be cleared with clear(), and you can
/// determine if it currently contains a value with good().  The
/// default constructor yields an empty value; other constructors
/// initialize with an internal stored value.
template<typename value_type>
class Value {
    template <typename T> friend class Value;
public:
    Value ();
    Value (const value_type& val);
//     Value (const Value<value_type>& val);
    template <typename T> Value (const Value<T>& val);
    Value<value_type>&       operator= (const Value<value_type>& v);
    Value<value_type>&       operator+= (const Value<value_type>& v);
    Value<value_type>&       operator-= (const Value<value_type>& v);
    Value<value_type>&       operator*= (const Value<value_type>& v);
    Value<value_type>&       operator/= (const Value<value_type>& v);

    /// @brief Clear the Value.
    void            clear ();

    /// @brief Check whether the Value currently is empty.
    ///
    /// This method is the opposite of Value::good().  It is
    /// provided for convenience and clarity of code, even though it
    /// may, strictly speaking, violate the principle that opposing
    /// methods should have manifestly opposing names.
    bool            empty () const;

    /// @brief Check whether the Value currently contains a value.
    bool            good () const;

    /// @brief Check whether the Value is NaN.
    bool            is_nan () const;

    /// @brief Get a pointer to the stored value.
    ///
    /// This method returns 0 (NULL) if there is no stored value, or
    /// the address of the stored value otherwise.
    value_type*     pointer ();

    /// @brief Get the stored value.
    ///
    /// This method throws an instance of MissingValueError if
    /// there is no stored value when it is called.
    value_type      value () const;

    /// @brief Implicitly convert the Value to the stored value
    /// type.
    ///
    /// This conversion throws an instance of MissingValueError if
    /// there is no stored value when it is called.
    operator value_type () const;

    template<typename T> friend T
    operator+ (const Value<T>& a, const Value<T>& b);
    template<typename T> friend T
    operator+ (const T& a, const Value<T>& b);
    template<typename T> friend T
    operator+ (const Value<T>& a, const T& b);

    template<typename T> friend T
    operator- (const Value<T>& a, const Value<T>& b);
    template<typename T> friend T
    operator- (const T& a, const Value<T>& b);
    template<typename T> friend T
    operator- (const Value<T>& a, const T& b);

    template<typename T> friend T
    operator* (const Value<T>& a, const Value<T>& b);
    template<typename T> friend T
    operator* (const T& a, const Value<T>& b);
    template<typename T> friend T
    operator* (const Value<T>& a, const T& b);

    template<typename T> friend T
    operator/ (const Value<T>& a, const Value<T>& b);
    template<typename T> friend T
    operator/ (const T& a, const Value<T>& b);
    template<typename T> friend T
    operator/ (const Value<T>& a, const T& b);

    template<typename T> friend bool
    operator< (const Value<T>& a, const Value<T>& b);
    template<typename T> friend bool
    operator< (const T& a, const Value<T>& b);
    template<typename T> friend bool
    operator< (const Value<T>& a, const T& b);

    template<typename T> friend bool
    operator<= (const Value<T>& a, const Value<T>& b);
    template<typename T> friend bool
    operator<= (const T& a, const Value<T>& b);
    template<typename T> friend bool
    operator<= (const Value<T>& a, const T& b);

    template<typename T> friend bool
    operator> (const Value<T>& a, const Value<T>& b);
    template<typename T> friend bool
    operator> (const T& a, const Value<T>& b);
    template<typename T> friend bool
    operator> (const Value<T>& a, const T& b);

    template<typename T> friend bool
    operator>= (const Value<T>& a, const Value<T>& b);
    template<typename T> friend bool
    operator>= (const T& a, const Value<T>& b);
    template<typename T> friend bool
    operator>= (const Value<T>& a, const T& b);

    template<typename T> friend bool
    operator== (const Value<T>& a, const Value<T>& b);
    template<typename T> friend bool
    operator== (const T& a, const Value<T>& b);
    template<typename T> friend bool
    operator== (const Value<T>& a, const T& b);

    template<typename T> friend bool
    operator!= (const Value<T>& a, const Value<T>& b);
    template<typename T> friend bool
    operator!= (const T& a, const Value<T>& b);
    template<typename T> friend bool
    operator!= (const Value<T>& a, const T& b);

    template<typename T> friend std::ostream&
    operator<< (std::ostream& os, const Value<T>& v);


private:

    void            check () const;
    value_type      _val;
    bool            _good;
};

/// @brief Indicate that an empty value was used.
class MissingValueError : public std::logic_error {
public:
    MissingValueError ()
        :   logic_error ("cannot use an empty Value")
    {
    }
};

} // namespace phys_tools

#include "value.tcpp"
#endif  /* VALUE_HPP */
