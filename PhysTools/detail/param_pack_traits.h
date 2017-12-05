#ifndef PARAM_PACK_TRAITS_H
#define PARAM_PACK_TRAITS_H

namespace phys_tools{
namespace detail{

///Meta-function which negates integral_constants
template<typename T>
class not_ : public std::integral_constant<bool,!T::value>{};

///Meta-function which conditionally ignores its input
///This is used to avoid 'evaluating' the second argumets to and_ and or_ when short circuiting
///When its first argument is false, ignore_if inherits its second argument,
///otherwise, it inherits from an integral_constant whose value is expected to be irrelevant
template<bool b, typename T>
class ignore_if : public T{};
template<typename T>
class ignore_if<true,T> : public std::true_type{};

///Meta-function which computes the conjunction of two integral constants, with short-circuiting
template<typename Left, typename Right>
class and_ : public std::integral_constant<bool,Left::value && ignore_if<!Left::value,Right>::value>{};

///Meta-function which computes the disjunction of two integral constants, with short-circuiting
template<typename Left, typename Right>
class or_ : public std::integral_constant<bool,Left::value || ignore_if<Left::value,Right>::value>{};

///Trait class that checks that all of the other types it is given are the same as ExpectedType
template<typename ExpectedType, typename FirstArg, typename... Args>
class all_args_are : public and_<std::is_same<ExpectedType,FirstArg>,all_args_are<ExpectedType,Args...>>{};
///Base case of all_args_are (just std::is_same)
template<typename ExpectedType, typename LastArg>
class all_args_are<ExpectedType,LastArg> : public std::is_same<ExpectedType,LastArg>{};
	
///Trait class that checks that all of the other types it is given are the same as ExpectedType,
///except for the last, which is allowed to be anything
template<typename ExpectedType, typename FirstArg, typename... Args>
class all_args_except_last_are : public and_<std::is_same<ExpectedType,FirstArg>,all_args_except_last_are<ExpectedType,Args...>>{};
///Base case of all_args_except_last_are, always true
template<typename ExpectedType, typename LastArg>
class all_args_except_last_are<ExpectedType,LastArg> : public std::true_type{};

///Trait class that identifies whether the final entry in Args is ExpectedType
template<typename ExpectedType, typename FirstArg, typename... Args>
class last_arg_is : public last_arg_is<ExpectedType, Args...>{};
///Failure case specialization of last_arg_is
template<typename ExpectedType, typename FinalArg>
class last_arg_is<ExpectedType, FinalArg> : public std::false_type{};
///Success case specialization of last_arg_is
template<typename ExpectedType>
class last_arg_is<ExpectedType, ExpectedType> : public std::true_type{};
	
	
///Trait class that checks that all of the other types it is given are convertible to ExpectedType
template<typename ExpectedType, typename FirstArg, typename... Args>
class all_args_convertible : public and_<std::is_convertible<FirstArg,ExpectedType>,all_args_convertible<ExpectedType,Args...>>{};
///Base case of all_args_convertible (just std::is_convertible)
template<typename ExpectedType, typename LastArg>
class all_args_convertible<ExpectedType,LastArg> : public std::is_convertible<LastArg,ExpectedType>{};

///Trait class that checks that all of the other types it is given are convertible to ExpectedType,
///except for the last, which is allowed to be anything
template<typename ExpectedType, typename FirstArg, typename... Args>
class all_args_except_last_convertible : public and_<std::is_convertible<FirstArg,ExpectedType>,all_args_except_last_convertible<ExpectedType,Args...>>{};
///Base case of all_args_except_last_convertible, always true
template<typename ExpectedType, typename LastArg>
class all_args_except_last_convertible<ExpectedType,LastArg> : public std::true_type{};

///Trait class that identifies whether the final entry in Args is convertible to ExpectedType
template<typename ExpectedType, typename FirstArg, typename... Args>
class last_arg_convertible : public last_arg_convertible<ExpectedType, Args...>{};
///Base case of last_arg_convertible (just std::is_convertible)
template<typename ExpectedType, typename LastArg>
class last_arg_convertible<ExpectedType, LastArg> : public std::is_convertible<LastArg,ExpectedType>{};

///Trait class that identifies whether the first entry in Args is convertible to ExpectedType
template<typename ExpectedType, typename FirstArg, typename... Args>
class first_arg_convertible : public std::is_convertible<FirstArg,ExpectedType>{};
	
///Trait class that identifies whether the first entry in Args is a pointer
template<typename FirstArg, typename... Args>
class first_is_pointer : public std::is_pointer<FirstArg>{};

}
}

#endif