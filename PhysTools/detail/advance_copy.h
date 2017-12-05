#ifndef ADVANCE_COPY_H
#define ADVANCE_COPY_H

#include <iterator>
#include <type_traits>

namespace phys_tools{
namespace detail{
	///Implementation of advance_copy for plain forward iterators
	template <class InputIterator, class Distance>
	InputIterator advance_copy_impl(InputIterator it, Distance n, std::forward_iterator_tag){
		for(; n > 0; --n) ++it;
		return(it);
	}
	///Implementation of advance_copy for bidirectional iterators
	///This is the same as for forward iterators, except that negative distances are handled
	template <class InputIterator, class Distance>
	InputIterator advance_impl(InputIterator it, Distance n, std::bidirectional_iterator_tag){
		if(n >= 0) for(; n > 0; --n) ++it;
		else for(; n < 0; ++n) --it;
		return(it);
	}
	///Implementation of advance_copy for random access iterators
	template <class InputIterator, class Distance>
	InputIterator advance_impl(InputIterator it, Distance n, std::random_access_iterator_tag){
		it+=n;
		return(it);
	}
	
	///Make a copy of an iterator advanced by a specified distance
	///\param it the iterator to copy and advance
	///\param n the distance to advance the copied iterator
	template <class InputIterator, class Distance>
	InputIterator advance_copy(InputIterator it, Distance n){
		static_assert(!std::is_same<typename std::iterator_traits<InputIterator>::iterator_category,
					                std::input_iterator_tag>(),
					  "advance_copy cannot be safely used on input iterators");
		return(advance_impl(it,n,typename std::iterator_traits<InputIterator>::iterator_category()));
	}
}
}

#endif