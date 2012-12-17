/*!
 * @file dsp++/algorithm.h
 * @brief Various unassorted algorithms loosely complementing these found in &lt;algorithm&gt;
 * @author Andrzej Ciarkowski <mailto:andrzej.ciarkowski@gmail.com>
 */
#ifndef DSP_ALGORITHM_H_INCLUDED
#define DSP_ALGORITHM_H_INCLUDED

#include <dsp++/config.h>

#include <iterator>
#include <functional>

#if !DSP_BOOST_CONCEPT_CHECKS_DISABLED
#include <boost/concept/requires.hpp>
#include <boost/concept_check.hpp>
#endif // DSP_BOOST_CONCEPT_CHECKS_DISABLED

namespace dsp {

/*!
 * @brief Copy at most n elements from the [src_begin, src_end) sequence into sequence starting at dest.
 * The source and destination sequences must not overlap.
 * @param src_begin start of sequence to copy.
 * @param src_end end of sequence to copy.
 * @param dest the destination of copy operation.
 * @param n maximum number of elements to copy.
 * @return number of copied elements, min(src_end - src_begin, n).
 * @tparam InputIterator type of the src_begin and src_end iterators, which must conform to
 * input iterator concept.
 * @tparam OutputIterator type of the dest iterator, which must conform to output iterator concept
 * and its value type must be convertible to the value type of InputIterator.
 */
template<class InputIterator, class OutputIterator> inline
#if !DSP_BOOST_CONCEPT_CHECKS_DISABLED
BOOST_CONCEPT_REQUIRES(((boost::InputIterator<InputIterator>))
		((boost::OutputIterator<OutputIterator, typename std::iterator_traits<InputIterator>::value_type>)),
		(size_t))
#else
	size_t
#endif
copy_at_most_n(InputIterator src_begin, InputIterator src_end, OutputIterator dest, size_t n)
{
	size_t i;
	for (i = 0; (i < n) && (src_begin != src_end); ++i, ++dest, ++src_begin)
		*dest = *src_begin;
	return i;
}

/*!
 * @brief Copy n elements from the [src, src + n) sequence into [dest, dest + n). The source and destination
 * sequences must not overlap.
 * @param src start of sequence to copy.
 * @param n number of elements to copy.
 * @param dest the destination of copy operation.
 * @return dest + n.
 * @tparam InputIterator type of the src iterator, which must conform to
 * input iterator concept.
 * @tparam OutputIterator type of the dest iterator, which must conform to output iterator concept
 * and its value type must be convertible to the value type of InputIterator.
 */
template<class InputIterator, class OutputIterator> inline
#if !DSP_BOOST_CONCEPT_CHECKS_DISABLED
BOOST_CONCEPT_REQUIRES(((boost::InputIterator<InputIterator>))
		((boost::OutputIterator<OutputIterator, typename std::iterator_traits<InputIterator>::value_type>)),
		(OutputIterator))
#else
	OutputIterator
#endif
copy_n(InputIterator src, size_t n, OutputIterator dest)
{
	for (size_t i = 0; i < n; ++i, ++dest, ++src)
		*dest = *src;
	return dest;
}

//template<class Sample> inline
//void unwrap_transform(Sample* vec, size_t length)
//{
//	Sample*
//}

/*!
 * @brief Base class for simple, sample-based algorithms which work by transforming input samples into
 * output samples in one-to-one manner. Such algorithm will have an <tt>sample_type operator()(sample_type in)</tt>
 * function which takes input sample and returns transformed output sample. This class is a specialization
 * of <tt>std::unary_function</tt> with input and output types being the same - {@link sample_type}.
 * @tparam Sample type of input and output sample.
 */
template<class Sample>
struct sample_based_transform: public std::unary_function<Sample, Sample>
{
	//!@brief Type of input/output sample.
	typedef Sample sample_type;
};

}

#endif /* DSP_ALGORITHM_H_INCLUDED */