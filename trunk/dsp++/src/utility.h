/*!
 * @author Andrzej Ciarkowski <mailto:andrzej.ciarkowski@gmail.com>
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <cstddef>
#include <functional>

#define DSP_ARRAY_SIZE(arr) (sizeof(arr)/sizeof(*arr))

namespace dsp { namespace detail {

template<class Object, class Member, class Test>
struct member_equal_t: public std::unary_function<Object, bool>
{
	Member Object::* member;
	Test test;

	bool operator ()(const Object& obj) const {return test == obj.*member;}

	member_equal_t(Member Object::* m, const Test& t): member(m), test(t) {}
};

template<class Object>
struct member_equal_t<Object, const char*, const char*>: public std::unary_function<Object, bool>
{
	const char* Object::* member;
	const char* test;

	bool operator ()(const Object& obj) const {return 0 == strcmp(test, obj.*member);}

	member_equal_t(const char* Object::* m, const char* t): member(m), test(t) {}
};


template<class Object, class Member, class Test>
member_equal_t<Object, Member, Test> member_equal(Member Object::* member, const Test& test)
{return member_equal_t<Object, Member, Test>(member, test);}

template<class Element, size_t N, class Member, class Test>
const Element* const match_element(const Element (& array)[N], Member Element::* member, const Test& test)
{
	const Element* end = array + N;
	const Element* res = std::find_if(array, end, member_equal_t<Element, Member, Test>(member, test));
	return (res == end ? NULL : res);
}

}}

#endif /* UTILITY_H_ */
