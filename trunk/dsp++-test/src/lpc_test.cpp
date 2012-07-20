/*!
 * @file lpc_test.cpp
 * 
 * @author Andrzej Ciarkowski <mailto:andrzej.ciarkowski@gmail.com>
 */

#include "lpc_test.h"
#include <dsp++/lpc.h>
#include <dsp++/float.h>

const float x[] = {0.5376671552658081,1.833885073661804,-2.258846759796143,0.862173318862915,0.318765252828598,-1.307688355445862,-0.4335920214653015};
const float a[] = {1,	0.532430375671096,	0.250609258485180,	-0.186377797109935,	0.0145011350307570,	0.192908951389234,	0.256205216566800};
const float e = 1.07965299964618;

const std::complex<double> cx[] = {0.342624466538650 + 0.714742903826096i, 3.57839693972576 - 0.204966058299775i,
		2.76943702988488 - 0.124144348216312i, -1.34988694015652 + 1.48969760778547i,
		3.03492346633185 + 1.40903448980048i, 0.725404224946106 + 1.41719241342961i,
		-0.0630548731896562 + 0.671497133608081i};

const std::complex<double> ca[] = {1.00000000000000 + 0.00000000000000i,	-0.160556055859024 + 0.0639303674052243i,
       -0.0593837867213146 - 0.137844464516607i, -0.395560957265486 - 0.272596618061747i,
       0.0321255729170954 - 0.0763324603694148i, -0.00768984884278122 + 0.0915604090332977i,
       0.0352739459789746 + 0.207142520125894i};

const std::complex<double> ce = 4.00753693110963;


void dsp::test::lpc_test::test_lpc_real()
{
	const size_t L = 7;
	dsp::lpc<double> lpc(L);
	float aa[L];
	float e = lpc(x, aa);
	CPPUNIT_ASSERT(std::equal(aa, aa + L - 1, a, dsp::within_range<float>(0.00001f)));
	CPPUNIT_ASSERT(dsp::within_range<float>(0.00001f)(e, ::e));
}

void dsp::test::lpc_test::test_lpc_complex()
{
	const size_t L = 7;
	dsp::lpc<std::complex<double> > lpc(L);
	std::complex<double> aa[L];
	std::complex<double> e = lpc(cx, aa);
	CPPUNIT_ASSERT(std::equal(aa, aa + L - 1, ca, dsp::within_range<std::complex<double> >(0.00001)));
	CPPUNIT_ASSERT(dsp::within_range<std::complex<double> >(0.00001)(e, ::ce));
}

