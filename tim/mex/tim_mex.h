#ifndef TIM_TIM_MEX_H
#define TIM_TIM_MEX_H

#include "mex.h"

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/array.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace Tim
{

	namespace Detail
	{

		enum
		{
			RealIsDouble = boost::is_same<real, double>::value
		};
		BOOST_STATIC_ASSERT(RealIsDouble);

	}

	template <typename SignalPtr_OutputIterator>
	void getSignals(const mxArray* input,
					SignalPtr_OutputIterator output);

	template <typename Integer_OutputIterator>
	void getIntegers(const mxArray* input,
					 Integer_OutputIterator output);

	integer getInteger(const mxArray* input);

	real getReal(const mxArray* input);

	void setNumberOfThreads(integer threads);

	void getSignalArray(
		const mxArray* signalSetArray, 
		Array<SignalPtr, 2> signalSet);
	
}

#include "tim/mex/tim_mex.hpp"

#endif
