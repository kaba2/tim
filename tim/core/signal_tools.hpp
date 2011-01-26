#ifndef TIM_SIGNAL_TOOLS_HPP
#define TIM_SIGNAL_TOOLS_HPP

#include "tim/core/signal_tools.h"

#include <pastel/sys/view_tools.h>
#include <pastel/sys/constantiterator.h>

#include <pastel/gfx/draw.h>

#include <iostream>

namespace Tim
{

	template <typename SignalPtr_Iterator>
	void constructPointSet(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer sampleBegin,
		integer sampleEnd,
		integer dimensionBegin,
		integer dimensionEnd,
		std::vector<const real*>& pointSet)
	{
		if (signalSet.empty())
		{
			pointSet.clear();
			return;
		}

		ENSURE(equalDimension(signalSet));

		const integer signalSamples = minSamples(signalSet);
		const integer signalDimension = signalSet.front()->dimension();

		ENSURE_OP(sampleBegin, <=, sampleEnd);
		ENSURE_OP(sampleBegin, >= , 0);
		ENSURE_OP(sampleEnd, <=, signalSamples);
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signalDimension);

		const integer dimension = dimensionEnd - dimensionBegin;
		const integer samples = sampleEnd - sampleBegin;
		const integer trials = signalSet.size();

		if (dimension == 0 ||
			samples == 0)
		{
			pointSet.clear();
			return;
		}

		pointSet.resize(samples * trials);

		SignalPtr_Iterator iter = signalSet.begin();
		const SignalPtr_Iterator iterEnd = signalSet.end();
		integer trial = 0;
		while(iter != iterEnd)
		{
			const SignalPtr signal = *iter;

			for (integer i = sampleBegin;i < sampleEnd;++i)
			{
				// The samples from the trials are interleaved.

				const integer index = (i - sampleBegin) * trials + trial;
				pointSet[index] = &signal->data()(i, dimensionBegin);
			}

			++trial;
			++iter;
		}
	}

}

#endif
