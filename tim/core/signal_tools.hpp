#ifndef TIM_SIGNAL_TOOLS_HPP
#define TIM_SIGNAL_TOOLS_HPP

#include "tim/core/signal_tools.h"

#include <pastel/sys/view_tools.h>
#include <pastel/sys/constantiterator.h>

#include <pastel/gfx/draw.h>

#include <iostream>

namespace Tim
{

	template <typename Image_View>
	void drawSignal(
		const SignalPtr& signal,
		const View<2, Color, Image_View>& image)
	{
		const integer dimension = signal->dimension();
		const integer samples = signal->samples();

		const integer width = image.width();
		const integer height = image.height();

		clear(Color(0), image);

		if (dimension == 1)
		{
			const real yMax = max(abs(signal->data()))[0];
			std::cout << yMax << std::endl;

			for (integer x = 0;x < samples;++x)
			{
				const real y = mabs(signal->data()(x)) / yMax;

				drawPixel(Vector2(x + 0.5, y * (height - 1)), 
					Color(0, 1, 0), image);
			}
		}
	}

	template <typename SignalPtr_Iterator>
	void constructPointSet(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
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
