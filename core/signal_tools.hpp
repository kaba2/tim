#ifndef TIM_SIGNAL_TOOLS_HPP
#define TIM_SIGNAL_TOOLS_HPP

#include "tim/core/signal_tools.h"

#include <pastel/sys/view_tools.h>
#include <pastel/sys/constantiterator.h>

#include <pastel/gfx/drawing.h>

#include <iostream>

namespace Tim
{

	template <typename Signal_ForwardIterator>
	integer minSamples(
		const Signal_ForwardIterator& signalBegin,
		integer signals,
		integer maxLag)
	{
		ENSURE_OP(signals, >=, 0);

		if (signals == 0)
		{
			return 0;
		}

		Signal_ForwardIterator iter = signalBegin;
		integer samples = (*iter)->samples() - maxLag;
		++iter;

		for (integer i = 1;i < samples;++i)
		{
			samples = std::max(samples, (*iter)->samples() - maxLag);

			++iter;
		}

		return samples;
	}

	template <typename Signal_ForwardIterator>
	bool equalDimension(
		const Signal_ForwardIterator& begin,
		integer signals)
	{
		ENSURE_OP(signals, >=, 0);

		if (signals == 0)
		{
			return true;
		}

		Signal_ForwardIterator iter = begin;
		integer dimension = (*begin)->dimension();
		++iter;

		for (integer i = 1;i < signals;++i)
		{
			if ((*iter)->dimension() != dimension)
			{
				return false;
			}

			++iter;
		}

		return true;
	}

	template <
		typename Signal_ForwardIterator,
		typename Lag_ForwardIterator>
	SignalPtr merge(
		const Signal_ForwardIterator& signalBegin,
		integer signals,
		const Lag_ForwardIterator& lagBegin)
	{
		ENSURE_OP(signals, >=, 0);

		if (signals == 0)
		{
			return SignalPtr();
		}

		const integer maxLag = *std::max_element(lagBegin, lagBegin + signals);
		const integer samples = minSamples(signalBegin, signals, maxLag);

		if  (samples <= 0)
		{
			return SignalPtr();
		}

		integer jointDimension = 0;
		Signal_ForwardIterator signalIter = signalBegin;
		for (integer i = 0;i < signals;++i)
		{
			const SignalPtr signal = *signalIter;
			jointDimension += signal->dimension();

			++signalIter;
		}

		// Allocate the joint signal.

		SignalPtr jointSignal(new Signal(samples, jointDimension));
		
		// Copy the signals as parts of the joint signal.

		integer dimensionOffset = 0;

		Lag_ForwardIterator lagIter = lagBegin;

		signalIter = signalBegin;
		for (integer i = 0;i < signals;++i)
		{
			const SignalPtr signal = *signalIter;
			const integer lagOffset = maxLag - *lagIter;
			const integer dimension = signal->dimension();

			for (integer i = 0;i < samples;++i)
			{
				std::copy(
					signal->data().rowBegin(i + lagOffset),
					signal->data().rowEnd(i + lagOffset),
					jointSignal->data().rowBegin(i) + dimensionOffset);
			}

			dimensionOffset += dimension;

			++signalIter;
			++lagIter;
		}

		return jointSignal;
	}

	template <typename Signal_ForwardIterator>
	SignalPtr merge(
		const Signal_ForwardIterator& signalBegin,
		integer signals)
	{
		return Tim::merge(signalBegin, signals,
			ConstantIterator<integer>(0, 0));
	}

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

				drawPixel(Point2(x + 0.5, y * (height - 1)), 
					Color(0, 1, 0), image);
			}
		}
	}

}

#endif
