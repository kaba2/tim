// Description: Algorithms for Signal's

#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"
#include "tim/core/signal_merge.h"
#include "tim/core/signal_split.h"
#include "tim/core/signal_properties.h"

#include <pastel/math/matrix.h>

#include <pastel/gfx/color.h>

#include <pastel/sys/iterator_range.h>
#include <pastel/sys/array.h>

#include <iostream>

namespace Tim
{

	//! Converts from nan-padded form to lagged form.
	/*!
	The output signal aliases the input signal such that
	the NaNs from the beginning of the signal are skipped.

	TIM Matlab allows two forms for signals, a _lagged form_, 
	and a _NaN-padded form_. In NaN-padded form the start of the 
	signal contains NaNs, whose role is only to keep the first
	sample associated to the time instant zero. When converting 
	from NaN-padded form to lagged form, the NaN part is removed 
	from the signal and the signal is interpreted to begin at
	the first non-NaN sample. To maintain time-correspondence,
	the first sample in the lagged form is taken to begin at a 
	time instant given by the number of NaN samples in the 
	NaN-padded form.
	*/
	TIM SignalPtr nanToLagged(const SignalPtr& signal);

	//! Prints the signal into an output stream.
	TIM std::ostream& operator<<(
		std::ostream& stream, const Signal& signal);

	//! Computes the covariance of the signal samples.
	TIM void computeCovariance(
		const SignalPtr& signal,
		MatrixD& result);

	//! Transforms the given signal to identity covariance.
	TIM void normalizeCovariance(
		const SignalPtr& signal,
		const MatrixD& covariance);

	template <typename Image_View>
	void drawSignal(
		const SignalPtr& signal,
		const View<2, Color, Image_View>& image);

	template <typename SignalPtr_Iterator>
	void constructPointSet(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer sampleBegin,
		integer sampleEnd,
		integer dimensionBegin,
		integer dimensionEnd,
		std::vector<const real*>& pointSet);

}

#include "tim/core/signal_generate.h"

#include "tim/core/signal_tools.hpp"

#endif
