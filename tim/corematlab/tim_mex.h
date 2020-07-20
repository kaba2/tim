// Description: Functions to ease interfacing with Matlab
// Documentation: tim_matlab_impl.txt

#ifndef TIM_TIM_MEX_H
#define TIM_TIM_MEX_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>
#include <pastel/sys/sequence/copy_n.h>
#include <pastelmatlab/matlab_argument.h>

#include <vector>

namespace Tim
{

	inline Signal asSignal(const MatrixView<dreal>& matrix) {
		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.

		integer samples = matrix.cols();
		integer dimension = matrix.rows();
		
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

		integer nans = 0;
		while(nans < samples)
		{
			if (!isNan(matrix(0, nans)))
			{
				break;
			}
			++nans;
		}

		return Signal(matrix.slicex(nans), nans);
	}

	template <ranges::forward_range Range>
	ranges::forward_range auto matlabMatricesAsSignals(const Range& matrices) 
	{
		return matrices | ranges::views::transform(
			[](auto&& matrix) -> Signal {return asSignal(matrix.view());}
		);
	}

	//! Retrieves an array of signals.
    template <typename Real>
	inline Array<Signal> asSignalArray(
		const Array<MatlabMatrix<Real>> matrices)
	{
		integer trials = matrices.width();
		integer signals = matrices.height();
		
		Array<Signal> signalSet(Vector2i(trials, signals));
		for (integer y = 0;y < signals;++y)
		{
			for (integer x = 0;x < trials;++x)
			{
				signalSet(x, y) = asSignal(matrices(x, y).view());
			}
		}

		return signalSet;
	}


}

#endif
