#include "tim/core/signal_properties.h"

#include <pastel/sys/isnan.h>

namespace Tim
{

	TIM integer nanSamples(const SignalPtr& signal)
	{
		const integer samples = signal->samples();
		const MatrixD& data = signal->data();

		// Count the number of rows which begin
		// with a NaN.

		integer i = 0;
		while(i < samples)
		{
			if (!isNan(data(i, 0)))
			{
				break;
			}
			++i;
		}
	
		return i;
	}

}
