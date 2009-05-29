#include "tim/core/embed.h"
#include "tim/core/signal_tools.h"

namespace Tim
{

	TIMCORE SignalPtr delayEmbed(
		const std::vector<real>& data,
		integer dimension,
		integer step)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(step >= 1, step);

		const integer stride = dimension * step;

		const integer points = 
			(data.size() + (stride - 1)) / stride;

		SignalPtr signal = SignalPtr(new Signal(points, dimension));

		SignalView view = signal->view();

		for (integer i = 0;i < points;++i)
		{
			for (integer j = 0;j < dimension;++j)
			{
				view(j, i) = data[(i * dimension + j) * step];
			}
		}

		return signal;
	}

}
