#ifndef TIM_RECONSTRUCTION_HPP
#define TIM_RECONSTRUCTION_HPP

#include "tim/core/reconstruction.h"
#include "tim/core/mytypes.h"

#include <pastel/sys/isnan.h>

#include <algorithm>

namespace Tim
{

	template <typename Real_InputIterator>
	void reconstruct(
		const ForwardIterator_Range<Real_InputIterator>& data)
	{
		real startValue = nan<real>();

		bool fill = false;
		bool nanRegion = false;
		integer nanRegionSize = 0;
		Real_InputIterator nanBegin = data.begin();

		Real_InputIterator iter = data.begin();
		const Real_InputIterator iterEnd = data.end();
		while(iter != iterEnd)
		{
			real value = *iter;
			if (nanRegion)
			{
				if (isNan(value))
				{
					++nanRegionSize;
				}
				else
				{
					nanRegion = false;
					fill = true;
				}
			}
			else
			{
				if (isNan(value))
				{
					nanRegionSize = 1;
					nanBegin = iter;
					nanRegion = true;
				}
				else
				{
					startValue = value;
				}
			}

			if (fill)
			{
				real endValue = value;
				if (isNan(startValue))
				{
					startValue = endValue;
				}
				if (isNan(endValue))
				{
					endValue = startValue;
				}

				const real valueAdd = 
					(endValue - startValue) / (nanRegionSize + 1);

				Real_InputIterator nanIter = nanBegin;

				real midValue = startValue;
				while(nanIter != iter)
				{
					midValue += valueAdd;
					*nanIter = midValue;
					++nanIter;
				}
				
				startValue = endValue;
				fill = false;
			}

			++iter;
		}

		if (nanRegion && !isNan(startValue))
		{
			std::fill(nanBegin, iterEnd, startValue);
		}
	}

}

#endif
