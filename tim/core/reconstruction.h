// Description: Linear reconstruction of NaNs.

#ifndef TIM_RECONSTRUCTION_H
#define TIM_RECONSTRUCTION_H

#include "tim/core/mytypes.h"
#include "pastel/sys/range.h"

#include <algorithm>

namespace Tim
{

	template <typename Real_Range>
	void reconstruct(Real_Range&& data)
	{
		dreal startValue = (dreal)Nan();

		bool fill = false;
		bool nanRegion = false;
		integer nanRegionSize = 0;
		auto nanBegin = std::begin(data);

		auto iter = std::begin(data);
		auto iterEnd = std::end(data);
		while(iter != iterEnd)
		{

			dreal value = *iter;
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
				dreal endValue = value;
				if (isNan(startValue))
				{
					startValue = endValue;
				}
				if (isNan(endValue))
				{
					endValue = startValue;
				}

				dreal valueAdd = 
					(endValue - startValue) / (nanRegionSize + 1);

				auto nanIter = nanBegin;

				dreal midValue = startValue;
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
