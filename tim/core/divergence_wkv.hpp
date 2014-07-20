#ifndef TIM_DIVERGENCE_WKV_HPP
#define TIM_DIVERGENCE_WKV_HPP

#include "tim/core/divergence_wkv.h"
#include "tim/core/signalpointset.h"

#include <pastel/geometry/pointkdtree_search_nearest.h>
#include <pastel/sys/predicate_indicator.h>
#include <pastel/sys/notequalto.h>

namespace Tim
{

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator>
	real divergenceWkv(
		const boost::iterator_range<SignalPtr_X_Iterator>& xSignalSet,
		const boost::iterator_range<SignalPtr_Y_Iterator>& ySignalSet)
	{
		// "A Nearest-Neighbor Approach to Estimating
		// Divergence between Continuous Random Vectors"
		// Qing Wang, Sanjeev R. Kulkarni, Sergio Verdu,
		// IEEE International Symposium on Information Theory (ISIT), 
		// 2006.

		if (xSignalSet.empty() || ySignalSet.empty())
		{
			return 0;
		}

		const integer xDimension = xSignalSet.front()->dimension();
		const integer yDimension = ySignalSet.front()->dimension();

		ENSURE_OP(xDimension, ==, yDimension);

		typedef SignalPointSet::Point_ConstIterator Point_ConstIterator;

		// Construct point-sets.

		SignalPointSetPtr xPointSet(new SignalPointSet(xSignalSet));
		SignalPointSetPtr yPointSet(new SignalPointSet(ySignalSet));

		const integer xSamples = xPointSet->samples();
		const integer ySamples = yPointSet->samples();

		integer acceptedSamples = 0;
		real estimate = 0;
#pragma omp parallel for reduction(+ : estimate, acceptedSamples)
		for (integer i = 0;i < xSamples;++i)
		{
			// Find out the nearest neighbor in X for a point in X.

			const Point_ConstIterator query = 
				*(xPointSet->begin() + i);

			const real xxDistance2 = 
				searchNearest(xPointSet->kdTree(), query,
					nullOutput(),
					predicateIndicator(query, NotEqualTo()));
			
			if (xxDistance2 > 0 && xxDistance2 < infinity<real>())
			{		
				// Find out the nearest neighbor in Y for a point in X.

				const Vector<real> queryPoint(
					ofDimension(xDimension), 
					withAliasing((real*)(query)->point()));

				const real xyDistance2 = 
					searchNearest(yPointSet->kdTree(), queryPoint);
				
				if (xyDistance2 > 0 && xyDistance2 < infinity<real>())
				{
					estimate += std::log(xyDistance2 / xxDistance2);
					++acceptedSamples;
				}
			}
		}

		if (acceptedSamples > 0)
		{
			// The factor 2 in the denominator is because 
			// 'xyDistance' and 'xxDistance' are squared distances
			// and thus need to be taken a square root. However,
			// this can be taken outside the logarithm with a 
			// division by 2.
			estimate *= (real)xDimension / (2 * acceptedSamples);
			estimate += std::log((real)ySamples / (xSamples - 1));
		}
		else
		{
			estimate = nan<real>();
		}

		return estimate;
	}

}

#endif
