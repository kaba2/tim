#include "tim/core/mutual_information.h"

#include <pastel/geometry/all_nearest_neighbors.h>

namespace Tim
{


	TIMCORE real mutualInformation(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal)
	{
		std::vector<SignalPtr> signalSet;
		signalSet.push_back(aSignal);
		signalSet.push_back(bSignal);
		
		return Pastel::mutualInformation(
			signalSet);
	}

}
