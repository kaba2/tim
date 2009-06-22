#include "tim/core/mutual_information.h"

namespace Tim
{

	TIMCORE real correlatedGaussianMutualInformation(
		real covarianceDeterminant)
	{
		return -0.5 * std::log(covarianceDeterminant);
	}

}

