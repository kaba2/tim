#include "tim/core/mutual_information.h"

#include <pastel/math/cholesky_decomposition_tools.h>

namespace Tim
{

	TIMCORE real correlatedGaussianMutualInformation(
		const DynamicMatrix& correlation)
	{
		const CholeskyDecomposition<Dynamic, real> cholesky(correlation);
		
		return -0.5 * std::log(determinant(cholesky));
	}

}

