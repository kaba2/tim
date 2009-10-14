#include "mex.h"

#include "tim/core/differential_entropy_nk.h"

#include <boost/static_assert.hpp>

#include <pastel/sys/pastelomp.h>

#include <pastel/math/euclidean_normbijection.h>

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		RealIsDouble = boost::is_same<real, double>::value
	};
	BOOST_STATIC_ASSERT(RealIsDouble);

	enum
	{
		xIndex,
		maxRelativeErrorIndex,
		threadsIndex
	};

	const integer trials = mxGetNumberOfElements(inputSet[xIndex]);

	std::vector<SignalPtr> xEnsemble;
	xEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[xIndex], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const integer samples = mxGetN(signalArray);
		const integer dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		xEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	const real maxRelativeError = *mxGetPr(inputSet[maxRelativeErrorIndex]);
	const integer threads = *mxGetPr(inputSet[threadsIndex]);

#if PASTEL_ENABLE_OMP != 0
	omp_set_num_threads(threads);
#endif
	
	integer intrinsicDimension = 0;

	const real entropy = 
		differentialEntropyNk(
		forwardRange(xEnsemble.begin(), xEnsemble.end()), 
		maxRelativeError,
		Euclidean_NormBijection<real>(),
		&intrinsicDimension);

	if (outputs > 0)
	{
		outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* outEntropy = mxGetPr(outputSet[0]);
		*outEntropy = entropy;
	}

	if (outputs > 1)
	{
		outputSet[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* outIntrinsicDimension = mxGetPr(outputSet[1]);
		*outIntrinsicDimension = intrinsicDimension;
	}
}
