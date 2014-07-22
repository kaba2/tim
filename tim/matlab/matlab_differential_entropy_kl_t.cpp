// Description: differential_entropy_kl_t
// DocumentationOf: differential_entropy_kl_t.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_kl.h"

void force_linking_differential_entropy_kl_t() {};

using namespace Tim;

namespace
{

	void matlabTemporalDifferentialEntropyKl(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			TimeWindowRadius,
			KNearest,
			FilterIndex,
			Inputs
		};

		enum Output
		{
			Estimate,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		std::vector<Signal> xEnsemble = getSignals(inputSet[X]);

		integer timeWindowRadius = asScalar<integer>(inputSet[TimeWindowRadius]);
		integer kNearest = asScalar<integer>(inputSet[KNearest]);

		std::vector<real> filter;
		getScalars(inputSet[FilterIndex], std::back_inserter(filter));

		Signal estimate = temporalDifferentialEntropyKl(
			range(xEnsemble.begin(), xEnsemble.end()), 
			timeWindowRadius, 
			kNearest,
			Default_NormBijection(),
			range(filter.begin(), filter.end()));

		integer nans = std::max(estimate.t(), (integer)0);
		integer skip = std::max(-estimate.t(), (integer)0); 
		integer samples = std::max(nans + estimate.samples() - skip, (integer)0);

		Array<real> result = createArray<real>(
			Vector2i(samples, 1), outputSet[Estimate]);
		std::fill_n(result.begin(), nans, nan<real>());
		std::copy(estimate.data().begin() + skip, 
			estimate.data().end(), result.begin() + nans);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_kl_t",
			matlabTemporalDifferentialEntropyKl);
	}

	CallFunction run(addFunction);

}
