// Description: Mutual information between marginals of a normal distribution
// DocumentationOf: mutual_information_normal.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/mutual_information_normal.h"

void force_linking_mutual_information_normal() {};

using namespace Tim;

namespace
{

	void matlabMutualInformationNormal(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			ProductOfCovarianceDeterminants,
			JointCovarianceDeterminant,
			Inputs
		};

		enum Output
		{
			MutualInformation,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		dreal productDetCov =
			matlabAsScalar<dreal>(inputSet[ProductOfCovarianceDeterminants]);
		dreal jointDetCov = 
			matlabAsScalar<dreal>(inputSet[JointCovarianceDeterminant]);

		dreal* mutualInformation = 
			matlabCreateScalar<dreal>(outputSet[MutualInformation]);
		*mutualInformation = 
			mutualInformationNormal<dreal>(productDetCov, jointDetCov);
	}

	void addFunction()
	{
		matlabAddFunction(
			"mutual_information_normal",
			matlabMutualInformationNormal);
	}

	CallFunction run(addFunction);

}
