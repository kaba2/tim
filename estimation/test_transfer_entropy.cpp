#include "estimation.h"

#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/embed.h"

#include <pastel/device/timer.h>

#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

#include <pastel/math/cholesky_decomposition_tools.h>
#include <pastel/math/matrix_tools.h>

using namespace Tim;

namespace
{

	void testTransferEntropy()
	{
		Timer timer;
		timer.setStart();

		log() << "Computing transfer entropies..." << logNewLine;
		log() << "Relative errors to correct analytic results shown in brackets." << logNewLine;

		const integer samples = 10000;
		const integer trials = 50;

		std::vector<SignalPtr> xEnsemble;
		xEnsemble.reserve(trials);

		std::vector<SignalPtr> xFutureEnsemble;
		xFutureEnsemble.reserve(trials);

		std::vector<SignalPtr> yEnsemble;
		yEnsemble.reserve(trials);

		for (integer i = 0;i < trials;++i)
		{
			const SignalPtr xSignal = generateGaussian(samples, 1);
			xEnsemble.push_back(xSignal);

			const SignalPtr xSignalFuture = delayEmbed(xSignal, 1, 1, 1);
			xFutureEnsemble.push_back(xSignalFuture);

			const SignalPtr ySignal = delayEmbed(xSignal, 1, 10, 1);
			yEnsemble.push_back(ySignal);
		}

		std::vector<real> estimateSet;
		transferEntropy(xEnsemble, xFutureEnsemble,
			yEnsemble, Array<2, SignalPtr>(), 5, 1, estimateSet);

		timer.store();
		log() << "Finished in " << timer.seconds() << " seconds." << logNewLine;
	}

	void testAdd()
	{
		timTestList().add("transfer_entropy", testTransferEntropy);
	}

	CallFunction run(testAdd);

}
