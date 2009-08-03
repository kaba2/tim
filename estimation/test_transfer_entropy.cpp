#include "estimation.h"

#include "tim/core/mutual_information.h"
#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/embed.h"

#include <pastel/gfx/drawing.h>
#include <pastel/gfx/pcx.h>

#include <pastel/device/timer.h>

#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

#include <pastel/math/cholesky_decomposition_tools.h>
#include <pastel/math/matrix_tools.h>

using namespace Tim;

namespace
{

	void generateGaussianTest(
		std::vector<SignalPtr>& xEnsemble,
		std::vector<SignalPtr>& xFutureEnsemble,
		std::vector<SignalPtr>& yEnsemble)
	{
		const integer samples = 10000;
		const integer trials = 50;

		for (integer i = 0;i < trials;++i)
		{
			const SignalPtr xSignal = generateGaussian(samples, 1);
			xEnsemble.push_back(xSignal);

			const SignalPtr xSignalFuture = delayEmbed(xSignal, 1, 1, 1);
			xFutureEnsemble.push_back(xSignalFuture);

			const SignalPtr ySignal = delayEmbed(xSignal, 1, 10, 1);
			yEnsemble.push_back(ySignal);
		}
	}

	void generateCouplingTest(
		Array<2, SignalPtr>& signalSet)
	{
		const integer samples = 1500;
		const integer trials = 50;
		const integer yxShift = 5;
		const integer zyShift = 10;
		signalSet.setExtent(trials, 3);

		for (integer i = 0;i < trials;++i)
		{
			const SignalPtr xSignal(new Signal(samples, 1));
			const SignalPtr ySignal(new Signal(samples, 1));
			const SignalPtr zSignal(new Signal(samples, 1));

			generateTimeVaryingCoupling(samples, yxShift, zyShift, 
				xSignal, ySignal, zSignal);

			signalSet(i, 0) = xSignal;
			signalSet(i, 1) = ySignal;
			signalSet(i, 2) = zSignal;

		}
	}

	template <
		typename Signal_ForwardIterator_A,
		typename Signal_ForwardIterator_B,
		typename OutputIterator>
	void mutualInformationLag(
		const Signal_ForwardIterator_A& aBegin,
		const Signal_ForwardIterator_A& aEnd,
		const Signal_ForwardIterator_B& bBegin,
		integer lagBegin,
		integer lagEnd,
		integer kNearest,
		integer maxRelativeError,
		const OutputIterator& outputBegin)
	{
		const integer trials = std::distance(aBegin, aEnd);

		Array<2, SignalPtr> signalSet(trials, 2);
		std::copy(aBegin, aEnd, signalSet.rowBegin(0));
		std::copy(bBegin, bBegin + trials, signalSet.rowBegin(1));

		OutputIterator outputIter = outputBegin;
		std::vector<integer> lagSet(2, 0);
		for (integer i = lagBegin;i < lagEnd;++i)
		{
			lagSet[1] = i;

			(*outputIter) = mutualInformation(
				signalSet, lagSet.begin(),
				kNearest, maxRelativeError);

			++outputIter;
		}
	}

	void testTransferEntropy()
	{
		Timer timer;
		timer.setStart();

		log() << "Computing transfer entropies..." << logNewLine;
		log() << "Relative errors to correct analytic results shown in brackets." << logNewLine;

		Array<2, SignalPtr> signalSet;
		std::vector<real> estimateSet;

		/*
		generateGaussianTest(
			xEnsemble, xFutureEnsemble, yEnsemble);

		transferEntropy(xEnsemble, xFutureEnsemble,
			yEnsemble, zEnsembleSet, 5, 20, estimateSet);

		drawTransferEntropy(estimateSet, "test_mte_gaussian.pcx");
		*/

		generateCouplingTest(signalSet);

		std::vector<real> miSet;
		miSet.reserve(50);

		mutualInformationLag(
			signalSet.rowBegin(0), signalSet.rowEnd(0),
			signalSet.rowBegin(1),
			0, 50,
			1, 0,
			std::back_inserter(miSet));

		const integer signals = signalSet.height();

		Array<2, SignalPtr> futureSet(signalSet.extent());

		for (integer i = 0;i < signals;++i)
		{
			delayEmbedFuture(
				signalSet.rowBegin(i), signalSet.rowEnd(i),
				futureSet.rowBegin(i), 1);
		}

		transferEntropy(
			signalSet,
			0, 2, 
			futureSet.rowBegin(0), futureSet.rowEnd(0),
			5, 20, 
			estimateSet);

		const SignalPtr estimate = SignalPtr(
			new Signal(estimateSet.size(), 1));
		
		std::copy(
			estimateSet.begin(),
			estimateSet.end(),
			estimate->data().begin());

		Array<2, Color> image(estimateSet.size(), 100);

		drawSignal(estimate, arrayView(image));
		savePcx(image, "test_mte_coupling.pcx");

		const SignalPtr miLag = SignalPtr(
			new Signal(miSet.size(), 1));
		
		std::copy(
			miSet.begin(),
			miSet.end(),
			miLag->data().begin());

		drawSignal(miLag, arrayView(image));
		savePcx(image, "test_mi_xy.pcx");

		drawSignal(signalSet(0, 0), arrayView(image));
		savePcx(image, "test_mte_x.pcx");

		drawSignal(signalSet(0, 1), arrayView(image));
		savePcx(image, "test_mte_y.pcx");

		drawSignal(futureSet(0, 0), arrayView(image));
		savePcx(image, "test_mte_xf.pcx");

		timer.store();
		log() << "Finished in " << timer.seconds() << " seconds." << logNewLine;
	}

	void testAdd()
	{
		timTestList().add("transfer_entropy", testTransferEntropy);
	}

	CallFunction run(testAdd);

}
