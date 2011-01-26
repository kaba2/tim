#include "estimation.h"

#include "tim/core/mutual_information.h"
#include "tim/core/partial_transfer_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/embed.h"

#include <pastel/gfx/draw.h>
#include <pastel/gfx/pcx.h>
#include <pastel/gfx/color.h>

#include <pastel/device/timer.h>

#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

#include <pastel/math/cholesky_decomposition_tools.h>
#include <pastel/math/matrix_tools.h>

using namespace Tim;

namespace
{

	template <typename Image_View>
	void drawSignal(
		const SignalPtr& signal,
		const View<2, Color, Image_View>& image)
	{
		const integer dimension = signal->dimension();
		const integer samples = signal->samples();

		const integer height = image.height();

		clear(Color(0), image);

		if (dimension == 1)
		{
			const real yMax = max(abs(signal->data()))[0];
			std::cout << yMax << std::endl;

			for (integer x = 0;x < samples;++x)
			{
				const real y = mabs(signal->data()(x)) / yMax;

				drawPixel(Vector2(x + 0.5, y * (height - 1)), 
					Color(0, 1, 0), image);
			}
		}
	}

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

			const SignalPtr xSignalFuture = delayEmbed(xSignal, 1, 1);
			xFutureEnsemble.push_back(xSignalFuture);

			const SignalPtr ySignal = delayEmbed(xSignal, 1, 1);
			yEnsemble.push_back(ySignal);
		}
	}

	void generateCouplingTest(
		Array<SignalPtr>& signalSet)
	{
		const integer samples = 1500;
		const integer trials = 5;
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

	void testTransferEntropy()
	{
		Timer timer;
		timer.setStart();

		log() << "Computing transfer entropies..." << logNewLine;
		log() << "Relative errors to correct analytic results shown in brackets." << logNewLine;

		const integer kNearest = 20;

		Array<SignalPtr> signalSet;

		//generateGaussianTest(
		//	xEnsemble, xFutureEnsemble, yEnsemble);

		//transferEntropy(xEnsemble, xFutureEnsemble,
		//	yEnsemble, zEnsembleSet, 5, 20, estimateSet);

		//drawTransferEntropy(estimateSet, "test_mte_gaussian.pcx");

		generateCouplingTest(signalSet);

		const integer miLags = 50;

		std::cout << "Computing mutual informations for X <-> Z..." << std::endl;

		std::vector<real> xzMiSet;
		xzMiSet.reserve(miLags);

		for (integer i = 0;i < miLags;++i)
		{
			std::cout << i << " ";

			xzMiSet.push_back(
				mutualInformation(
				range(signalSet.rowBegin(0), signalSet.rowEnd(0)),
				range(signalSet.rowBegin(2), signalSet.rowEnd(2)),
				i, kNearest));
		}
		std::cout << std::endl;

		{
			const SignalPtr miLag = SignalPtr(
				new Signal(xzMiSet.size(), 1));

			Array<Color> image(xzMiSet.size(), 100);
			
			std::copy(
				xzMiSet.begin(),
				xzMiSet.end(),
				miLag->data().begin());

			drawSignal(miLag, arrayView(image));
			savePcx(image, "test_mi_xz.pcx");
		}

		std::vector<SignalPtr> futureSet;

		delayEmbedFuture(
			range(signalSet.rowBegin(0), signalSet.rowEnd(0)),
			std::back_inserter(futureSet), 1);

		const SignalPtr estimate = temporalPartialTransferEntropy(
			range(signalSet.rowBegin(0), signalSet.rowEnd(0)),
			range(signalSet.rowBegin(1), signalSet.rowEnd(1)),
			range(signalSet.rowBegin(2), signalSet.rowEnd(2)),
			range(futureSet.begin(), futureSet.end()),
			5,
			0, 10, 0, 0);

		Array<Color> image(estimate->samples(), 100);

		drawSignal(estimate, arrayView(image));
		savePcx(image, "test_mte_coupling.pcx");

		drawSignal(signalSet(0, 0), arrayView(image));
		savePcx(image, "test_mte_x.pcx");

		drawSignal(signalSet(0, 1), arrayView(image));
		savePcx(image, "test_mte_y.pcx");

		drawSignal(futureSet.front(), arrayView(image));
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
