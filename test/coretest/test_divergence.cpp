#include "estimation.h"

#include "tim/core/signal_tools.h"
#include "tim/core/delay_embed.h"

#include <pastel/sys/random.h>

#include "tim/core/divergence_wkv.h"

using namespace Tim;

namespace
{

	void testDivergence()
	{
		SignalPtr xSignal = generateGaussian(10000, 10);
		SignalPtr ySignal = generateGaussian(10000, 10);

		const real div = divergenceWkv(xSignal, ySignal);
		log() << "Divergence = " << div << logNewLine;
	}

	void addTest()
	{
		timTestList().add("Divergence", testDivergence);
	}

	CallFunction run(addTest);

}
