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
		SignalData xSignal = generateGaussian(10, 10000);
		SignalData ySignal = generateGaussian(10, 10000);

		const dreal div = divergenceWkv(xSignal, ySignal);
		log() << "Divergence = " << div << logNewLine;
	}

	void addTest()
	{
		timTestList().add("Divergence", testDivergence);
	}

	CallFunction run(addTest);

}
