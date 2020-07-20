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
		Signal xSignal = generateGaussian(10000, 10);
		Signal ySignal = generateGaussian(10000, 10);

		const dreal div = divergenceWkv(xSignal, ySignal);
		log() << "Divergence = " << div << logNewLine;
	}

	void addTest()
	{
		timTestList().add("Divergence", testDivergence);
	}

	CallFunction run(addTest);

}
