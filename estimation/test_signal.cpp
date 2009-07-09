#include "estimation.h"

#include "tim/core/signal_tools.h"
#include "tim/core/mutual_information.h"
#include "tim/core/embed.h"

#include <pastel/sys/view_all.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

using namespace Tim;

namespace
{

	class SignalTest
		: public TestSuite
	{
	public:
		SignalTest()
			: TestSuite(&timTestReport())
		{
		}

		virtual void run()
		{
			testEmbed();
			testSplit();
		}

		void testSplit()
		{
			SignalPtr xy = 
				SignalPtr(new Signal(3, 5));

			xy->data() |=
				0, 1, 2, 3, 4,
				5, 6, 7, 8, 9,
				10, 11, 12, 13, 14;

			SignalPtr x = slice(xy, 0, 1);
			SignalPtr xCorrect =
				SignalPtr(new Signal(1, 5));
			xCorrect->data() |=
				0, 1, 2, 3, 4;

			TEST_ENSURE(x->data() == xCorrect->data());

			SignalPtr y = slice(xy, 1, 3);
			SignalPtr yCorrect =
				SignalPtr(new Signal(2, 5));
			yCorrect->data() |=
				5, 6, 7, 8, 9,
				10, 11, 12, 13, 14;

			TEST_ENSURE(y->data() == yCorrect->data());
		}

		void testEmbed()
		{
			SignalPtr x = 
				SignalPtr(new Signal(1, 15));
			x->data() |= 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;

			{
				SignalPtr y = delayEmbed(x, 3);
				SignalPtr yCorrect =
					SignalPtr(new Signal(3, 13));
				yCorrect->data() |=
					0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
					1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
					2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;

				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 1);
				SignalPtr yCorrect =
					SignalPtr(new Signal(3, 12));
				yCorrect->data() |=
					1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
					2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
					3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;

				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 2);
				SignalPtr yCorrect =
					SignalPtr(new Signal(3, 11));
				yCorrect->data() |=
					2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
					3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
					4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;
		
				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 0, 2);
				SignalPtr yCorrect =
					SignalPtr(new Signal(3, 11));
				yCorrect->data() |=
					0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
					2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
					4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;
				/*
				yCorrect->data() |=
					0, 2, 4, 6, 8, 10,
					2, 4, 6, 8, 10, 12,
					4, 6, 8, 10, 12, 14;
				*/
		
				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 1, 2);
				SignalPtr yCorrect =
					SignalPtr(new Signal(3, 10));
				yCorrect->data() |=
					1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
					3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
					5, 6, 7, 8, 9, 10, 11, 12, 13, 14;
		
				TEST_ENSURE(y->data() == yCorrect->data());
			}
		}
	};

	void testSignal()
	{
		SignalTest test;
		test.run();
	}

	void addTest()
	{
		timTestList().add("signal", testSignal);
	}

	CallFunction run(addTest);

}

