#include "estimation.h"

#include "tim/core/reconstruction.h"

using namespace Tim;

namespace
{

	class ReconstructionTest
		: public TestSuite
	{
	public:
		ReconstructionTest()
			: TestSuite(&timTestReport())
		{
		}

		virtual void run()
		{
			testReconstruction();
		}

		template <typename Real_Iterator, typename Real_Iterator2>
		void testCase(const ForwardIterator_Range<Real_Iterator>& data,
			const ForwardIterator_Range<Real_Iterator2>& correct)
		{
			reconstruct(data);

			TEST_ENSURE(std::equal(data.begin(), data.end(), 
				correct.begin()));

			//std::copy(data.begin(), data.end(), std::ostream_iterator<real>(std::cout, " "));
			//std::cout << std::endl;
		}

		void testReconstruction()
		{
			{
				real data[] = {1, nan<real>(), nan<real>(), nan<real>(), 5};
				const real correct[] = {1, 2, 3, 4, 5};

				testCase(range(data), 
					range(correct));
			}
		
			{
				real data[] = {1, nan<real>(), nan<real>(), nan<real>(), nan<real>()};
				const real correct[] = {1, 1, 1, 1, 1};

				testCase(range(data), 
					range(correct));
			}

			{
				real data[] = {nan<real>(), nan<real>(), nan<real>(), nan<real>(), 5};
				const real correct[] = {5, 5, 5, 5, 5};

				testCase(range(data), 
					range(correct));
			}

			{
				real data[] = {nan<real>(), nan<real>(), 3, 4, nan<real>()};
				const real correct[] = {3, 3, 3, 4, 4};

				testCase(range(data), 
					range(correct));
			}

			{
				real data[] = {1, nan<real>(), 3, nan<real>(), 5};
				const real correct[] = {1, 2, 3, 4, 5};

				testCase(range(data), 
					range(correct));
			}

			{
				real data[] = {nan<real>(), nan<real>(), 3, nan<real>(), nan<real>()};
				const real correct[] = {3, 3, 3, 3, 3};

				testCase(range(data), 
					range(correct));
			}

			{
				real data[] = {nan<real>(), nan<real>(), 3, 4, 5};
				const real correct[] = {3, 3, 3, 4, 5};

				testCase(range(data), 
					range(correct));
			}

			{
				real data[] = {1, 2, 3, nan<real>(), nan<real>()};
				const real correct[] = {1, 2, 3, 3, 3};

				testCase(range(data), 
					range(correct));
			}
		}
	};

	void testReconstruction()
	{
		ReconstructionTest test;
		test.run();
	}

	void addTest()
	{
		timTestList().add("Reconstruction", testReconstruction);
	}

	CallFunction run(addTest);

}

