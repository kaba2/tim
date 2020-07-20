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
		void testCase(const boost::iterator_range<Real_Iterator>& data,
			const boost::iterator_range<Real_Iterator2>& correct)
		{
			reconstruct(data);

			TEST_ENSURE(std::equal(data.begin(), data.end(), 
				correct.begin()));

			//std::copy(data.begin(), data.end(), std::ostream_iterator<dreal>(std::cout, " "));
			//std::cout << std::endl;
		}

		void testReconstruction()
		{
			{
				dreal data[] = {1, (dreal)Nan(), (dreal)Nan(), (dreal)Nan(), 5};
				const dreal correct[] = {1, 2, 3, 4, 5};

				testCase(range(data), 
					range(correct));
			}
		
			{
				dreal data[] = {1, (dreal)Nan(), (dreal)Nan(), (dreal)Nan(), (dreal)Nan()};
				const dreal correct[] = {1, 1, 1, 1, 1};

				testCase(range(data), 
					range(correct));
			}

			{
				dreal data[] = {(dreal)Nan(), (dreal)Nan(), (dreal)Nan(), (dreal)Nan(), 5};
				const dreal correct[] = {5, 5, 5, 5, 5};

				testCase(range(data), 
					range(correct));
			}

			{
				dreal data[] = {(dreal)Nan(), (dreal)Nan(), 3, 4, (dreal)Nan()};
				const dreal correct[] = {3, 3, 3, 4, 4};

				testCase(range(data), 
					range(correct));
			}

			{
				dreal data[] = {1, (dreal)Nan(), 3, (dreal)Nan(), 5};
				const dreal correct[] = {1, 2, 3, 4, 5};

				testCase(range(data), 
					range(correct));
			}

			{
				dreal data[] = {(dreal)Nan(), (dreal)Nan(), 3, (dreal)Nan(), (dreal)Nan()};
				const dreal correct[] = {3, 3, 3, 3, 3};

				testCase(range(data), 
					range(correct));
			}

			{
				dreal data[] = {(dreal)Nan(), (dreal)Nan(), 3, 4, 5};
				const dreal correct[] = {3, 3, 3, 4, 5};

				testCase(range(data), 
					range(correct));
			}

			{
				dreal data[] = {1, 2, 3, (dreal)Nan(), (dreal)Nan()};
				const dreal correct[] = {1, 2, 3, 3, 3};

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

