#include "tim/console/functions.h"
#include "tim/console/console_parser.h"
#include "tim/console/errorlog.h"

#include <tim/core/differential_entropy.h>
#include <tim/core/mutual_information.h>
#include <tim/core/partial_mutual_information.h>
#include <tim/core/transfer_entropy.h>
#include <tim/core/partial_transfer_entropy.h>
#include <tim/core/divergence_wkv.h>

#include <pastel/sys/string_tools.h>

#include <fstream>
#include <map>
#include <string>

#include <boost/function.hpp>

namespace Tim
{

	namespace
	{

		boost::any read_csv(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			const std::string fileName = boost::any_cast<std::string>(argSet[0]);
			integer samples = boost::any_cast<integer>(argSet[1]);
			const integer dimension = boost::any_cast<integer>(argSet[2]);
			const integer trials = boost::any_cast<integer>(argSet[3]);
			const integer series = boost::any_cast<integer>(argSet[4]);
			const std::string separatorSet = boost::any_cast<std::string>(argSet[5]);
			const SignalPtr order = boost::any_cast<SignalPtr>(argSet[6]);

			// Check parameters.

			bool error = false;
			if (samples < 0)
			{
				reportError("'samples' must be non-negative.");
				error = true;
			}
			if (dimension < 0)
			{
				reportError("'dimension' must be non-negative.");
				error = true;
			}
			if (trials < 0)
			{
				reportError("'trials' must be non-negative.");
				error = true;
			}
			if (series < 0)
			{
				reportError("'series' must be non-negative.");
				error = true;
			}

			const Tuple<integer, 4> permutation(
				order->data()(0),
				order->data()(1),
				order->data()(2),
				order->data()(3));

			Tuple<integer, 4> check(permutation);
			std::sort(check.begin(), check.end());
			if (check[0] != 0 || check[1] != 1 ||
				check[2] != 2 || check[3] != 3)
			{
				reportError("'order' must contain exactly the integers 0, 1, 2, and 3 in some order.");
				throw FunctionCall_Exception();
			}

			std::ifstream file(fileName.c_str());
			if (!file.is_open())
			{
				reportError("Could not open data file " + fileName + ".");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			const bool samplesUnknown = (passedArgs == 1);

			if (!samplesUnknown && 
				(samples == 0 || dimension == 0 || series == 0 || trials == 0))
			{
				// Nothing to read.
				return boost::any(SignalPtr(new Signal));
			}

			Vector<integer, 4> width = permute(
				Vector<integer, 4>(samples, dimension, trials, series),
				permutation);

			// Read the values into a continuous array.

			std::vector<real> data;
			real value = 0;

			const integer samplesToRead = product(width);
			if (!samplesUnknown)
			{
				data.reserve(samplesToRead);
			}

			integer readSamples = 0;

			while(true)
			{
				// Read as many white-space separated values
				// as possible.
				while(file >> value)
				{
					data.push_back(value);
					++readSamples;
				}
				if (!samplesUnknown && readSamples == samplesToRead)
				{
					// No need to read more samples since
					// we already got all we need to form the
					// signals.
					break;
				}
				
				// Clear the failbit, since we just failed
				// to read a floating point number.
				file.clear();

				// See if we can read a separator instead.
				char separator = 0;
				file >> separator;
				if (!file)
				{
					// No luck. Seems like this is the
					// end of the data.
					break;
				}
				else if (separatorSet.find(separator) == std::string::npos)
				{
					// There was something, but it was not a separator either.
					reportError("Data file " + fileName + " had an invalid separator " +
						separator + ".");
					throw FunctionCall_Exception();
				}
			}

			// Unpack the flattened data to separate signals.

			if (samplesUnknown)
			{
				// We now know the number of samples.
				samples = data.size();
				width[permutation[0]] = samples;
			}

			Vector<integer, 4> stride;
			for (integer i = 0;i < 4;++i)
			{
				const integer position = permutation[i];
				integer accu = 1;
				for (integer j = 0;j < position;++j)
				{
					accu *= width[j];
				}
				stride[i] = accu;
			}

			Array<SignalPtr>* cellArray = new Array<SignalPtr>(trials, series);

			for (integer y = 0;y < series;++y)
			{
				for (integer x = 0;x < trials;++x)
				{
					SignalPtr signal = SignalPtr(new Signal(samples, dimension));
					for (integer i = 0;i < dimension;++i)
					{
						for (integer j = 0;j < samples;++j)
						{
							const integer offset = dot(stride, Vector<integer, 4>(j, i, x, y));
							signal->data()(j, i) = data[offset];
						}
					}
					(*cellArray)(x, y) = signal;
				}
			}

			
			return boost::any(cellArray);
		}

		boost::any differential_entropy_kl(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* cell = boost::any_cast<Cell*>(argSet[0]);
			real maxRelativeError = boost::any_cast<real>(argSet[1]);
			integer kNearest = boost::any_cast<integer>(argSet[2]);

			// Check parameters.

			bool error = false;
			if (maxRelativeError < 0)
			{
				reportError("maxRelativeError must be non-negative.");
				error = true;
			}
			
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.

			const real de = differentialEntropyKl(
				forwardRange(cell->begin(), cell->end()),
				maxRelativeError, kNearest);
				
			return boost::any(de);
		}
		
		boost::any differential_entropy_kl_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* cell = boost::any_cast<Cell*>(argSet[0]);
			integer timeWindowRadius = boost::any_cast<integer>(argSet[1]);;
			real maxRelativeError = boost::any_cast<real>(argSet[2]);
			integer kNearest = boost::any_cast<integer>(argSet[3]);
			
			// Check parameters.

			bool error = false;
			if (timeWindowRadius < 0)
			{
				reportError("timeWindowRadius must be non-negative.");
				error = true;
			}
			
			if (maxRelativeError < 0)
			{
				reportError("maxRelativeError must be non-negative.");
				error = true;
			}
			
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			std::vector<real> deSet;
			
			temporalDifferentialEntropyKl(
				forwardRange(cell->begin(), cell->end()),
				timeWindowRadius,
				std::back_inserter(deSet),
				maxRelativeError, kNearest);
				
			SignalPtr signal = SignalPtr(new Signal(deSet.size(), 1));
			std::copy(deSet.begin(), deSet.end(),
				signal->data().begin());
				
			return boost::any(signal);
		}

		boost::any differential_entropy_nk(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* cell = boost::any_cast<Cell*>(argSet[0]);;
			real maxRelativeError = boost::any_cast<real>(argSet[1]);

			// Check parameters.

			bool error = false;
			if (maxRelativeError < 0)
			{
				reportError("maxRelativeError must be non-negative.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			integer dimension = 0;

			const real de = differentialEntropyNk(
				forwardRange(cell->begin(), cell->end()),
				maxRelativeError, 
				Default_NormBijection(),
				&dimension);
				
			SignalPtr signal = SignalPtr(new Signal(2, 1));
			signal->data()(0) = de;
			signal->data()(1) = dimension;
				
			return boost::any(signal);
		}

		boost::any mutual_information_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			integer timeWindowRadius = boost::any_cast<integer>(argSet[2]);
			integer xLag = boost::any_cast<integer>(argSet[3]);
			integer yLag = boost::any_cast<integer>(argSet[4]);
			integer kNearest = boost::any_cast<integer>(argSet[5]);

			// Check parameters.
			
			bool error = false;
			if (timeWindowRadius < 0)
			{
				reportError("timeWindowRadius must be non-negative.");
				error = true;
			}
			
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			std::vector<real> miSet;
			
			temporalMutualInformation(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				timeWindowRadius,
				std::back_inserter(miSet),
				xLag, yLag,
				kNearest);
				
			SignalPtr signal = SignalPtr(new Signal(miSet.size(), 1));
			std::copy(miSet.begin(), miSet.end(),
				signal->data().begin());
				
			return boost::any(signal);
		}

		boost::any mutual_information(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			integer xLag = boost::any_cast<integer>(argSet[2]);
			integer yLag = boost::any_cast<integer>(argSet[3]);
			integer kNearest = boost::any_cast<integer>(argSet[4]);

			// Check parameters.

			bool error = false;
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.

			const real mi = mutualInformation(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				xLag, yLag,
				kNearest);
				
			return boost::any(mi);
		}

		boost::any mutual_information_pt(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			Cell* zCell = boost::any_cast<Cell*>(argSet[2]);
			integer timeWindowRadius = boost::any_cast<integer>(argSet[3]);
			integer xLag = boost::any_cast<integer>(argSet[4]);
			integer yLag = boost::any_cast<integer>(argSet[5]);
			integer zLag = boost::any_cast<integer>(argSet[6]);
			integer kNearest = boost::any_cast<integer>(argSet[7]);
			
			// Check parameters.

			bool error = false;
			if (timeWindowRadius < 0)
			{
				reportError("timeWindowRadius must be non-negative.");
				error = true;
			}
			
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			std::vector<real> miSet;
			
			temporalPartialMutualInformation(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				forwardRange(zCell->begin(), zCell->end()),
				timeWindowRadius,
				std::back_inserter(miSet),
				xLag, yLag, zLag,
				kNearest);
				
			SignalPtr signal = SignalPtr(new Signal(miSet.size(), 1));
			std::copy(miSet.begin(), miSet.end(),
				signal->data().begin());
				
			return boost::any(signal);
		}

		boost::any mutual_information_p(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			Cell* zCell = boost::any_cast<Cell*>(argSet[2]);
			integer xLag = boost::any_cast<integer>(argSet[3]);
			integer yLag = boost::any_cast<integer>(argSet[4]);
			integer zLag = boost::any_cast<integer>(argSet[5]);
			integer kNearest = boost::any_cast<integer>(argSet[6]);
			
			// Check parameters.

			bool error = false;
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const real pmi = partialMutualInformation(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				forwardRange(zCell->begin(), zCell->end()),
				xLag, yLag, zLag,
				kNearest);
				
			return boost::any(pmi);
		}

		boost::any transfer_entropy_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			Cell* wCell = boost::any_cast<Cell*>(argSet[2]);
			integer timeWindowRadius = boost::any_cast<integer>(argSet[3]);
			integer xLag = boost::any_cast<integer>(argSet[4]);
			integer yLag = boost::any_cast<integer>(argSet[5]);
			integer wLag = boost::any_cast<integer>(argSet[6]);
			integer kNearest = boost::any_cast<integer>(argSet[7]);
			
			// Check parameters.

			bool error = false;
			if (timeWindowRadius < 0)
			{
				reportError("timeWindowRadius must be non-negative.");
				error = true;
			}
			
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			std::vector<real> teSet;
			
			temporalTransferEntropy(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				forwardRange(wCell->begin(), wCell->end()),
				timeWindowRadius,
				std::back_inserter(teSet),
				xLag, yLag, wLag,
				kNearest);
				
			SignalPtr signal = SignalPtr(new Signal(teSet.size(), 1));
			std::copy(teSet.begin(), teSet.end(),
				signal->data().begin());
				
			return boost::any(signal);
		}

		boost::any transfer_entropy_pt(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			Cell* zCell = boost::any_cast<Cell*>(argSet[2]);
			Cell* wCell = boost::any_cast<Cell*>(argSet[3]);
			integer timeWindowRadius = boost::any_cast<integer>(argSet[4]);
			integer xLag = boost::any_cast<integer>(argSet[5]);
			integer yLag = boost::any_cast<integer>(argSet[6]);
			integer zLag = boost::any_cast<integer>(argSet[7]);
			integer wLag = boost::any_cast<integer>(argSet[8]);
			integer kNearest = boost::any_cast<integer>(argSet[9]);

			// Check parameters.

			bool error = false;
			if (timeWindowRadius < 0)
			{
				reportError("timeWindowRadius must be non-negative.");
				error = true;
			}
			
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			std::vector<real> teSet;
			
			temporalPartialTransferEntropy(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				forwardRange(zCell->begin(), zCell->end()),
				forwardRange(wCell->begin(), wCell->end()),
				timeWindowRadius,
				std::back_inserter(teSet),
				xLag, yLag, zLag, wLag,
				kNearest);
				
			SignalPtr signal = SignalPtr(new Signal(teSet.size(), 1));
			std::copy(teSet.begin(), teSet.end(),
				signal->data().begin());
				
			return boost::any(signal);
		}

		boost::any transfer_entropy(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			Cell* wCell = boost::any_cast<Cell*>(argSet[2]);
			integer xLag = boost::any_cast<integer>(argSet[3]);
			integer yLag = boost::any_cast<integer>(argSet[4]);
			integer wLag = boost::any_cast<integer>(argSet[5]);
			integer kNearest = boost::any_cast<integer>(argSet[6]);

			// Check parameters.

			bool error = false;
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const real te = transferEntropy(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				forwardRange(wCell->begin(), wCell->end()),
				xLag, yLag, wLag,
				kNearest);
				
			return boost::any(te);
		}

		boost::any transfer_entropy_p(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.
			
			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);
			Cell* zCell = boost::any_cast<Cell*>(argSet[2]);
			Cell* wCell = boost::any_cast<Cell*>(argSet[3]);
			integer xLag = boost::any_cast<integer>(argSet[4]);
			integer yLag = boost::any_cast<integer>(argSet[5]);
			integer zLag = boost::any_cast<integer>(argSet[6]);
			integer wLag = boost::any_cast<integer>(argSet[7]);
			integer kNearest = boost::any_cast<integer>(argSet[8]);

			// Check parameters.

			bool error = false;
			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const real pte = partialTransferEntropy(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()),
				forwardRange(zCell->begin(), zCell->end()),
				forwardRange(wCell->begin(), wCell->end()),
				xLag, yLag, zLag, wLag,
				kNearest);
				
			return boost::any(pte);
		}

		boost::any divergence_wkv(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			Cell* xCell = boost::any_cast<Cell*>(argSet[0]);
			Cell* yCell = boost::any_cast<Cell*>(argSet[1]);

			// Compute.

			const real div = divergenceWkv(
				forwardRange(xCell->begin(), xCell->end()),
				forwardRange(yCell->begin(), yCell->end()));
				
			return boost::any(div);
		}

		struct FunctionInfo
		{
			typedef boost::function<boost::any(const AnySet& argSet, integer passedArgs)> Callback;

			FunctionInfo()
			: callback()
			, minArgs(0)
			, parameterSet()
			{
			}

			template <typename Any_Iterator>
			FunctionInfo(
				Callback callback_, 
				const ForwardRange<Any_Iterator>& parameterSet_,
				integer minArgs_)
			: callback(callback_)
			, parameterSet(parameterSet_.begin(), parameterSet_.end())
			, minArgs(minArgs_)
			{
				ENSURE_OP(minArgs, <=, parameterSet.size());
			}

			Callback callback;
			std::vector<boost::any> parameterSet;
			integer minArgs;
		};

		typedef std::map<std::string, FunctionInfo> FunctionMap;
		typedef FunctionMap::const_iterator FunctionIterator;

		FunctionMap functionMap;

		void initializeFunctions()
		{
			static bool initialized = false;

			if (initialized)
			{
				return;
			}

			// read_csv
			{
				SignalPtr order = SignalPtr(new Signal(4, 1));
				order->data()(0) = 0;
				order->data()(1) = 1;
				order->data()(2) = 2;
				order->data()(3) = 3;

				boost::any parameterSet[] =
				{
					// filename
					boost::any(std::string()),
					// samples
					boost::any(integer(0)),
					// dimension
					boost::any(integer(1)),
					// trials
					boost::any(integer(1)),
					// series
					boost::any(integer(1)),
					// separatorSet
					boost::any(std::string(",;")),
					// order
					boost::any(order)
				};
			
				functionMap.insert(
					std::make_pair(
					"read_csv", 
					FunctionInfo(
						read_csv, 
						forwardRange(parameterSet), 1)));
			}

			// differential_entropy_kl
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any((Cell*)0),
					// maxRelativeError
					boost::any((real)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"differential_entropy_kl", 
					FunctionInfo(
						differential_entropy_kl, 
						forwardRange(parameterSet), 1)));
			}
			
			// differential_entropy_kl_t
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any((Cell*)0),
					// timeWindowRadius
					boost::any((integer)0),
					// maxRelativeError
					boost::any((real)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"differential_entropy_kl_t", 
					FunctionInfo(
						differential_entropy_kl_t, 
						forwardRange(parameterSet), 2)));
			}
				
			// differential_entropy_nk
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any((Cell*)0),
					// maxRelativeError
					boost::any((real)0)
				};
			
				functionMap.insert(
					std::make_pair(
					"differential_entropy_nk", 
					FunctionInfo(
						differential_entropy_nk, 
						forwardRange(parameterSet), 1)));
			}

			// divergence_wkv
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0)
				};
			
				functionMap.insert(
					std::make_pair(
					"divergence_wkv", 
					FunctionInfo(
						divergence_wkv, 
						forwardRange(parameterSet), 2)));
			}

			// mutual_information_t
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// timeWindowRadius
					boost::any((integer)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"mutual_information_t", 
					FunctionInfo(
						mutual_information_t, 
						forwardRange(parameterSet), 3)));
			}

			// mutual_information
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"mutual_information", 
					FunctionInfo(
						mutual_information, 
						forwardRange(parameterSet), 2)));
			}

			// mutual_information_pt
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// zData
					boost::any((Cell*)0),
					// timeWindowRadius
					boost::any((integer)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// zLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"mutual_information_pt", 
					FunctionInfo(
						mutual_information_pt, 
						forwardRange(parameterSet), 4)));
			}

			// mutual_information_p
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// zData
					boost::any((Cell*)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// zLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"mutual_information_p", 
					FunctionInfo(
						mutual_information_p, 
						forwardRange(parameterSet), 3)));
			}
				
			// transfer_entropy_t
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// wData
					boost::any((Cell*)0),
					// timeWindowRadius
					boost::any((integer)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// wLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"transfer_entropy_t", 
					FunctionInfo(
						transfer_entropy_t, 
						forwardRange(parameterSet), 4)));
			}

			// transfer_entropy_pt
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// zData
					boost::any((Cell*)0),
					// wData
					boost::any((Cell*)0),
					// timeWindowRadius
					boost::any((integer)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// zLag
					boost::any((integer)0),
					// wLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"transfer_entropy_pt", 
					FunctionInfo(
						transfer_entropy_pt, 
						forwardRange(parameterSet), 5)));
			}

			// transfer_entropy
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// wData
					boost::any((Cell*)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// wLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"transfer_entropy", 
					FunctionInfo(
						transfer_entropy, 
						forwardRange(parameterSet), 3)));
			}

			// transfer_entropy_p
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any((Cell*)0),
					// yData
					boost::any((Cell*)0),
					// zData
					boost::any((Cell*)0),
					// wData
					boost::any((Cell*)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// zLag
					boost::any((integer)0),
					// wLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"transfer_entropy_p", 
					FunctionInfo(
						transfer_entropy_p, 
						forwardRange(parameterSet), 4)));
			}

			initialized = true;
		}

	}

	boost::any functionCall(const std::string& name, const AnySet& argSet)
	{
		initializeFunctions();
	
		FunctionIterator iter = functionMap.find(name);
		if (iter == functionMap.end())
		{
			reportError("Undefined function '" + name + "'.");
			throw FunctionCall_Exception();
		}

		const FunctionInfo& info = iter->second;
		const integer callArgs = info.parameterSet.size();
		const integer inputArgs = argSet.size();

		if (inputArgs > callArgs)
		{
			reportError(name + "(): Too many arguments.");
			throw FunctionCall_Exception();
		}
		if (inputArgs < info.minArgs)
		{
			reportError(name + "(): Not enough arguments (min " + integerToString(info.minArgs) + ").");
			throw FunctionCall_Exception();
		}

		bool error = false;
		for (integer i = 0;i < inputArgs;++i)
		{
			if (argSet[i].type() != info.parameterSet[i].type())
			{
				reportError(name + "(): Parameter " + integerToString(i) + " is of the wrong type.");
				error = true;
			}
		}

		if (error)
		{
			throw FunctionCall_Exception();
		}

		AnySet callSet = argSet;
		for (integer i = inputArgs;i < callArgs;++i)
		{
			callSet.push_back(info.parameterSet[i]);
		}

		return info.callback(callSet, inputArgs);
	}

}
