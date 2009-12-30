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

		struct FunctionInfo
		{
			typedef boost::function<boost::any(const AnySet& argSet)> Callback;

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

			// load
			{
				boost::any parameterSet[] =
				{
					// name
					boost::any(std::string()),
					// separatorSet
					boost::any(std::string(",;")),
				};
			
				functionMap.insert(
					std::make_pair(
					"load", 
					FunctionInfo(
						load, 
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

		return info.callback(callSet);
	}

	boost::any load(const AnySet& argSet)
	{
		std::string fileName = boost::any_cast<std::string>(argSet[0]);
		std::string separatorSet = boost::any_cast<std::string>(argSet[1]);

		std::vector<real> data;
		real value = 0;

		std::ifstream file(fileName.c_str());
		if (!file.is_open())
		{
			std::cerr << "Error: Could not open data file " << fileName << "." << std::endl;
			return boost::any();
		}

		while(file >> value)
		{
			data.push_back(value);
			char separator = 0;
			if (!(file >> separator))
			{
				break;
			}
			if (separatorSet.find(separator) == std::string::npos)
			{
				std::cerr << "Error: Data file " << fileName << " had an invalid separator " 
					<< separator << "." << std::endl;
				break;
			}
		}

		SignalPtr signal = SignalPtr(new Signal(data.size(), 1));
		std::copy(data.begin(), data.end(),
			signal->data().begin());
			
		return boost::any(signal);
	}

	boost::any differential_entropy_kl(const AnySet& argSet)
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
			return boost::any();
		}

		// Compute.

		const real de = differentialEntropyKl(
			forwardRange(cell->begin(), cell->end()),
			maxRelativeError, kNearest);
			
		return boost::any(de);
	}
	
	boost::any differential_entropy_kl_t(const AnySet& argSet)
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
			return boost::any();
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

	boost::any differential_entropy_nk(const AnySet& argSet)
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
			return boost::any();
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

	boost::any mutual_information_t(const AnySet& argSet)
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
			return boost::any();
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

	boost::any mutual_information(const AnySet& argSet)
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
			return boost::any();
		}

		// Compute.

		const real mi = mutualInformation(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			xLag, yLag,
			kNearest);
			
		return boost::any(mi);
	}

	boost::any mutual_information_pt(const AnySet& argSet)
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
			return boost::any();
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

	boost::any mutual_information_p(const AnySet& argSet)
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
			return boost::any();
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

	boost::any transfer_entropy_t(const AnySet& argSet)
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
			return boost::any();
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

	boost::any transfer_entropy_pt(const AnySet& argSet)
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
			return boost::any();
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

	boost::any transfer_entropy(const AnySet& argSet)
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
			return boost::any();
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

	boost::any transfer_entropy_p(const AnySet& argSet)
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
			return boost::any();
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

	boost::any divergence_wkv(const AnySet& argSet)
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

}
