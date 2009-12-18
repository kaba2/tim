#include "tim/console/functions.h"
#include "tim/console/console_parser.h"

#include <tim/core/differential_entropy.h>
#include <tim/core/mutual_information.h>
#include <tim/core/partial_mutual_information.h>
#include <tim/core/transfer_entropy.h>
#include <tim/core/partial_transfer_entropy.h>
#include <tim/core/divergence_wkv.h>

#include <fstream>

namespace Tim
{

	boost::any* load(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		std::string separatorSet = ",;";
		std::string fileName;

		bool error = false;
				
		try
		{
			fileName = boost::any_cast<std::string>(*argSet[0]);
			if (args > 1)
			{
				separatorSet = boost::any_cast<std::string>(*argSet[1]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}

		std::vector<real> data;
		real value = 0;

		std::ifstream file(fileName.c_str());
		if (!file.is_open())
		{
			std::cerr << "Error: Could not open data file " << fileName << "." << std::endl;
			return new boost::any;
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

		Signal* signal = new Signal(data.size(), 1);
		std::copy(data.begin(), data.end(),
			signal->data().begin());
			
		return new boost::any(signal);
	}

	boost::any* differential_entropy_kl(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* cell = 0;
		real maxRelativeError = 0;
		integer kNearest = 1;
		bool error = false;
				
		try
		{
			cell = boost::any_cast<Cell*>(*argSet[0]);
			if (args > 1)
			{
				maxRelativeError = boost::any_cast<real>(*argSet[1]);
			}
			if (args > 2)
			{
				kNearest = boost::any_cast<integer>(*argSet[2]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (maxRelativeError < 0)
		{
			reportError(location, "maxRelativeError must be non-negative.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}

		const real de = differentialEntropyKl(
			forwardRange(cell->begin(), cell->end()),
			maxRelativeError, kNearest);
			
		return new boost::any(de);
	}
	
	boost::any* differential_entropy_kl_t(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* cell = 0;
		real maxRelativeError = 0;
		integer kNearest = 1;
		integer timeWindowRadius = 0;
		
		bool error = false;
		try
		{
			cell = boost::any_cast<Cell*>(*argSet[0]);
			if (args > 1)
			{
				timeWindowRadius = boost::any_cast<integer>(*argSet[1]);
			}
			if (args > 2)
			{
				maxRelativeError = boost::any_cast<real>(*argSet[2]);
			}
			if (args > 3)
			{
				kNearest = boost::any_cast<integer>(*argSet[3]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (timeWindowRadius < 0)
		{
			reportError(location, "timeWindowRadius must be non-negative.");
			error = true;
		}
		
		if (maxRelativeError < 0)
		{
			reportError(location, "maxRelativeError must be non-negative.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		std::vector<real> deSet;
		
		temporalDifferentialEntropyKl(
			forwardRange(cell->begin(), cell->end()),
			timeWindowRadius,
			std::back_inserter(deSet),
			maxRelativeError, kNearest);
			
		Signal* signal = new Signal(deSet.size(), 1);
		std::copy(deSet.begin(), deSet.end(),
			signal->data().begin());
			
		return new boost::any(signal);
	}

	boost::any* differential_entropy_nk(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* cell = 0;
		real maxRelativeError = 0;
		bool error = false;
				
		try
		{
			cell = boost::any_cast<Cell*>(*argSet[0]);
			if (args > 1)
			{
				maxRelativeError = boost::any_cast<real>(*argSet[1]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (maxRelativeError < 0)
		{
			reportError(location, "maxRelativeError must be non-negative.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		integer dimension = 0;

		const real de = differentialEntropyNk(
			forwardRange(cell->begin(), cell->end()),
			maxRelativeError, 
			Default_NormBijection(),
			&dimension);
			
		Signal* signal = new Signal(2, 1);
		signal->data()(0) = de;
		signal->data()(1) = dimension;
			
		return new boost::any(signal);
	}

	boost::any* mutual_information_t(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		integer timeWindowRadius = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			if (args > 2)
			{
				timeWindowRadius = boost::any_cast<integer>(*argSet[2]);
			}
			if (args > 3)
			{
				xLag = boost::any_cast<integer>(*argSet[3]);
			}
			if (args > 4)
			{
				yLag = boost::any_cast<integer>(*argSet[4]);
			}
			if (args > 5)
			{
				kNearest = boost::any_cast<integer>(*argSet[5]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (timeWindowRadius < 0)
		{
			reportError(location, "timeWindowRadius must be non-negative.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		std::vector<real> miSet;
		
		temporalMutualInformation(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			timeWindowRadius,
			std::back_inserter(miSet),
			xLag, yLag,
			kNearest);
			
		Signal* signal = new Signal(miSet.size(), 1);
		std::copy(miSet.begin(), miSet.end(),
			signal->data().begin());
			
		return new boost::any(signal);
	}

	boost::any* mutual_information(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			if (args > 2)
			{
				xLag = boost::any_cast<integer>(*argSet[2]);
			}
			if (args > 3)
			{
				yLag = boost::any_cast<integer>(*argSet[3]);
			}
			if (args > 4)
			{
				kNearest = boost::any_cast<integer>(*argSet[4]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		const real mi = mutualInformation(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			xLag, yLag,
			kNearest);
			
		return new boost::any(mi);
	}

	boost::any* mutual_information_pt(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		Cell* zCell = 0;
		integer timeWindowRadius = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer zLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			zCell = boost::any_cast<Cell*>(*argSet[2]);
			if (args > 3)
			{
				timeWindowRadius = boost::any_cast<integer>(*argSet[3]);
			}
			if (args > 4)
			{
				xLag = boost::any_cast<integer>(*argSet[4]);
			}
			if (args > 5)
			{
				yLag = boost::any_cast<integer>(*argSet[5]);
			}
			if (args > 6)
			{
				zLag = boost::any_cast<integer>(*argSet[6]);
			}
			if (args > 7)
			{
				kNearest = boost::any_cast<integer>(*argSet[7]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (timeWindowRadius < 0)
		{
			reportError(location, "timeWindowRadius must be non-negative.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		std::vector<real> miSet;
		
		temporalPartialMutualInformation(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			forwardRange(zCell->begin(), zCell->end()),
			timeWindowRadius,
			std::back_inserter(miSet),
			xLag, yLag, zLag,
			kNearest);
			
		Signal* signal = new Signal(miSet.size(), 1);
		std::copy(miSet.begin(), miSet.end(),
			signal->data().begin());
			
		return new boost::any(signal);
	}

	boost::any* mutual_information_p(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		Cell* zCell = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer zLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			zCell = boost::any_cast<Cell*>(*argSet[2]);
			if (args > 3)
			{
				xLag = boost::any_cast<integer>(*argSet[3]);
			}
			if (args > 4)
			{
				yLag = boost::any_cast<integer>(*argSet[4]);
			}
			if (args > 5)
			{
				zLag = boost::any_cast<integer>(*argSet[5]);
			}
			if (args > 6)
			{
				kNearest = boost::any_cast<integer>(*argSet[6]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		const real pmi = partialMutualInformation(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			forwardRange(zCell->begin(), zCell->end()),
			xLag, yLag, zLag,
			kNearest);
			
		return new boost::any(pmi);
	}

	boost::any* transfer_entropy_t(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		Cell* wCell = 0;
		integer timeWindowRadius = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer wLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			wCell = boost::any_cast<Cell*>(*argSet[2]);
			if (args > 3)
			{
				timeWindowRadius = boost::any_cast<integer>(*argSet[3]);
			}
			if (args > 4)
			{
				xLag = boost::any_cast<integer>(*argSet[4]);
			}
			if (args > 5)
			{
				yLag = boost::any_cast<integer>(*argSet[5]);
			}
			if (args > 6)
			{
				wLag = boost::any_cast<integer>(*argSet[6]);
			}
			if (args > 7)
			{
				kNearest = boost::any_cast<integer>(*argSet[7]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (timeWindowRadius < 0)
		{
			reportError(location, "timeWindowRadius must be non-negative.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		std::vector<real> teSet;
		
		temporalTransferEntropy(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			forwardRange(wCell->begin(), wCell->end()),
			timeWindowRadius,
			std::back_inserter(teSet),
			xLag, yLag, wLag,
			kNearest);
			
		Signal* signal = new Signal(teSet.size(), 1);
		std::copy(teSet.begin(), teSet.end(),
			signal->data().begin());
			
		return new boost::any(signal);
	}

	boost::any* transfer_entropy_pt(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		Cell* zCell = 0;
		Cell* wCell = 0;
		integer timeWindowRadius = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer zLag = 0;
		integer wLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			zCell = boost::any_cast<Cell*>(*argSet[2]);
			wCell = boost::any_cast<Cell*>(*argSet[3]);
			if (args > 4)
			{
				timeWindowRadius = boost::any_cast<integer>(*argSet[4]);
			}
			if (args > 5)
			{
				xLag = boost::any_cast<integer>(*argSet[5]);
			}
			if (args > 6)
			{
				yLag = boost::any_cast<integer>(*argSet[6]);
			}
			if (args > 7)
			{
				zLag = boost::any_cast<integer>(*argSet[7]);
			}
			if (args > 8)
			{
				wLag = boost::any_cast<integer>(*argSet[8]);
			}
			if (args > 9)
			{
				kNearest = boost::any_cast<integer>(*argSet[9]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (timeWindowRadius < 0)
		{
			reportError(location, "timeWindowRadius must be non-negative.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
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
			
		Signal* signal = new Signal(teSet.size(), 1);
		std::copy(teSet.begin(), teSet.end(),
			signal->data().begin());
			
		return new boost::any(signal);
	}

	boost::any* transfer_entropy(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		Cell* wCell = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer wLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			wCell = boost::any_cast<Cell*>(*argSet[2]);
			if (args > 3)
			{
				xLag = boost::any_cast<integer>(*argSet[3]);
			}
			if (args > 4)
			{
				yLag = boost::any_cast<integer>(*argSet[4]);
			}
			if (args > 5)
			{
				wLag = boost::any_cast<integer>(*argSet[5]);
			}
			if (args > 6)
			{
				kNearest = boost::any_cast<integer>(*argSet[6]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		const real te = transferEntropy(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			forwardRange(wCell->begin(), wCell->end()),
			xLag, yLag, wLag,
			kNearest);
			
		return new boost::any(te);
	}

	boost::any* transfer_entropy_p(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		Cell* zCell = 0;
		Cell* wCell = 0;
		integer xLag = 0;
		integer yLag = 0;
		integer zLag = 0;
		integer wLag = 0;
		integer kNearest = 1;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
			zCell = boost::any_cast<Cell*>(*argSet[2]);
			wCell = boost::any_cast<Cell*>(*argSet[3]);
			if (args > 4)
			{
				xLag = boost::any_cast<integer>(*argSet[4]);
			}
			if (args > 5)
			{
				yLag = boost::any_cast<integer>(*argSet[5]);
			}
			if (args > 6)
			{
				zLag = boost::any_cast<integer>(*argSet[6]);
			}
			if (args > 7)
			{
				wLag = boost::any_cast<integer>(*argSet[7]);
			}
			if (args > 8)
			{
				kNearest = boost::any_cast<integer>(*argSet[8]);
			}
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (kNearest < 1)
		{
			reportError(location, "kNearest must be at least 1.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		const real pte = partialTransferEntropy(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()),
			forwardRange(zCell->begin(), zCell->end()),
			forwardRange(wCell->begin(), wCell->end()),
			xLag, yLag, zLag, wLag,
			kNearest);
			
		return new boost::any(pte);
	}

	boost::any* divergence_wkv(const YYLTYPE& location, const AnySet& argSet)
	{
		const integer args = argSet.size();
		
		Cell* xCell = 0;
		Cell* yCell = 0;
		
		bool error = false;
		try
		{
			xCell = boost::any_cast<Cell*>(*argSet[0]);
			yCell = boost::any_cast<Cell*>(*argSet[1]);
		}
		catch(const boost::bad_any_cast&)
		{
			reportError(location, "Invalid argument type.");
			error = true;
		}
		
		if (error)
		{
			return new boost::any;
		}
		
		const real div = divergenceWkv(
			forwardRange(xCell->begin(), xCell->end()),
			forwardRange(yCell->begin(), yCell->end()));
			
		return new boost::any(div);
	}

}
