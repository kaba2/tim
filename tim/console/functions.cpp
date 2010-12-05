#include "tim/console/functions.h"
#include "tim/console/console_parser.h"
#include "tim/console/errorlog.h"

#include <tim/core/differential_entropy.h>
#include <tim/core/renyi_entropy.h>
#include <tim/core/tsallis_entropy.h>
#include <tim/core/mutual_information.h>
#include <tim/core/partial_mutual_information.h>
#include <tim/core/transfer_entropy.h>
#include <tim/core/partial_transfer_entropy.h>
#include <tim/core/divergence_wkv.h>
#include <tim/core/embed.h>

#include <pastel/sys/string_tools.h>

#include <fstream>
#include <map>
#include <string>

#include <boost/function.hpp>

namespace Tim
{

	namespace
	{

		boost::any delay_embed(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				signalIndex,
				kIndex,
				dtIndex
			};

			const SignalPtr signal = 
				boost::any_cast<SignalPtr>(argSet[signalIndex]);
			const integer k = 
				boost::any_cast<integer>(argSet[kIndex]);
			const integer dt = 
				boost::any_cast<integer>(argSet[dtIndex]);

			// Check parameters.

			bool error = false;
			if (k <= 0)
			{
				reportError("'k' must be positive.");
				error = true;
			}
			if (dt <= 0)
			{
				reportError("'dt' must be positive.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			return boost::any(delayEmbed(signal, k, dt));
		}

		boost::any delay_embed_future(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				signalIndex,
				kIndex,
				dtIndex
			};

			const SignalPtr signal = 
				boost::any_cast<SignalPtr>(argSet[signalIndex]);
			const integer k = 
				boost::any_cast<integer>(argSet[kIndex]);
			const integer dt = 
				boost::any_cast<integer>(argSet[dtIndex]);

			// Check parameters.

			bool error = false;
			if (k <= 0)
			{
				reportError("'k' must be positive.");
				error = true;
			}
			if (dt <= 0)
			{
				reportError("'dt' must be positive.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			return boost::any(delayEmbedFuture(signal, k, dt));
		}

		boost::any read_csv(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				fileNameIndex,
				samplesIndex,
				dimensionIndex,
				trialsIndex,
				seriesIndex,
				separatorSetIndex,
				orderIndex
			};

			const std::string fileName = 
				boost::any_cast<std::string>(argSet[fileNameIndex]);
			integer samples = 
				boost::any_cast<integer>(argSet[samplesIndex]);
			const integer dimension = 
				boost::any_cast<integer>(argSet[dimensionIndex]);
			const integer trials = 
				boost::any_cast<integer>(argSet[trialsIndex]);
			const integer series = 
				boost::any_cast<integer>(argSet[seriesIndex]);
			const std::string separatorSet = 
				boost::any_cast<std::string>(argSet[separatorSetIndex]);
			const SignalPtr order = 
				boost::any_cast<SignalPtr>(argSet[orderIndex]);

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
			if (order->data().size() != 4)
			{
				reportError("'order' must have exactly 4 elements.");
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
				error = true;
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

			const integer samplesToRead = product(width);
			if (!samplesUnknown)
			{
				data.reserve(samplesToRead);
			}

			integer readSamples = 0;
			bool finished = false;
			while(!finished)
			{
				// Read one word.
				std::string text;
				file >> text;
				if (!file)
				{
					// There is no more data in the file.
					finished = true;
					break;
				}
				// We found a word.
				while(!text.empty())
				{
					// From the word, try to read a real
					// number.
					integer indexEnd = 0;
					const real value = stringToReal(text, &indexEnd);
					if (indexEnd > 0)
					{
						// Found some real value.
						data.push_back(value);
						++readSamples;
						if (!samplesUnknown && readSamples == samplesToRead)
						{
							// No need to read more samples since
							// we already got all we need to form the
							// signals.
							finished = true;
							break;
						}
						// Concentrate on the rest of the word.
						text = text.substr(indexEnd);
					}
					else
					{
						// See if we can read a separator instead.
						char separator = text[0];
						if (separatorSet.find(separator) == std::string::npos)
						{
							// There was something, but it was not a separator either.
							reportError(std::string("Invalid separator ") + separator + ".");
							throw FunctionCall_Exception();
						}
						text = text.substr(1);
					}
				}
			}

			if (!samplesUnknown && readSamples != samplesToRead)
			{
				reportError(std::string("Not enough samples could be read.") + 
					" Needed " + integerToString(samplesToRead) + 
					", read " + integerToString(readSamples) + ".");
				throw FunctionCall_Exception();
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

			CellPtr cellArray(new Cell(trials, series));

			for (integer y = 0;y < series;++y)
			{
				for (integer x = 0;x < trials;++x)
				{
					integer nans = 0;
					for (integer j = 0;j < samples;++j)
					{
						const integer offset = dot(stride, Vector<integer, 4>(j, 0, x, y));
						if (!StdExt::isNan(data[offset]))
						{
							break;
						}
						++nans;
					}
					SignalPtr signal = SignalPtr(new Signal(samples - nans, dimension, nans));
					for (integer i = 0;i < dimension;++i)
					{
						for (integer j = nans;j < samples;++j)
						{
							const integer offset = dot(stride, Vector<integer, 4>(j, i, x, y));
							signal->data()(j - nans, i) = data[offset];
						}
					}
					(*cellArray)(x, y) = signal;
				}
			}

			
			return boost::any(cellArray);
		}

		boost::any write_csv(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			const std::string fileName = 
				boost::any_cast<std::string>(argSet[0]);
			CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[1]);
			const SignalPtr order = 
				boost::any_cast<SignalPtr>(argSet[2]);
			std::string separator0 = 
				boost::any_cast<std::string>(argSet[3]);
			const std::string separator1 = 
				boost::any_cast<std::string>(argSet[4]);
			const std::string separator2 = 
				boost::any_cast<std::string>(argSet[5]);
			const std::string separator3 = 
				boost::any_cast<std::string>(argSet[6]);
			const std::string prefix0 = 
				boost::any_cast<std::string>(argSet[7]);
			const std::string suffix0 = 
				boost::any_cast<std::string>(argSet[8]);
			const std::string prefix1 = 
				boost::any_cast<std::string>(argSet[9]);
			const std::string suffix1 = 
				boost::any_cast<std::string>(argSet[10]);
			const std::string prefix2 = 
				boost::any_cast<std::string>(argSet[11]);
			const std::string suffix2 = 
				boost::any_cast<std::string>(argSet[12]);
			const std::string prefix3 = 
				boost::any_cast<std::string>(argSet[13]);
			const std::string suffix3 = 
				boost::any_cast<std::string>(argSet[14]);
			const std::string prefix4 = 
				boost::any_cast<std::string>(argSet[15]);
			const std::string suffix4 = 
				boost::any_cast<std::string>(argSet[16]);

			// Check parameters.

			bool error = false;
			if (order->data().size() != 4)
			{
				reportError("'order' must have exactly 4 elements.");
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
				error = true;
			}

			std::ofstream file(fileName.c_str());
			if (!file.is_open())
			{
				reportError("Could not open data file " + fileName + " for writing.");
				error = true;
			}

			if (!equalDimension(range(cell->begin(), cell->end())))
			{
				reportError("The signals in the cell-array must have equal dimension.");
				error = true;
			}

			if (!equalSamples(range(cell->begin(), cell->end())))
			{
				reportError("The signals in the cell-array must have equal number of samples.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			const integer trials = cell->width();
			const integer series = cell->height();
			const integer nans = (*cell)(0)->t();
			const integer samples = (*cell)(0)->samples() + nans;
			const integer dimension = (*cell)(0)->dimension();

			Vector<integer, 4> width = permute(
				Vector<integer, 4>(samples, dimension, trials, series),
				permutation);

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

			const integer dataSize = trials * series * samples * dimension;

			std::vector<real> data(dataSize);

			// Write the values into a continuous array.

			for (integer y = 0;y < series;++y)
			{
				for (integer x = 0;x < trials;++x)
				{
					const SignalPtr signal = (*cell)(x, y);
					for (integer i = 0;i < dimension;++i)
					{
						for (integer j = 0;j < nans;++j)
						{
							const integer offset = dot(stride, Vector<integer, 4>(j, i, x, y));
							data[offset] = nan<real>();
						}
						for (integer j = nans;j < samples;++j)
						{
							const integer offset = dot(stride, Vector<integer, 4>(j, i, x, y));
							data[offset] = signal->data()(j - nans, i);
						}
					}
				}
			}

			// Write the array into the file.

			if (separator0.empty())
			{
				separator0 = " ";
			}

			integer index = 0;

			file << prefix4;
			for (integer i3 = 0;i3 < width[3];++i3)
			{
				if (i3 > 0)
				{
					file << separator3;
					file << std::endl;
				}
				file << prefix3;
				for (integer i2 = 0;i2 < width[2];++i2)
				{
					if (i2 > 0)
					{
						file << separator2;
						file << std::endl;
					}
					file << prefix2;
					for (integer i1 = 0;i1 < width[1];++i1)
					{
						if (i1 > 0)
						{
							file << separator1;
							file << std::endl;
						}
						file << prefix1;
						for (integer i0 = 0;i0 < width[0];++i0)
						{
							if (i0 > 0)
							{
								file << separator0;
							}
							file << prefix0;
							file << realToString(data[index]);
							file << suffix0;
							++index;
						}
						file << suffix1;
					}
					file << suffix2;
				}
				file << suffix3;
			}
			file << suffix4;

			return boost::any((integer)0);
		}

		boost::any differential_entropy_kl(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				cellIndex,
				kNearestIndex
			};

			const CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[cellIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);

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

			const real de = differentialEntropyKl(
				range(cell->begin(), cell->end()),
				kNearest);
				
			return boost::any(de);
		}
		
		boost::any differential_entropy_kl_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				cellIndex,
				timeWindowRadiusIndex,
				kNearestIndex,
				filterIndex
			};

			const CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[cellIndex]);
			const integer timeWindowRadius = 
				boost::any_cast<integer>(argSet[timeWindowRadiusIndex]);;
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			const SignalPtr filter =
				boost::any_cast<SignalPtr>(argSet[filterIndex]);
			
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

			if (!odd(filter->data().size()))
			{
				reportError("The number of elements in filter must be odd.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const SignalPtr de = temporalDifferentialEntropyKl(
				range(cell->begin(), cell->end()),
				timeWindowRadius,
				kNearest,
				Default_NormBijection(),
				range(filter->data().begin(), filter->data().end()));
				
			return boost::any(de);
		}

		boost::any differential_entropy_nk(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				cellIndex
			};

			const CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[cellIndex]);

			// Check parameters.

			bool error = false;
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			integer dimension = 0;

			const real de = differentialEntropyNk(
				range(cell->begin(), cell->end()),
				Default_NormBijection(),
				&dimension);
				
			SignalPtr signal = SignalPtr(new Signal(2, 1));
			signal->data()(0) = de;
			signal->data()(1) = dimension;
				
			return boost::any(signal);
		}

		boost::any renyi_entropy_lps(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				cellIndex,
				qIndex,
				kNearestIndex
			};

			const CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[cellIndex]);
			const real q = 
				boost::any_cast<real>(argSet[qIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);

			// Check parameters.

			bool error = false;
			if (q <= 0)
			{
				reportError("q must be positive.");
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

			const real de = renyiEntropyLps(
				range(cell->begin(), cell->end()),
				q,
				kNearest);
				
			return boost::any(de);
		}
		
		boost::any renyi_entropy_lps_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				cellIndex,
				timeWindowRadiusIndex,
				qIndex,
				kNearestIndex,
				filterIndex
			};

			const CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[cellIndex]);
			const integer timeWindowRadius = 
				boost::any_cast<integer>(argSet[timeWindowRadiusIndex]);
			const real q = 
				boost::any_cast<real>(argSet[qIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			const SignalPtr filter =
				boost::any_cast<SignalPtr>(argSet[filterIndex]);
			
			// Check parameters.

			bool error = false;
			if (timeWindowRadius < 0)
			{
				reportError("timeWindowRadius must be non-negative.");
				error = true;
			}
			
			if (q <= 0)
			{
				reportError("q must be positive.");
				error = true;
			}

			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}
			
			if (!odd(filter->data().size()))
			{
				reportError("The number of elements in filter must be odd.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const SignalPtr re = temporalRenyiEntropyLps(
				range(cell->begin(), cell->end()),
				timeWindowRadius,
				q,
				kNearest,
				range(filter->data().begin(), filter->data().end()));
				
			return boost::any(re);
		}

		boost::any tsallis_entropy_lps(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.
			
			enum
			{
				cellIndex,
				qIndex,
				kNearestIndex
			};

			const CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[cellIndex]);
			const real q = 
				boost::any_cast<real>(argSet[qIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);

			// Check parameters.

			bool error = false;
			if (q <= 0)
			{
				reportError("q must be positive.");
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

			const real de = tsallisEntropyLps(
				range(cell->begin(), cell->end()),
				q,
				kNearest);
				
			return boost::any(de);
		}
		
		boost::any tsallis_entropy_lps_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				cellIndex,
				timeWindowRadiusIndex,
				qIndex,
				kNearestIndex,
				filterIndex
			};

			const CellPtr cell = 
				boost::any_cast<CellPtr>(argSet[cellIndex]);
			const integer timeWindowRadius = 
				boost::any_cast<integer>(argSet[timeWindowRadiusIndex]);
			const real q = 
				boost::any_cast<real>(argSet[qIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			const SignalPtr filter =
				boost::any_cast<SignalPtr>(argSet[filterIndex]);
			
			// Check parameters.

			bool error = false;
			if (timeWindowRadius < 0)
			{
				reportError("timeWindowRadius must be non-negative.");
				error = true;
			}
			
			if (q <= 0)
			{
				reportError("q must be positive.");
				error = true;
			}

			if (kNearest < 1)
			{
				reportError("kNearest must be at least 1.");
				error = true;
			}

			if (!odd(filter->data().size()))
			{
				reportError("The number of elements in filter must be odd.");
				error = true;
			}
			
			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const SignalPtr estimate = temporalTsallisEntropyLps(
				range(cell->begin(), cell->end()),
				timeWindowRadius,
				q,
				kNearest,
				range(filter->data().begin(), filter->data().end()));
				
			return boost::any(estimate);
		}

		boost::any mutual_information_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				timeWindowRadiusIndex,
				xLagIndex,
				yLagIndex,
				kNearestIndex,
				filterIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const integer timeWindowRadius = 
				boost::any_cast<integer>(argSet[timeWindowRadiusIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			const SignalPtr filter =
				boost::any_cast<SignalPtr>(argSet[filterIndex]);

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
			
			if (!odd(filter->data().size()))
			{
				reportError("The number of elements in filter must be odd.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const SignalPtr mi = temporalMutualInformation(
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				timeWindowRadius,
				xLag, yLag,
				kNearest,
				range(filter->data().begin(), filter->data().end()));
				
			return boost::any(mi);
		}

		boost::any mutual_information(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				xLagIndex,
				yLagIndex,
				kNearestIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);

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
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				xLag, yLag,
				kNearest);
				
			return boost::any(mi);
		}

		boost::any mutual_information_pt(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				zCellIndex,
				timeWindowRadiusIndex,
				xLagIndex,
				yLagIndex,
				zLagIndex,
				kNearestIndex,
				filterIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const CellPtr zCell = 
				boost::any_cast<CellPtr>(argSet[zCellIndex]);
			const integer timeWindowRadius = 
				boost::any_cast<integer>(argSet[timeWindowRadiusIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer zLag = 
				boost::any_cast<integer>(argSet[zLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			const SignalPtr filter =
				boost::any_cast<SignalPtr>(argSet[filterIndex]);
			
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
			
			if (!odd(filter->data().size()))
			{
				reportError("The number of elements in filter must be odd.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const SignalPtr mi = temporalPartialMutualInformation(
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				range(zCell->begin(), zCell->end()),
				timeWindowRadius,
				xLag, yLag, zLag,
				kNearest,
				range(filter->data().begin(), filter->data().end()));
				
			return boost::any(mi);
		}

		boost::any mutual_information_p(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				zCellIndex,
				xLagIndex,
				yLagIndex,
				zLagIndex,
				kNearestIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const CellPtr zCell = 
				boost::any_cast<CellPtr>(argSet[zCellIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer zLag = 
				boost::any_cast<integer>(argSet[zLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			
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
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				range(zCell->begin(), zCell->end()),
				xLag, yLag, zLag,
				kNearest);
				
			return boost::any(pmi);
		}

		boost::any transfer_entropy_t(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				wCellIndex,
				timeWindowRadiusIndex,
				xLagIndex,
				yLagIndex,
				wLagIndex,
				kNearestIndex,
				filterIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const CellPtr wCell = 
				boost::any_cast<CellPtr>(argSet[wCellIndex]);
			const integer timeWindowRadius = 
				boost::any_cast<integer>(argSet[timeWindowRadiusIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer wLag = 
				boost::any_cast<integer>(argSet[wLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			const SignalPtr filter =
				boost::any_cast<SignalPtr>(argSet[filterIndex]);
			
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
			
			if (!odd(filter->data().size()))
			{
				reportError("The number of elements in filter must be odd.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const SignalPtr te = temporalTransferEntropy(
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				range(wCell->begin(), wCell->end()),
				timeWindowRadius,
				xLag, yLag, wLag,
				kNearest,
				range(filter->data().begin(), filter->data().end()));
				
			return boost::any(te);
		}

		boost::any transfer_entropy_pt(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				zCellIndex,
				wCellIndex,
				timeWindowRadiusIndex,
				xLagIndex,
				yLagIndex,
				zLagIndex,
				wLagIndex,
				kNearestIndex,
				filterIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const CellPtr zCell = 
				boost::any_cast<CellPtr>(argSet[zCellIndex]);
			const CellPtr wCell = 
				boost::any_cast<CellPtr>(argSet[wCellIndex]);
			const integer timeWindowRadius = 
				boost::any_cast<integer>(argSet[timeWindowRadiusIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer zLag = 
				boost::any_cast<integer>(argSet[zLagIndex]);
			const integer wLag = 
				boost::any_cast<integer>(argSet[wLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);
			const SignalPtr filter =
				boost::any_cast<SignalPtr>(argSet[filterIndex]);

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
			
			if (!odd(filter->data().size()))
			{
				reportError("The number of elements in filter must be odd.");
				error = true;
			}

			if (error)
			{
				throw FunctionCall_Exception();
			}

			// Compute.
			
			const SignalPtr pte = temporalPartialTransferEntropy(
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				range(zCell->begin(), zCell->end()),
				range(wCell->begin(), wCell->end()),
				timeWindowRadius,
				xLag, yLag, zLag, wLag,
				kNearest,
				range(filter->data().begin(), filter->data().end()));
				
			return boost::any(pte);
		}

		boost::any transfer_entropy(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				wCellIndex,
				xLagIndex,
				yLagIndex,
				wLagIndex,
				kNearestIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const CellPtr wCell = 
				boost::any_cast<CellPtr>(argSet[wCellIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer wLag = 
				boost::any_cast<integer>(argSet[wLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);

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
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				range(wCell->begin(), wCell->end()),
				xLag, yLag, wLag,
				kNearest);
				
			return boost::any(te);
		}

		boost::any transfer_entropy_p(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex,
				zCellIndex,
				wCellIndex,
				xLagIndex,
				yLagIndex,
				zLagIndex,
				wLagIndex,
				kNearestIndex
			};
			
			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);
			const CellPtr zCell = 
				boost::any_cast<CellPtr>(argSet[zCellIndex]);
			const CellPtr wCell = 
				boost::any_cast<CellPtr>(argSet[wCellIndex]);
			const integer xLag = 
				boost::any_cast<integer>(argSet[xLagIndex]);
			const integer yLag = 
				boost::any_cast<integer>(argSet[yLagIndex]);
			const integer zLag = 
				boost::any_cast<integer>(argSet[zLagIndex]);
			const integer wLag = 
				boost::any_cast<integer>(argSet[wLagIndex]);
			const integer kNearest = 
				boost::any_cast<integer>(argSet[kNearestIndex]);

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
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()),
				range(zCell->begin(), zCell->end()),
				range(wCell->begin(), wCell->end()),
				xLag, yLag, zLag, wLag,
				kNearest);
				
			return boost::any(pte);
		}

		boost::any divergence_wkv(const AnySet& argSet, integer passedArgs)
		{
			// Retrieve parameters.

			enum
			{
				xCellIndex,
				yCellIndex
			};

			const CellPtr xCell = 
				boost::any_cast<CellPtr>(argSet[xCellIndex]);
			const CellPtr yCell = 
				boost::any_cast<CellPtr>(argSet[yCellIndex]);

			// Compute.

			const real div = divergenceWkv(
				range(xCell->begin(), xCell->end()),
				range(yCell->begin(), yCell->end()));
				
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

			SignalPtr filter(new Signal(1, 1));
			filter->data()(0) = 1;

			// delay_embed
			{
				boost::any parameterSet[] =
				{
					// signal
					boost::any(SignalPtr()),
					// k
					boost::any(integer(1)),
					// dt
					boost::any(integer(1))
				};

				functionMap.insert(
					std::make_pair(
					"delay_embed", 
					FunctionInfo(
						delay_embed, 
						range(parameterSet), 2)));
			}

			// delay_embed_future
			{
				boost::any parameterSet[] =
				{
					// signal
					boost::any(SignalPtr()),
					// k
					boost::any(integer(1)),
					// dt
					boost::any(integer(1))
				};

				functionMap.insert(
					std::make_pair(
					"delay_embed_future", 
					FunctionInfo(
						delay_embed_future, 
						range(parameterSet), 2)));
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
						range(parameterSet), 1)));
			}

			// write_csv
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
					// data
					boost::any(CellPtr()),
					// order
					boost::any(order),
					// separatorLevel0
					boost::any(std::string(",")),
					// separatorLevel1
					boost::any(std::string(",")),
					// separatorLevel2
					boost::any(std::string(",")),
					// separatorLevel3
					boost::any(std::string(",")),
					// prefixLevel0
					boost::any(std::string("")),
					// suffixLevel0
					boost::any(std::string("")),
					// prefixLevel1
					boost::any(std::string("")),
					// suffixLevel1
					boost::any(std::string("")),
					// prefixLevel2
					boost::any(std::string("")),
					// suffixLevel2
					boost::any(std::string("")),
					// prefixLevel3
					boost::any(std::string("")),
					// suffixLevel3
					boost::any(std::string("")),
					// prefixLevel4
					boost::any(std::string("")),
					// suffixLevel4
					boost::any(std::string(""))
				};
			
				functionMap.insert(
					std::make_pair(
					"write_csv", 
					FunctionInfo(
						write_csv, 
						range(parameterSet), 2)));
			}

			// differential_entropy_kl
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any(CellPtr()),
					// kNearest
					boost::any((integer)1)
				};
			
				functionMap.insert(
					std::make_pair(
					"differential_entropy_kl", 
					FunctionInfo(
						differential_entropy_kl, 
						range(parameterSet), 1)));
			}
			
			// differential_entropy_kl_t
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any(CellPtr()),
					// timeWindowRadius
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1),
					// filter
					boost::any(filter)
				};
			
				functionMap.insert(
					std::make_pair(
					"differential_entropy_kl_t", 
					FunctionInfo(
						differential_entropy_kl_t, 
						range(parameterSet), 2)));
			}
				
			// differential_entropy_nk
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any(CellPtr()),
				};
			
				functionMap.insert(
					std::make_pair(
					"differential_entropy_nk", 
					FunctionInfo(
						differential_entropy_nk, 
						range(parameterSet), 1)));
			}

			// renyi_entropy_lps
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any(CellPtr()),
					// q
					boost::any((real)2),
					// kNearestSuggestion
					boost::any((integer)0)
				};
			
				functionMap.insert(
					std::make_pair(
					"renyi_entropy_lps", 
					FunctionInfo(
						renyi_entropy_lps, 
						range(parameterSet), 1)));
			}
			
			// renyi_entropy_lps_t
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any(CellPtr()),
					// timeWindowRadius
					boost::any((integer)0),
					// q
					boost::any((real)2),
					// kNearestSuggestion
					boost::any((integer)0),
					// filter
					boost::any(filter)
				};
			
				functionMap.insert(
					std::make_pair(
					"renyi_entropy_lps_t", 
					FunctionInfo(
						renyi_entropy_lps_t, 
						range(parameterSet), 2)));
			}

			// tsallis_entropy_lps
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any(CellPtr()),
					// q
					boost::any((real)2),
					// kNearestSuggestion
					boost::any((integer)0)
				};
			
				functionMap.insert(
					std::make_pair(
					"tsallis_entropy_lps", 
					FunctionInfo(
						tsallis_entropy_lps, 
						range(parameterSet), 1)));
			}
			
			// tsallis_entropy_lps_t
			{
				boost::any parameterSet[] =
				{
					// data
					boost::any(CellPtr()),
					// timeWindowRadius
					boost::any((integer)0),
					// q
					boost::any((real)2),
					// kNearestSuggestion
					boost::any((integer)0),
					// filter
					boost::any(filter)
				};
			
				functionMap.insert(
					std::make_pair(
					"tsallis_entropy_lps_t", 
					FunctionInfo(
						tsallis_entropy_lps_t, 
						range(parameterSet), 2)));
			}

			// divergence_wkv
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr())
				};
			
				functionMap.insert(
					std::make_pair(
					"divergence_wkv", 
					FunctionInfo(
						divergence_wkv, 
						range(parameterSet), 2)));
			}

			// mutual_information_t
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
					// timeWindowRadius
					boost::any((integer)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1),
					// filter
					boost::any(filter)
				};
			
				functionMap.insert(
					std::make_pair(
					"mutual_information_t", 
					FunctionInfo(
						mutual_information_t, 
						range(parameterSet), 3)));
			}

			// mutual_information
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
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
						range(parameterSet), 2)));
			}

			// mutual_information_pt
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
					// zData
					boost::any(CellPtr()),
					// timeWindowRadius
					boost::any((integer)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// zLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1),
					// filter
					boost::any(filter)
				};
			
				functionMap.insert(
					std::make_pair(
					"mutual_information_pt", 
					FunctionInfo(
						mutual_information_pt, 
						range(parameterSet), 4)));
			}

			// mutual_information_p
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
					// zData
					boost::any(CellPtr()),
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
						range(parameterSet), 3)));
			}
				
			// transfer_entropy_t
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
					// wData
					boost::any(CellPtr()),
					// timeWindowRadius
					boost::any((integer)0),
					// xLag
					boost::any((integer)0),
					// yLag
					boost::any((integer)0),
					// wLag
					boost::any((integer)0),
					// kNearest
					boost::any((integer)1),
					// filter
					boost::any(filter)
				};
			
				functionMap.insert(
					std::make_pair(
					"transfer_entropy_t", 
					FunctionInfo(
						transfer_entropy_t, 
						range(parameterSet), 4)));
			}

			// transfer_entropy_pt
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
					// zData
					boost::any(CellPtr()),
					// wData
					boost::any(CellPtr()),
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
					boost::any((integer)1),
					// filter
					boost::any(filter)
				};
			
				functionMap.insert(
					std::make_pair(
					"transfer_entropy_pt", 
					FunctionInfo(
						transfer_entropy_pt, 
						range(parameterSet), 5)));
			}

			// transfer_entropy
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
					// wData
					boost::any(CellPtr()),
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
						range(parameterSet), 3)));
			}

			// transfer_entropy_p
			{
				boost::any parameterSet[] =
				{
					// xData
					boost::any(CellPtr()),
					// yData
					boost::any(CellPtr()),
					// zData
					boost::any(CellPtr()),
					// wData
					boost::any(CellPtr()),
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
						range(parameterSet), 4)));
			}

			initialized = true;
		}

		std::string typeToString(const std::type_info& that)
		{
			if (that == typeid(std::string))
			{
				return "string";
			}
			if (that == typeid(integer))
			{
				return "integer";
			}
			if (that == typeid(real))
			{
				return "real";
			}
			if (that == typeid(SignalPtr))
			{
				return "signal";
			}
			if (that == typeid(CellPtr))
			{
				return "cell-array";
			}
			
			return "unknown";
		}

	}

	boost::any functionCall(const std::string& name, const AnySet& argSet)
	{
		initializeFunctions();

		ErrorLog_Namespace errorName(name + "(): ");
	
		FunctionIterator iter = functionMap.find(name);
		if (iter == functionMap.end())
		{
			reportError("Undefined function.");
			throw FunctionCall_Exception();
		}

		const FunctionInfo& info = iter->second;
		const integer callArgs = info.parameterSet.size();
		const integer inputArgs = argSet.size();

		if (inputArgs > callArgs)
		{
			reportError("Too many arguments.");
			throw FunctionCall_Exception();
		}
		if (inputArgs < info.minArgs)
		{
			reportError("Not enough arguments (min " + integerToString(info.minArgs) + ").");
			throw FunctionCall_Exception();
		}

		bool error = false;
		for (integer i = 0;i < inputArgs;++i)
		{
			if (argSet[i].type() != info.parameterSet[i].type())
			{
				reportError(
					"Argument " + integerToString(i) + " is of the wrong type. " + 
					"Expected " + typeToString(info.parameterSet[i].type()) + 
					", got " + typeToString(argSet[i].type()) + ".");
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

		boost::any result = info.callback(callSet, inputArgs);

		return result;
	}

}
