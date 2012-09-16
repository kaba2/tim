// Description: Linear reconstruction of NaNs.

#ifndef TIM_RECONSTRUCTION_H
#define TIM_RECONSTRUCTION_H

#include "pastel/sys/range.h"

namespace Tim
{

	template <typename Real_InputIterator>
	void reconstruct(
		const boost::iterator_range<Real_InputIterator>& data);

}

#include "tim/core/reconstruction.hpp"

#endif
