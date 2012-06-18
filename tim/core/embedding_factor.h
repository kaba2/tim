// Description: Embedding dimension via false neighbors

// TIM 1.2.0
// Kalle Rutanen
// http://kaba.hilvi.org
// Copyright (c) 2009 - 2011
//
// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published 
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef TIM_EMBEDDING_FACTOR_FN_H
#define TIM_EMBEDDING_FACTOR_FN_H

#include "pastel/sys/mytypes.h"
#include "pastel/sys/forwarditerator_range.h"

namespace Tim
{

	template <typename SignalPtr_Iterator>
	integer embeddingFactorFn(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet);

}

#include "tim/core/embedding_factor_fn.hpp"

#endif
