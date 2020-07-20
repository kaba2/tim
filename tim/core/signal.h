// Description: Time series

#ifndef TIM_SIGNAL_H
#define TIM_SIGNAL_H

#include <pastel/math/matrix/matrix.h>

#include <pastel/sys/range.h>

#include "tim/core/mytypes.h"

#include <vector>

namespace Tim
{

	//! Time series
	class TIM Signal
	{
	public:
		//! Construct an empty signal.
		Signal() = default;
		Signal(const Signal& that) = default;
		Signal(Signal&& that) = default;

		//! Construct a zero signal.
		/*!
		Preconditions:
		samples >= 0
		dimension >= 0
		*/
		Signal(MatrixView<dreal> view, integer t = 0)
			: data_(view)
			, t_(t)
		{
		}

		//! Swaps two signals.
		void swap(Signal& that)
		{
			data_.swap(that.data_);
			std::swap(t_, that.t_);
		}

		//! Assigns from another signal.
		Signal& operator=(Signal that)
		{
			swap(that);
			return *this;
		}

		//! Returns the dimension of the signal.
		integer dimension() const
		{
			return data_.rows();
		}

		//! Returns the number of samples in the signal.
		integer samples() const
		{
			return data_.cols();
		}

		//! Returns the samples in a matrix.
		/*!
		Each row contains one sample.
		*/
		MatrixView<dreal> data() const
		{
			return data_;
		}

		auto matrix() const
		{
			return asMatrix(data());
		}

		//! Sets the time position of the first sample.
		void setT(integer t)
		{
			t_ = t;
		}

		//! Returns the time position of the first sample.
		integer t() const
		{
			return t_;
		}

		//! Returns an iterator to the beginning of the point set.
		/*!
		Here the signal data is interpreted so that each sample 
		is a point which can be accessed by the Array_PointPolicy<dreal>
		point policy. See 'pastel/sys/pointpolicy.txt'.
		*/
		ranges::random_access_range auto pointRange(integer dimensionBegin = 0) {
			return sparseRange(intervalRange(
					data_.data() + dimensionBegin, data_.data() + dimensionBegin + samples() * dimension()), 
				dimension());
		}

		//! Returns an iterator to the beginning of the point set.
		/*!
		Here the signal data is interpreted so that each sample 
		is a point which can be accessed by the Array_PointPolicy<dreal>
		point policy. See 'pastel/sys/pointpolicy.txt'.
		*/
		ranges::random_access_range auto pointRange(integer dimensionBegin = 0) const {
			return sparseRange(intervalRange(
					data_.data() + dimensionBegin, data_.data() + dimensionBegin + samples() * dimension()), 
				dimension());
		}
	private:
		//! The signal data.
		/*!
		Each row is one sample, and thus the width
		of this matrix is the dimensionality. The number
		of samples is given by the number of rows.
		*/
		MatrixView<dreal> data_;
		
		//! The time position of the first sample.
		integer t_ = 0;
	};

	class SignalData {
	public:
		SignalData() = default;
		
		SignalData(integer dimension, integer samples, integer t = 0) 
		: data_(dimension, samples)
		, t_(t)
		{}

		SignalData(const SignalData& that) = default;
		SignalData(SignalData&& that) = default;
		SignalData& operator=(const SignalData& that) = default;
		SignalData& operator=(SignalData&& that) = default;

		operator Signal() const {
			return data();
		}

		integer samples() const {
			return data_.cols();
		}

		integer dimension() const {
			return data_.rows();
		}

		integer t() const {
			return t_;
		}

		MatrixView<dreal> data() const {
			return removeConst(data_).view();
		}

	private:
		MatrixData<dreal> data_;
		integer t_ = 0;
	};

}

#endif
