// Description: tim_matlab() function registration
// Documentation: tim_matlab_cpp.txt

#ifndef TIM_TIM_MATLAB_H
#define TIM_TIM_MATLAB_H

#include "tim/matlab/tim_mex.h"

namespace Tim
{

	//! Registers a function to be callable from Matlab.
	/*!
	The functions that are registered with this function 
	are made callable via the tim_matlab() mex function
	(whose entry point is located in tim_matlab.cpp).
	*/
	void matlabAddFunction(
		const std::string& name,
		MatlabFunction* function);

}

#endif
