// Description: Built-in functions for the interpreter

#ifndef TIM_FUNCTIONS_H
#define TIM_FUNCTIONS_H

#include "tim/console/console_scanner.h"

#include <string>

#include <boost/any.hpp>

namespace Tim
{

	class FunctionCall_Exception
	{
	};

	boost::any functionCall(const std::string& name, const AnySet& argSet);

	boost::any load(const AnySet& argSet);

	boost::any differential_entropy_kl(const AnySet& argSet);
	boost::any differential_entropy_kl_t(const AnySet& argSet);

	boost::any differential_entropy_nk(const AnySet& argSet);

	boost::any mutual_information_t(const AnySet& argSet);
	boost::any mutual_information(const AnySet& argSet);
	boost::any mutual_information_pt(const AnySet& argSet);
	boost::any mutual_information_p(const AnySet& argSet);

	boost::any transfer_entropy_t(const AnySet& argSet);
	boost::any transfer_entropy_pt(const AnySet& argSet);
	boost::any transfer_entropy(const AnySet& argSet);
	boost::any transfer_entropy_p(const AnySet& argSet);

	boost::any divergence_wkv(const AnySet& argSet);

}

#endif
