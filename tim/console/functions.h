#ifndef TIM_FUNCTIONS_H
#define TIM_FUNCTIONS_H

#include "tim/console/scanner.h"

#include <boost/any.hpp>

namespace Tim
{

	boost::any* differential_entropy_kl(const YYLTYPE& location, const AnySet& argSet);
	boost::any* differential_entropy_kl_t(const YYLTYPE& location, const AnySet& argSet);

	boost::any* differential_entropy_nk(const YYLTYPE& location, const AnySet& argSet);

	boost::any* mutual_information_t(const YYLTYPE& location, const AnySet& argSet);
	boost::any* mutual_information(const YYLTYPE& location, const AnySet& argSet);
	boost::any* mutual_information_pt(const YYLTYPE& location, const AnySet& argSet);
	boost::any* mutual_information_p(const YYLTYPE& location, const AnySet& argSet);

	boost::any* transfer_entropy_t(const YYLTYPE& location, const AnySet& argSet);
	boost::any* transfer_entropy_pt(const YYLTYPE& location, const AnySet& argSet);
	boost::any* transfer_entropy(const YYLTYPE& location, const AnySet& argSet);
	boost::any* transfer_entropy_p(const YYLTYPE& location, const AnySet& argSet);

	boost::any* divergence_wkv(const YYLTYPE& location, const AnySet& argSet);

}

#endif
