// Description: Dll-keywords for the TIM library
// Documentation: basics.txt

#ifndef TIM_TIM_LIBRARY_H
#define TIM_TIM_LIBRARY_H

#include <pastel/sys/environment.h>

#ifdef _USRDLL
#	ifdef TIM_EXPORTS
#		define TIM PASTEL_DLLEXPORT
#	else
#		define TIM PASTEL_DLLIMPORT
#	endif
#else
#	define TIM
#endif

#endif
