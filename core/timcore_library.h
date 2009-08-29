// Description: Dll-keywords for the TimCore library

#ifndef TIM_TIMCORE_LIBRARY_H
#define TIM_TIMCORE_LIBRARY_H

#include <pastel/sys/environment.h>

#if defined TIMCORE_EXPORTS
#   define TIMCORE PASTEL_DLLEXPORT
#else
#   define TIMCORE PASTEL_DLLIMPORT
#endif

#endif
