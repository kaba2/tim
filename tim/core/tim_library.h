// Description: Dll-keywords for the Tim library

#ifndef TIM_TIM_LIBRARY_H
#define TIM_TIM_LIBRARY_H

#include <pastel/sys/environment.h>

#if defined TIM_EXPORTS
#   define TIM PASTEL_DLLEXPORT
#else
#   define TIM PASTEL_DLLIMPORT
#endif

#endif
