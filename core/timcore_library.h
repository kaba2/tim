#ifndef TIMCORE_TIMCORE_LIBRARY_H
#define TIMCORE_TIMCORE_LIBRARY_H

#include <pastel/sys/environment.h>

#if defined TIMCORE_EXPORTS
#   define TIMCORE PASTEL_DLLEXPORT
#else
#   define TIMCORE PASTEL_DLLIMPORT
#endif

#endif
