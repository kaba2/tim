# Description: Pastel configuration

# PastelDirectory (string):
#    The directory to the Pastel library.
#
# returns
# -------
#
# PastelIncludeDirectory (string):
#    A directory to add to include directories, such that
#    #include <pastel/sys/mytypes.h>
#    becomes valid.
#
# PastelLibraryPath (string):
#    Path to the Pastel library.
#
# PastelLibraryFilename (string):
#    Filename-part of ${PastelLibraryPath} (e.g. libpastelsys.a).
#
# PastelLibraryDirectory (string):
#    Directory-part of ${PastelLibraryPath}.
#
# PastelLibraryName (string):
#    The name of the library (e.g. pastelsys).

string (TOLOWER "${CMAKE_BUILD_TYPE}" LOWER_CMAKE_BUILD_TYPE)

set (PastelIncludeDirectory "${PastelDirectory}")
set (PastelLibraryDirectory "${PastelDirectory}/lib/${ToolSet}/${LOWER_CMAKE_BUILD_TYPE}")

EcCheckPathExists("Pastel (include)" "${PastelIncludeDirectory}")
EcCheckPathExists("Pastel (library)" "${PastelLibraryDirectory}")


