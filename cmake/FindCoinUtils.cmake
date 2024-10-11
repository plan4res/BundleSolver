# --------------------------------------------------------------------------- #
#    CMake find module for CoinUtils                                          #
#                                                                             #
#    This module finds CoinUtils include directories and libraries.           #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(CoinUtils [version] [EXACT] [REQUIRED])                 #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        CoinUtils_FOUND         - True if headers are found                  #
#        CoinUtils_INCLUDE_DIRS  - Include directories                        #
#        CoinUtils_LIBRARIES     - Libraries to be linked                     #
#        CoinUtils_VERSION       - Version number                             #
#                                                                             #
#    The search results are saved in these persistent cache entries:          #
#                                                                             #
#        CoinUtils_INCLUDE_DIR   - Directory containing headers               #
#        CoinUtils_LIBRARY       - The found library                          #
#                                                                             #
#    This module can read a search path from the variable:                    #
#                                                                             #
#        CoinUtils_ROOT          - Preferred CoinUtils location               #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        Coin::CoinUtils                                                      #
#                                                                             #
#    This find module is provided because CoinUtils does not provide          #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
find_package(BZip2 REQUIRED QUIET)
find_package(LAPACK REQUIRED QUIET)

# Check if already in cache
if (CoinUtils_INCLUDE_DIR AND CoinUtils_LIBRARY)
    set(CoinUtils_FOUND TRUE)
else ()

    # ----- Find the library ------------------------------------------------ #
    # Note that find_path() creates a cache entry
    find_path(CoinUtils_INCLUDE_DIR
              NAMES CoinUtilsConfig.h
              HINTS ${CoinUtils_ROOT}/include
              PATH_SUFFIXES coin coinutils/coin coin-or
              DOC "CoinUtils include directory.")

    find_library(CoinUtils_LIBRARY
                 NAMES CoinUtils
                 HINTS ${CoinUtils_ROOT}/lib
                 DOC "CoinUtils library.")

    # ----- Parse the version ----------------------------------------------- #
    if (CoinUtils_INCLUDE_DIR)
        file(STRINGS
                "${CoinUtils_INCLUDE_DIR}/CoinUtilsConfig.h"
                _coinutils_version_lines REGEX "#define COINUTILS_VERSION_(MAJOR|MINOR|RELEASE)")

        string(REGEX REPLACE ".*COINUTILS_VERSION_MAJOR *\([0-9]*\).*" "\\1" _coinutils_version_major "${_coinutils_version_lines}")
        string(REGEX REPLACE ".*COINUTILS_VERSION_MINOR *\([0-9]*\).*" "\\1" _coinutils_version_minor "${_coinutils_version_lines}")
        string(REGEX REPLACE ".*COINUTILS_VERSION_RELEASE *\([0-9]*\).*" "\\1" _coinutils_version_release "${_coinutils_version_lines}")

        set(CoinUtils_VERSION "${_coinutils_version_major}.${_coinutils_version_minor}.${_coinutils_version_release}")
        unset(_coinutils_version_lines)
        unset(_coinutils_version_major)
        unset(_coinutils_version_minor)
        unset(_coinutils_version_release)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            CoinUtils
            REQUIRED_VARS CoinUtils_LIBRARY CoinUtils_INCLUDE_DIR
            VERSION_VAR CoinUtils_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (CoinUtils_FOUND)
    set(CoinUtils_INCLUDE_DIRS "${CoinUtils_INCLUDE_DIR}")
    set(CoinUtils_LIBRARIES "${CoinUtils_LIBRARY}")

    if (NOT TARGET Coin::CoinUtils)
        add_library(Coin::CoinUtils UNKNOWN IMPORTED)
        set_target_properties(
                Coin::CoinUtils PROPERTIES
                IMPORTED_LOCATION "${CoinUtils_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${CoinUtils_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "BZip2::BZip2")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(CoinUtils_INCLUDE_DIR
                 CoinUtils_LIBRARY
                 CoinUtils_VERSION)

# --------------------------------------------------------------------------- #
