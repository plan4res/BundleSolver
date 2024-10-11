# --------------------------------------------------------------------------- #
#    CMake find module for Clp                                                #
#                                                                             #
#    This module finds Clp include directories and libraries.                 #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(Clp [version] [EXACT] [REQUIRED])                       #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        Clp_FOUND         - True if headers are found                        #
#        Clp_INCLUDE_DIRS  - Include directories                              #
#        Clp_LIBRARIES     - Libraries to be linked                           #
#        Clp_VERSION       - Version number                                   #
#                                                                             #
#    The search results are saved in these persistent cache entries:          #
#                                                                             #
#        Clp_INCLUDE_DIR   - Directory containing headers                     #
#        Clp_LIBRARY       - The found library                                #
#                                                                             #
#    This module can read a search path from the variable:                    #
#                                                                             #
#        Clp_ROOT          - Preferred Clp location                           #
#                                                                             #
#    The following IMPORTED targets are also defined:                         #
#                                                                             #
#        Coin::Clp                                                            #
#        Coin::ClpSolver                                                      #
#        Coin::OsiClp                                                         #
#                                                                             #
#    This find module is provided because Clp does not provide                #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
find_package(CoinUtils REQUIRED)

# ----- Find the library ---------------------------------------------------- #
# Note that find_path() creates a cache entry
find_path(Clp_INCLUDE_DIR
          NAMES ClpConfig.h
          HINTS ${Clp_ROOT}/include
          PATH_SUFFIXES coin clp/coin coin-or
          DOC "Clp include directory.")

find_library(Clp_LIBRARY
             NAMES Clp
             HINTS ${Clp_ROOT}/lib
             DOC "Clp library.")

# ----- ClpSolver component ------------------------------------------------- #
# Note that find_path() creates a cache entry
find_path(Clp_ClpSolver_INCLUDE_DIR
          NAMES ClpSolve.hpp
          HINTS ${Clp_ROOT}/include
          PATH_SUFFIXES coin clp/coin coin-or
          DOC "ClpSolver include directory.")

find_library(Clp_ClpSolver_LIBRARY
             NAMES ClpSolver
             HINTS ${Clp_ROOT}/lib
             DOC "ClpSolver library.")

if (Clp_ClpSolver_INCLUDE_DIR AND Clp_ClpSolver_LIBRARY)
    set(Clp_ClpSolver_FOUND TRUE)
endif ()

# ----- OsiClp component ---------------------------------------------------- #
# Note that find_path() creates a cache entry
find_path(Clp_OsiClp_INCLUDE_DIR
          NAMES OsiClpSolverInterface.hpp
          HINTS ${Clp_ROOT}/include
          PATH_SUFFIXES coin clp/coin coin-or
          DOC "OsiClp include directory.")

find_library(Clp_OsiClp_LIBRARY
             NAMES OsiClp
             HINTS ${Clp_ROOT}/lib
             DOC "OsiClp library.")

if (Clp_OsiClp_INCLUDE_DIR AND Clp_OsiClp_LIBRARY)
    set(Clp_OsiClp_FOUND TRUE)
endif ()

# ----- Parse the version --------------------------------------------------- #
if (Clp_INCLUDE_DIR)
    file(STRINGS
            "${Clp_INCLUDE_DIR}/ClpConfig.h"
            _clp_version_lines REGEX "#define CLP_VERSION_(MAJOR|MINOR|RELEASE)")

    string(REGEX REPLACE ".*CLP_VERSION_MAJOR *\([0-9]*\).*" "\\1" _clp_version_major "${_clp_version_lines}")
    string(REGEX REPLACE ".*CLP_VERSION_MINOR *\([0-9]*\).*" "\\1" _clp_version_minor "${_clp_version_lines}")
    string(REGEX REPLACE ".*CLP_VERSION_RELEASE *\([0-9]*\).*" "\\1" _clp_version_release "${_clp_version_lines}")

    set(Clp_VERSION "${_clp_version_major}.${_clp_version_minor}.${_clp_version_release}")
    unset(_clp_version_lines)
    unset(_clp_version_major)
    unset(_clp_version_minor)
    unset(_clp_version_release)
endif ()

# ----- Handle the standard arguments --------------------------------------- #
# The following macro manages the QUIET, REQUIRED and version-related
# options passed to find_package(). It also sets <PackageName>_FOUND if
# REQUIRED_VARS are set.
# REQUIRED_VARS should be cache entries and not output variables. See:
# https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
find_package_handle_standard_args(
        Clp
        REQUIRED_VARS Clp_LIBRARY Clp_INCLUDE_DIR
        VERSION_VAR Clp_VERSION
        HANDLE_COMPONENTS)

# ----- Export the targets -------------------------------------------------- #
if (Clp_FOUND)
    set(Clp_INCLUDE_DIRS "${Clp_INCLUDE_DIR}")
    set(Clp_LIBRARIES "${Clp_LIBRARY}")

    if (NOT TARGET Coin::Clp)
        add_library(Coin::Clp UNKNOWN IMPORTED)
        set_target_properties(
                Coin::Clp PROPERTIES
                IMPORTED_LOCATION "${Clp_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${Clp_INCLUDE_DIR}")
    endif ()
endif ()

if (Clp_ClpSolver_FOUND)
    set(Clp_ClpSolver_INCLUDE_DIRS "${Clp_ClpSolver_INCLUDE_DIR}")
    set(Clp_ClpSolver_LIBRARIES "${Clp_ClpSolver_LIBRARY}")

    if (NOT TARGET Coin::ClpSolver)
        add_library(Coin::ClpSolver UNKNOWN IMPORTED)
        set_target_properties(
                Coin::ClpSolver PROPERTIES
                IMPORTED_LOCATION "${Clp_ClpSolver_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${Clp_ClpSolver_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "Coin::Clp")
    endif ()
endif ()

if (Clp_OsiClp_FOUND)
    set(Clp_OsiClp_INCLUDE_DIRS "${Clp_OsiClp_INCLUDE_DIR}")
    set(Clp_OsiClp_LIBRARIES "${Clp_OsiClp_LIBRARY}")

    if (NOT TARGET Coin::OsiClp)
        add_library(Coin::OsiClp UNKNOWN IMPORTED)
        set_target_properties(
                Coin::OsiClp PROPERTIES
                IMPORTED_LOCATION "${Clp_OsiClp_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${Clp_OsiClp_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "Coin::Clp;Coin::Osi")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(Clp_INCLUDE_DIR Clp_ClpSolver_INCLUDE_DIR Clp_OsiClp_INCLUDE_DIR
                 Clp_LIBRARY Clp_ClpSolver_LIBRARY Clp_OsiClp_LIBRARY
                 Clp_VERSION)

# --------------------------------------------------------------------------- #
