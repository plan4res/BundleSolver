# --------------------------------------------------------------------------- #
#    CMake find module for Open Solver Interface                              #
#                                                                             #
#    This module finds Osi include directories and libraries.                 #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(Osi [version] [EXACT] [REQUIRED])                       #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        Osi_FOUND         - True if headers are found                        #
#        Osi_INCLUDE_DIRS  - Include directories                              #
#        Osi_LIBRARIES     - Libraries to be linked                           #
#        Osi_VERSION       - Version number                                   #
#                                                                             #
#    The search results are saved in these persistent cache entries:          #
#                                                                             #
#        Osi_INCLUDE_DIR   - Directory containing headers                     #
#        Osi_LIBRARY       - The found library                                #
#                                                                             #
#    This module can read a search path from the variable:                    #
#                                                                             #
#        Osi_ROOT          - Preferred Osi location                           #
#                                                                             #
#    The following IMPORTED targets are also defined:                         #
#                                                                             #
#        Coin::Osi                                                            #
#        Coin::OsiCpx                                                         #
#        Coin::OsiGrb                                                         #
#                                                                             #
#    This find module is provided because Osi does not provide                #
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
find_package(CPLEX)
find_package(GUROBI)

# ----- Find the library ---------------------------------------------------- #
# Note that find_path() creates a cache entry
find_path(Osi_INCLUDE_DIR
          NAMES OsiConfig.h
          HINTS ${Osi_ROOT}/include
          PATH_SUFFIXES coin osi/coin coin-or
          DOC "Osi include directory.")

find_library(Osi_LIBRARY
             NAMES Osi
             HINTS ${Osi_ROOT}/lib
             DOC "Osi library.")

# ----- OsiCpx component ---------------------------------------------------- #
if (CPLEX_FOUND)
    # Note that find_path() creates a cache entry
    find_path(Osi_OsiCpx_INCLUDE_DIR
              NAMES OsiCpxSolverInterface.hpp
              HINTS ${Osi_ROOT}/include
              PATH_SUFFIXES coin osi/coin coin-or
              DOC "OsiCpx include directory.")

    find_library(Osi_OsiCpx_LIBRARY
                 NAMES OsiCpx
                 HINTS ${Osi_ROOT}/lib
                 DOC "OsiCpx library.")

    if (Osi_OsiCpx_INCLUDE_DIR AND Osi_OsiCpx_LIBRARY)
        set(Osi_OsiCpx_FOUND TRUE)
    endif ()
endif ()

# ----- OsiGrb component ---------------------------------------------------- #
if (GUROBI_FOUND)
    # Note that find_path() creates a cache entry
    find_path(Osi_OsiGrb_INCLUDE_DIR
              NAMES OsiGrbSolverInterface.hpp
              HINTS ${Osi_ROOT}/include
              PATH_SUFFIXES coin osi/coin coin-or
              DOC "OsiGrb include directory.")

    find_library(Osi_OsiGrb_LIBRARY
                 NAMES OsiGrb
                 HINTS ${Osi_ROOT}/lib
                 DOC "OsiGrb library.")

    if (Osi_OsiGrb_INCLUDE_DIR AND Osi_OsiGrb_LIBRARY)
        set(Osi_OsiGrb_FOUND TRUE)
    endif ()
endif ()

# ----- Parse the version --------------------------------------------------- #
if (Osi_INCLUDE_DIR)
    file(STRINGS
            "${Osi_INCLUDE_DIR}/OsiConfig.h"
            _osi_version_lines REGEX "#define OSI_VERSION_(MAJOR|MINOR|RELEASE)")

    string(REGEX REPLACE ".*OSI_VERSION_MAJOR *\([0-9]*\).*" "\\1" _osi_version_major "${_osi_version_lines}")
    string(REGEX REPLACE ".*OSI_VERSION_MINOR *\([0-9]*\).*" "\\1" _osi_version_minor "${_osi_version_lines}")
    string(REGEX REPLACE ".*OSI_VERSION_RELEASE *\([0-9]*\).*" "\\1" _osi_version_release "${_osi_version_lines}")

    set(Osi_VERSION "${_osi_version_major}.${_osi_version_minor}.${_osi_version_release}")
    unset(_osi_version_lines)
    unset(_osi_version_major)
    unset(_osi_version_minor)
    unset(_osi_version_release)
endif ()

# ----- Handle the standard arguments --------------------------------------- #
# The following macro manages the QUIET, REQUIRED and version-related
# options passed to find_package(). It also sets <PackageName>_FOUND if
# REQUIRED_VARS are set.
# REQUIRED_VARS should be cache entries and not output variables. See:
# https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
find_package_handle_standard_args(
        Osi
        REQUIRED_VARS Osi_LIBRARY Osi_INCLUDE_DIR
        VERSION_VAR Osi_VERSION
        HANDLE_COMPONENTS)

# ----- Export the targets -------------------------------------------------- #
if (Osi_FOUND)
    set(Osi_INCLUDE_DIRS "${Osi_INCLUDE_DIR}")
    set(Osi_LIBRARIES "${Osi_LIBRARY}")

    if (NOT TARGET Coin::Osi)
        add_library(Coin::Osi UNKNOWN IMPORTED)
        set_target_properties(
                Coin::Osi PROPERTIES
                IMPORTED_LOCATION "${Osi_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${Osi_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES};Coin::CoinUtils")
    endif ()
endif ()

if (Osi_OsiCpx_FOUND)
    set(Osi_OsiCpx_INCLUDE_DIRS "${Osi_OsiCpx_INCLUDE_DIR}")
    set(Osi_OsiCpx_LIBRARIES "${Osi_OsiCpx_LIBRARY}")

    if (NOT TARGET Coin::OsiCpx)
        add_library(Coin::OsiCpx UNKNOWN IMPORTED)
        set_target_properties(
                Coin::OsiCpx PROPERTIES
                IMPORTED_LOCATION "${Osi_OsiCpx_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${Osi_OsiCpx_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "Coin::Osi;CPLEX::Cplex")
    endif ()
endif ()

if (Osi_OsiGrb_FOUND)
    set(Osi_OsiGrb_INCLUDE_DIRS "${Osi_OsiGrb_INCLUDE_DIR}")
    set(Osi_OsiGrb_LIBRARIES "${Osi_OsiGrb_LIBRARY}")

    if (NOT TARGET Coin::OsiGrb)
        add_library(Coin::OsiGrb UNKNOWN IMPORTED)
        set_target_properties(
                Coin::OsiGrb PROPERTIES
                IMPORTED_LOCATION "${Osi_OsiGrb_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${Osi_OsiGrb_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "Coin::Osi;GUROBI::Gurobi")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(Osi_INCLUDE_DIR Osi_OsiCpx_INCLUDE_DIR Osi_OsiGrb_INCLUDE_DIR
                 Osi_LIBRARY Osi_OsiCpx_LIBRARY Osi_OsiGrb_LIBRARY
                 Osi_VERSION)

# --------------------------------------------------------------------------- #
