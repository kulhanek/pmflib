#ifndef PMFMainHeaderH
#define PMFMainHeaderH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301  USA
// ===============================================================================

/** \defgroup   CPMFLib
    \brief      C++ part of PMFLibrary
*/
/** \file       PMFMainHeader.hpp
*   \brief      Main header file for CPMF library.
*   \ingroup    CPMFLib
*/

//------------------------------------------------------------------------------

extern const char* LibBuildVersion_PMF;

//------------------------------------------------------------------------------
// 32bit vs 64bit integers

#ifdef HAVE_MKL_ILP64
typedef long int FTINT;
typedef unsigned long int UFTINT;
#else
typedef int FTINT;
typedef unsigned int UFTINT;
#endif

//------------------------------------------------------------------------------

#if defined _WIN32 || defined __CYGWIN__
#ifdef PMFLIB_BUILDING_DLL
#ifdef __GNUC__
#define PMFLIB_DLL_PUBLIC __attribute__((dllexport))
#else
#define PMFLIB_DLL_PUBLIC __declspec(dllexport)
#endif
#else
#ifdef __GNUC__
#define PMFLIB_DLL_PUBLIC __attribute__((dllimport))
#else
#define PMFLIB_DLL_PUBLIC __declspec(dllimport)
#endif
#define PMFLIB_DLL_LOCAL
#endif
#else
#if __GNUC__ >= 4
#define PMFLIB_DLL_PUBLIC __attribute__ ((visibility("default")))
#define PMFLIB_DLL_LOCAL  __attribute__ ((visibility("hidden")))
#else
#define PMFLIB_DLL_PUBLIC
#define PMFLIB_DLL_LOCAL
#endif
#endif

#define PMF_PACKAGE PMFLIB_DLL_PUBLIC

//------------------------------------------------------------------------------

#endif
