#ifndef PMFOperationH
#define PMFOperationH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2006      Petr Kulhanek, kulhanek@chemi.muni.cz
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
// =============================================================================

#include <PMFMainHeader.hpp>
#include <Operation.hpp>

//------------------------------------------------------------------------------

/// get initial data
DECLARE_OPERATION(OperationPMF_GetInitialData);

/// exchange data
DECLARE_OPERATION(OperationPMF_ExchangeData);

/// get data
DECLARE_OPERATION(OperationPMF_GetData);

/// get data
DECLARE_OPERATION(OperationPMF_FlushServerData);

//------------------------------------------------------------------------------

#endif
