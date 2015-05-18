#ifndef ABFAdminH
#define ABFAdminH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// ===============================================================================

#include "ABFAdmOptions.hpp"
#include <AdminClient.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>

//------------------------------------------------------------------------------

class CABFAdmin : private CAdminClient {
public:
    // constructor
    CABFAdmin(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CABFAdmOptions      Options;            // program options
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// get ABF accumulator - this method is specific
    bool GetABFAccumulator(void);
};

//------------------------------------------------------------------------------

#endif
