#ifndef LauncherH
#define LauncherH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <SmartThread.hpp>
#include <vector>
#include <PrmFile.hpp>
#include <ostream>
#include <VerboseStr.hpp>
#include <SmallString.hpp>

//------------------------------------------------------------------------------

class CLauncherJob {
public:
    CLauncherJob(void);
public:
    CSmallString    HostName;
    CSmallString    JobDir;
    CSmallString    JobName;
    int             BeadID;
    CSmallString    JobID;
    int             SerialID;   // serial id used for psubmit suffix
    bool            Submitted;  // job was submitted
};

//------------------------------------------------------------------------------

class CLauncher : public CSmartThread {
public:
    // constructor
    CLauncher(void);

    /// read control file
    bool ReadControl(CPrmFile& confile,std::ostream& vout);

    /// is enabled?
    bool IsEnabled(void);

    /// start launcher
    bool StartLauncher(std::ostream& vout);

// section of private data -----------------------------------------------------
private:
    // setup
    bool                        Enabled;

    // output stream
    CSmallString                LogFile;
    std::ofstream               lout;
    CSmallString                JobFile;
    std::vector<CLauncherJob>   Jobs;

    // wrappers
    CSmallString    DistributeKeyWrapper;
    CSmallString    SubmitJobWrapper;
    CSmallString    JobStatusWrapper;

    // sleep times
    int             DistributeKeySleepTime;
    int             CheckServerSleepTime;

// read controls ------------------------
    /// read luncher [setup]
    bool ReadSetup(CPrmFile& confile,std::ostream& vout);

    /// read [wrappers]
    bool ReadWrappers(CPrmFile& confile,std::ostream& vout);

    /// read job specification from {JOBS} or individual file {MAIN}
    bool ReadJobs(CPrmFile& confile,std::ostream& vout);

// main execution
    /// main launcher thread
    virtual void ExecuteThread(void);

// specific actions----------------------
    /// distribute server key to all clients
    bool DistributeKey(void);

    /// distribute key for given job
    bool DistributeKeyForJob(const CLauncherJob& job);

    /// submit all client jobs until rendezvous
    bool SubmitAllJobsAndWaitForRendezvous(void);

    /// submit all client jobs
    bool SubmitAllJobs(void);

    /// wait for all jobs
    bool WaitForAllJobs(void);

    /// wait for all jobs
    int WaitForJobs(void);

    /// submit individual job
    bool SubmitJob(CLauncherJob& job,CSmallString& id);

    /// is it job terminated?
    bool IsJobFinished(CLauncherJob& job);
};

//------------------------------------------------------------------------------

#endif
