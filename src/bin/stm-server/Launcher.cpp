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

#include <stdio.h>
#include <ErrorSystem.hpp>
#include <PMFOperation.hpp>
#include <ExtraOperation.hpp>
#include "Launcher.hpp"
#include "StringServer.hpp"
#include <ServerCommand.hpp>
#include <iomanip>
#include "StringServer.hpp"
#include <FileSystem.hpp>
#include <errno.h>
#include <unistd.h>

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CLauncherJob::CLauncherJob(void)
{
    BeadID = -1;
    SerialID = 1;
    Submitted = false;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CLauncher::CLauncher(void)
{
    Enabled = false;
    JobFile = "jobs.in";
    LogFile = "launcher.log";

    DistributeKeyWrapper = "distribute-key";
    SubmitJobWrapper = "submit-job";
    JobStatusWrapper = "job-status";

    DistributeKeySleepTime = 5;
    CheckServerSleepTime = 5;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CLauncher::ReadControl(CPrmFile& confile,ostream& vout)
{
    vout << endl;
    vout << ":::::::::::::::::::::::::::::::::::: {LAUNCHER} ::::::::::::::::::::::::::::::::" << endl;

    if( confile.OpenGroup("LAUNCHER") == false ) {
        vout << ">> Info: No {LAUNCHER} group is specified - disabling job launcher ..." << endl;
        Enabled = false;
        return(true);
    }

    Enabled = true;
    vout << ">> Info: Job launcher is enabled ..." << endl;

    // [setup]
    if( ReadSetup(confile,vout) == false ){
        ES_TRACE_ERROR("unable to read setup");
        return(false);
    }

    // [wrappers]
    if( ReadWrappers(confile,vout) == false ){
        ES_TRACE_ERROR("unable to read wrappers");
        return(false);
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::::: {JOBS} :::::::::::::::::::::::::::::::::::" << endl;

    // read jobs
    if( (JobFile.GetLength() > 2) && (JobFile[0] == '{') && (JobFile[JobFile.GetLength()-1] == '}') ){
        CSmallString group = JobFile.GetSubStringFromTo(1,JobFile.GetLength()-2);
        vout << ">> Info: Job specification is read from the {" << group << "} group." << endl;
        if( confile.OpenGroup(group) == false ) {
            CSmallString error;
            error << "unable to open group '" << group << "' with job specification for reading";
            ES_ERROR(error);
            return(false);
        }
        if( ReadJobs(confile,vout) == false ){
            ES_TRACE_ERROR("unable to read jobs");
            return(false);
        }
    } else {
        vout << ">> Info: Job specification is read from the '" << JobFile << "' file." << endl;
        CPrmFile prmfile;
        if( prmfile.Read(JobFile) == false) {
            CSmallString error;
            error << "unable to load job specification file '" << JobFile << "'";
            ES_ERROR(error);
            return(false);
        }
        if( ReadJobs(prmfile,vout) == false ){
            ES_TRACE_ERROR("unable to read jobs");
            return(false);
        }
        if( prmfile.CountULines() > 0 ){
            vout << endl;
            vout << "::::::::::::::::::::::::::::::: Unprocessed items ::::::::::::::::::::::::::::::" << endl;
            ES_ERROR("unprocessed items found in job specification file");
            prmfile.Dump(stderr,true);
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CLauncher::ReadSetup(CPrmFile& confile,ostream& vout)
{
    vout << endl;
    vout << "=== [setup] ====================================================================" << endl;
    if( confile.OpenSection("setup") == false ) {
        vout << "Job specification file name (jobs)             = " << left << setw(20) << JobFile << "  (default)" << endl;
        vout << "Log file name (log)                            = " << left << setw(20) << LogFile << "  (default)" << endl;
        return(true);
    }

    if( confile.GetStringByKey("jobs",JobFile) == true  ) {
        vout << "Job specification file name (jobs)             = " << left << setw(20) << JobFile << endl;
    } else {
        vout << "Job specification file name (jobs)             = " << left << setw(20) << JobFile << "  (default)" << endl;
    }

    if( confile.GetStringByKey("log",LogFile) == true  ) {
        vout << "Log file name (log)                            = " << left << setw(20) << LogFile << endl;
    } else {
        vout << "Log file name (log)                            = " << left << setw(20) << LogFile << "  (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CLauncher::ReadWrappers(CPrmFile& confile,ostream& vout)
{
    vout << endl;
    vout << "=== [wrappers] =================================================================" << endl;
    if( confile.OpenSection("wrappers") == false ) {
        vout << "Distribute key wrapper (distribute-key)        = " << left << setw(20) << DistributeKeyWrapper << "  (default)" << endl;
        vout << "Submit job wrapper (submit-job)                = " << left << setw(20) << SubmitJobWrapper << "  (default)" << endl;
        vout << "Job status wrapper (job-status)                = " << left << setw(20) << JobStatusWrapper << "  (default)" << endl;
        return(true);
    }

    if( confile.GetStringByKey("distribute-key",DistributeKeyWrapper) == true  ) {
        vout << "Distribute key wrapper (distribute-key)        = " << left << setw(20) << DistributeKeyWrapper << endl;
    } else {
        vout << "Distribute key wrapper (distribute-key)        = " << left << setw(20) << DistributeKeyWrapper << "  (default)" << endl;
    }

    if( confile.GetStringByKey("submit-job",SubmitJobWrapper) == true  ) {
        vout << "Submit job wrapper (submit-job)                = " << left << setw(20) << SubmitJobWrapper << endl;
    } else {
        vout << "Submit job wrapper (submit-job)                = " << left << setw(20) << SubmitJobWrapper << "  (default)" << endl;
    }

    if( confile.GetStringByKey("job-status",JobStatusWrapper) == true  ) {
        vout << "Job status wrapper (job-status)                = " << left << setw(20) << JobStatusWrapper << endl;
    } else {
        vout << "Job status wrapper (job-status)                = " << left << setw(20) << JobStatusWrapper << "  (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CLauncher::ReadJobs(CPrmFile& confile,ostream& vout)
{
    CSmallString default_hostname;
    CSmallString default_name;

    // read default setup
    vout << endl;
    vout << "=== [default] ==================================================================" << endl;
    if( confile.OpenSection("default") == false ) {
        vout << ">> Info: No default hostname and job name defined." << endl;
    } else {
        if( confile.GetStringByKey("hostname",default_hostname) == true  ) {
            vout << "Default job hostname (hostname)            = " << left << default_hostname << endl;
        }

        if( confile.GetStringByKey("jobname",default_name) == true  ) {
            vout << "Default job name (jobname)                 = " << left << default_name << endl;
        }
    }

    // count number of jobs
    int count = 0;
    confile.FirstSection();
    do {
        if( confile.GetSectionName() == "job" ){
            count++;
        }

    } while( confile.NextSection() );

    vout << endl;
    vout << "Number of job specifications = " << count << endl;

    if( count == 0 ){
        ES_ERROR("no job is specified");
        return(false);
    }

    vout << endl;
    vout << "# BeadID     Name          HostName     Path                                    " << endl;
    vout << "# ------ ------------- ---------------- ----------------------------------------" << endl;

    confile.FirstSection();
    int id = 1;
    do {
        if( confile.GetSectionName() == "job" ){
            CLauncherJob job;
            if( ! default_hostname.IsEmpty() ){
                job.HostName = default_hostname;
            } else {
                if( confile.GetStringByKey("hostname",job.HostName) == false  ) {
                    CSmallString error;
                    error << "default hostname is not specified thus it must be explicitly provided for job id " << id;
                    ES_ERROR(error);
                }
            }
            if( ! default_name.IsEmpty() ){
                job.JobName = default_name;
            } else {
                if( confile.GetStringByKey("jobname",job.JobName) == false  ) {
                    CSmallString error;
                    error << "default job name is not specified thus it must be explicitly provided for job id " << id;
                    ES_ERROR(error);
                }
            }
            if( confile.GetStringByKey("path",job.JobDir) == false  ) {
                CSmallString error;
                error << "job path must be explicitly provided for job id " << id;
                ES_ERROR(error);
            }
            if( confile.GetIntegerByKey("bead_id",job.BeadID) == false  ) {
                CSmallString error;
                error << "bead ID must be explicitly provided for job id " << id;
                ES_ERROR(error);
            }

            Jobs.push_back(job);
            vout << right << setw(8) << job.BeadID << " ";
            vout << left << setw(13) << job.JobName << " ";
            vout << left << setw(16) << job.HostName << " ";
            vout << left << job.JobDir << endl;
            id++;
        }

    } while( confile.NextSection() );

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CLauncher::IsEnabled(void)
{
    return(Enabled);
}

//------------------------------------------------------------------------------

bool CLauncher::StartLauncher(ostream& vout)
{
    vout << endl;
    vout << "Starting job launcher server ... " << endl;

    // do we have wrappers
    CFileName cwd;
    CFileSystem::GetCurrentDir(cwd);

    if( CFileSystem::IsFile(cwd / DistributeKeyWrapper ) == false ){
        CSmallString error;
        error << "distribute-key wrapper (" << DistributeKeyWrapper << ") was not found";
        ES_ERROR(error);
        return(false);
    }
    CFileSystem::SetPosixMode(DistributeKeyWrapper,0744);

    if( CFileSystem::IsFile(cwd / SubmitJobWrapper ) == false ){
        CSmallString error;
        error << "submit-job wrapper (" << SubmitJobWrapper << ") was not found";
        ES_ERROR(error);
        return(false);
    }
    CFileSystem::SetPosixMode(SubmitJobWrapper,0744);

    if( CFileSystem::IsFile(cwd / JobStatusWrapper ) == false ){
        CSmallString error;
        error << "job-status wrapper (" << JobStatusWrapper << ") was not found";
        ES_ERROR(error);
        return(false);
    }
    CFileSystem::SetPosixMode(JobStatusWrapper,0744);

    // open log file
    lout.open(LogFile);
    if( ! lout ){
        CSmallString error;
        error << "unable to open launcher server log file '" << LogFile << "'";
        ES_ERROR(error);
        return(false);
    }

    if( StartThread() == false ) return(false);

    return(true);
}

//------------------------------------------------------------------------------

void CLauncher::ExecuteThread(void)
{
    lout << "# ------------------------------------------------------------------------------" << endl;
    lout << "# LAUNCHER SERVER STARTED" << endl;
    lout << "# ------------------------------------------------------------------------------" << endl;

// distribute server key
    lout << endl;
    lout << "1) Distributing the server key to all clients ..." << endl;
    if( DistributeKey() == false ){
        ES_ERROR("unable to distribute server key");
        StringServer.TerminateServer();
        return;
    }

    StringServer.Beads.BeginAsynchronousMode();

    int iter = 1;
    lout << endl;
    lout << "# ::::::::::::::::::::::::::::::: ENTERING STM LOOP ::::::::::::::::::::::::::::" << endl;

    for(;;){
        lout << endl;
        lout << "##### ITERATION: " << iter << endl;
        if( ThreadTerminated ){
            StringServer.TerminateServer();
            return;
        }

        lout << endl;
        lout << "2) Submitting client jobs until they all are in rendezvous state ..." << endl;
        if( SubmitAllJobsAndWaitForRendezvous() == false ){
            ES_ERROR("unable to submit all jobs");
            StringServer.TerminateServer();
            return;
        }

        // process path data
        lout << endl;
        lout << "3) Processing STM path data ..." << endl;
        StringServer.Beads.ProcessPathAsynchronously();

        lout << endl;
        lout << "4) Wait until all jobs finish ..." << endl;
        if( WaitForAllJobs() == false ){
            ES_ERROR("unable to wait for all jobs");
            StringServer.TerminateServer();
            return;
        }

        if( (StringServer.Beads.GetSTMStatus() == ESTMS_MAX_STEPS_REACHED) ||
            (StringServer.Beads.GetSTMStatus() == ESTMS_COMPLETED) ) {
            break; // exit STM loop
        }

        iter++;
    }

    lout << endl;
    lout << "5) Terminating STM server ..." << endl;
    StringServer.TerminateServer();

    lout << endl;
    lout << "# ------------------------------------------------------------------------------" << endl;
    lout << "# LAUNCHER SERVER FINISHED" << endl;
    lout << "# ------------------------------------------------------------------------------" << endl;
}

//------------------------------------------------------------------------------

bool CLauncher::DistributeKey(void)
{
    // wait for serverkey
    lout << "   Waiting until the server key is ready ..." << endl;
    lout << "       Testing every " << DistributeKeySleepTime << " seconds" << endl;
    while( (StringServer.IsServerKeyReady() == false) && (ThreadTerminated == false) ){
        sleep(DistributeKeySleepTime);
    }
    if( ThreadTerminated ) {
        lout << ">>> INFO: Terminated upon external request ..." << endl;
        return(false);
    }
    lout << "   The server key is ready." << endl;

    // distribute key
    lout << endl;
    lout << "   Now distributing the key by the " << DistributeKeyWrapper << " wrapper ..." << endl;

    vector<CLauncherJob>::iterator  it = Jobs.begin();
    vector<CLauncherJob>::iterator  ie = Jobs.end();

    int id = 1;
    while( (it != ie) && (ThreadTerminated == false)  ){
        lout << "        # " << setfill('0') << setw(3) << id << setfill(' ') << " job ...";
        lout.flush();
        if( DistributeKeyForJob(*it) == false ){
            lout << "FAILED" << endl;
            return(false);
        } else {
            lout << " done." << endl;
        }
        it++;
        id++;
    }
    if( ThreadTerminated ) {
        lout << ">>> INFO: Terminated upon external request ..." << endl;
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CLauncher::DistributeKeyForJob(const CLauncherJob& job)
{
    // assemble command
    CFileName cmd;
    CFileSystem::GetCurrentDir(cmd);
    cmd = cmd / DistributeKeyWrapper
        + " " + StringServer.GetServerKeyName()
        + " " + job.HostName
        + " " + job.JobDir;

    stringstream str;

    FILE* p_file = popen(cmd,"r");
    if( p_file == NULL ){
        lout << ">>> WRAPPER ERROR: " << endl;
        lout << "popen failed, errno(" << strerror(errno) <<  ")" << endl;
        return( false );
    }

    // copy output to stream
    char line[100];
    while( fgets(line, 99, p_file) ){
        str << line;
    }
    int retcode = pclose(p_file);

    if( retcode != 0 ){
        lout << endl;
        lout << ">>> WRAPPER ERROR OUTPUT: " << endl;
        lout << str.str(); // show output only if an error
        lout << "Return code: " << retcode;
        if( retcode < 0 ){
            lout << " errno(" << strerror(errno) <<  ")";
        }
        lout << endl;
    }

    return( retcode == 0 );
}

//------------------------------------------------------------------------------

bool CLauncher::SubmitAllJobsAndWaitForRendezvous(void)
{
    lout << "   Checking STM server/client status every " << CheckServerSleepTime << " seconds." << endl;
    SubmitAllJobs();
    while( (StringServer.Beads.GetNumOfBeadsInRendezvousState() != StringServer.Beads.GetNumOfBeads())
        && (ThreadTerminated == false) ){
        sleep(CheckServerSleepTime);
        SubmitAllJobs();
    }
    if( ThreadTerminated ) {
        lout << ">>> INFO: Terminated upon external request ..." << endl;
        return(false);
    }

    lout << "   All beads are in rendezvous state ..." << endl;

    return(true);
}

//------------------------------------------------------------------------------

bool CLauncher::SubmitAllJobs(void)
{
    // submit individual jobs
    vector<CLauncherJob>::iterator  it = Jobs.begin();
    vector<CLauncherJob>::iterator  ie = Jobs.end();

    while( (it != ie) && (ThreadTerminated == false)  ){
        CLauncherJob& job = *it;
        CBead* p_bead = StringServer.Beads.GetBead(job.BeadID);
        it++;

        if( (p_bead->GetModeStatus() == BMS_PREPARED) &&
            ( (p_bead->GetMode() == BMO_INITIALIZATION) ||
              (p_bead->GetMode() == BMO_EQUILIBRATION) ||
              (p_bead->GetMode() == BMO_ACCUMULATION) ||
              (p_bead->GetMode() == BMO_PRODUCTION) ) ) {
            if( job.Submitted ) continue; // already submitted
            // submit new job
            lout << "      # " << setfill('0') << setw(3) << job.BeadID << setfill(' ') << " bead ";
            lout << "(" << p_bead->GetModeString() << ") ... ";
            CSmallString job_id;
            if( SubmitJob(job,job_id) == false ){
                lout << "FAILED" << endl;
                job.Submitted = false; // we will try it later
            } else {
                lout << job_id << endl;
                job.JobID = job_id;
                job.Submitted = true;
            }
        }

        if( (p_bead->GetModeStatus() == BMS_FINISHED) &&
            ( (p_bead->GetMode() == BMO_INITIALIZATION) ||
              (p_bead->GetMode() == BMO_EQUILIBRATION) ||
              (p_bead->GetMode() == BMO_ACCUMULATION) ||
              (p_bead->GetMode() == BMO_PRODUCTION) ) ){
            if( job.Submitted == false ) continue; // process only submitted jobs
            // is job finished?
            if( IsJobFinished(job) == true ){
                lout << "      # " << setfill('0') << setw(3) << job.BeadID << setfill(' ') << " bead ";
                lout << "(" << p_bead->GetModeString() << ") ... finished" << endl;
                p_bead->MoveToNextMode();
                job.Submitted = false;
            }
        }
    }
    if( ThreadTerminated ) return(false);

    return(true);
}

//------------------------------------------------------------------------------

bool CLauncher::WaitForAllJobs(void)
{
    lout << "   Checking job status every " << CheckServerSleepTime << " seconds." << endl;
    while( (WaitForJobs() > 0) && (ThreadTerminated == false) ){
        sleep(CheckServerSleepTime);
    }
    if( ThreadTerminated ) {
        lout << ">>> INFO: Terminated upon external request ..." << endl;
        return(false);
    }

    lout << "   All bead clients were finished." << endl;

    return(true);
}

//------------------------------------------------------------------------------

int CLauncher::WaitForJobs(void)
{
    // submit individual jobs
    vector<CLauncherJob>::iterator  it = Jobs.begin();
    vector<CLauncherJob>::iterator  ie = Jobs.end();

    int count = 0;

    while( (it != ie) && (ThreadTerminated == false)  ){
        CLauncherJob& job = *it;
        it++;

        if( job.Submitted == false ){
            continue; // process only submitted jobs
        }

        // is job finished?
        if( IsJobFinished(job) == true ){
            lout << "      # " << setfill('0') << setw(3) << job.BeadID << setfill(' ') << " bead";
            lout << " ... finished" << endl;
            job.Submitted = false;
        } else {
            count++;
        }

    }
    if( ThreadTerminated ) return(-1);

    return(count);
}

//------------------------------------------------------------------------------

bool CLauncher::SubmitJob(CLauncherJob& job,CSmallString& id)
{
    // assemble command
    CFileName cmd;
    CFileSystem::GetCurrentDir(cmd);
    cmd = cmd / SubmitJobWrapper
        + " " + job.HostName
        + " " + job.JobDir
        + " " + job.JobName
        + " " + CSmallString(job.SerialID);

    stringstream str;

    FILE* p_file = popen(cmd,"r");
    if( p_file == NULL ){
        lout << ">>> WRAPPER ERROR: " << endl;
        lout << "popen failed, errno(" << strerror(errno) <<  ")" << endl;
        return( false );
    }

    // copy output to stream
    char line[100];
    while( fgets(line, 99, p_file) ){
        str << line;
    }
    int retcode = pclose(p_file);

    if( retcode != 0 ){
        lout << endl;
        lout << ">>> WRAPPER ERROR OUTPUT: " << endl;
        lout << str.str(); // show output only if an error
        lout << "Return code: " << retcode;
        if( retcode < 0 ){
            lout << " errno(" << strerror(errno) <<  ")";
        }
        lout << endl;
        return(false);
    }

    // extract job id
    stringstream idstr(str.str());
    string sid;
    getline(idstr,sid);
    id = sid;

    // increase serial number of job
    job.SerialID++;

    return(true);
}

//------------------------------------------------------------------------------

bool CLauncher::IsJobFinished(CLauncherJob& job)
{
    // assemble command
    CFileName cmd;
    CFileSystem::GetCurrentDir(cmd);
    cmd = cmd / JobStatusWrapper
        + " " + job.JobID;

    stringstream str;

    FILE* p_file = popen(cmd,"r");
    if( p_file == NULL ){
        lout << ">>> WRAPPER ERROR: " << endl;
        lout << "popen failed, errno(" << strerror(errno) <<  ")" << endl;
        return( false );
    }

    // copy output to stream
    char line[100];
    while( fgets(line, 99, p_file) ){
        str << line;
    }
    int retcode = pclose(p_file);

    if( retcode != 0 ){
        lout << endl;
        lout << ">>> WRAPPER ERROR OUTPUT: " << endl;
        lout << str.str(); // show output only if an error
        lout << "Return code: " << retcode;
        if( retcode < 0 ){
            lout << " errno(" << strerror(errno) <<  ")";
        }
        lout << endl;
        return(false);
    }

    // extract job id
    stringstream idstr(str.str());
    char status;
    idstr >> status;

    if( status == 'C' ) return(true);   // finished
    if( status == 'U' ) return(true);   // unknown

    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
