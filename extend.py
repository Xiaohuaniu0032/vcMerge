#!/usr/bin/python
import os
import subprocess
import pprint
import sys
import time
import requests

 
def fileExist(bucket):
    """
    Simple function that checks
    whether given path on the server is present
    Input: File Path
    Return: bool
    """
    path = bucket["request_get"].get("path", '')
    return(os.path.isfile(path))
    

def reportAvailable(bucket):
    """
    Checks the status or initializes the run if it is not queued yet.
    If the run is still in progress, it returns the status.
    If it has not run before, it queues using the drmaa see: def start
    """
    report = bucket["request_get"].get("report", '')
    resultDir = bucket["request_get"].get("results_dir", '')
    vid = bucket["request_get"].get("vid", '')
    cid = bucket["request_get"].get("cid", '')
    fid = bucket["request_get"].get("fid", '')
    pluginBase = "/results/plugins"
    pluginName = bucket["name"]
    
    reportFile = "%s/%s.report" % (resultDir, report)
    
    out = {}
    out['status'] = "NA";
    
    if os.path.isfile(reportFile):
        # return 1
        with open(reportFile) as f: s = f.read()
        if s == '1':
            out['status'] = 'Completed'
        else:
            out['status'] = checkStatus(s)
        return out
    else:
        
        try:
            pluginDir = "%s/%s" % (pluginBase, pluginName)
            script = "%s/pipeline.py" % (pluginDir)
            jobid = start(resultDir,"python", [script,resultDir,vid,cid,fid,pluginDir])
            
            if int(jobid) > 0 :
                f=open(reportFile,"w+")
                f.write(jobid)
                f.close()
                out['jobid'] = jobid
                out['status'] = 'job is running'
            else:
                out['error'] = 'DRMAA FAILED'
        
            return out
        except Exception, e:
            out['error'] = str(e)
            return out
            
    return out

def createDrmaaSession():
    """
    Common environment set for drmaa and the session initialization
    """
    import drmaa

    os.environ["SGE_ROOT"] = "/var/lib/gridengine"
    os.environ["SGE_CELL"] = "iontorrent"
    os.environ["SGE_CLUSTER_NAME"] = "p6444"
    os.environ["SGE_QMASTER_PORT"] = "6444"
    os.environ["SGE_EXECD_PORT"] = "6445"
    os.environ["SGE_ENABLED"] = "True"
    os.environ["DRMAA_LIBRARY_PATH"] = "/usr/lib/libdrmaa.so.1.0"

    _session = drmaa.Session()
    try:
        _session.initialize()
    except:
        _session.exit()
        #maybe need to return False?
        return False
    return _session

def start(path, script, args):
    """
    Queue the job using drmaa
    Input:
    Return: jobid
    """
    _session = createDrmaaSession()
    # Prepare drmaa Job - SGE/gridengine only
    jt = _session.createJobTemplate()
    jt.nativeSpecification = " -q %s" % ("plugin.q")

    jt.workingDirectory = path
    jt.outputPath = ":" + os.path.join(path, "AST_drmaa_stdout.txt")
    jt.joinFiles = True # Merge stdout and stderr

    jt.remoteCommand = script
    jt.args = args    

    # Submit the job to drmaa
    jobid = _session.runJob(jt)

    _session.deleteJobTemplate(jt)
    _session.exit()
    return jobid

def checkStatus(jobid):
    """
    Checks the status of the job and returns the result in human readable format.
    It also handles some rare drmaa related exceptions as well.
    """
    import drmaa
    
    _session = createDrmaaSession()
    
    decodestatus = {drmaa.JobState.UNDETERMINED: 'process status cannot be determined',
                        drmaa.JobState.QUEUED_ACTIVE: 'job is queued and active',
                        drmaa.JobState.SYSTEM_ON_HOLD: 'job is queued and in system hold',
                        drmaa.JobState.USER_ON_HOLD: 'job is queued and in user hold',
                        drmaa.JobState.USER_SYSTEM_ON_HOLD: 'job is queued and in user and system hold',
                        drmaa.JobState.RUNNING: 'job is running',
                        drmaa.JobState.SYSTEM_SUSPENDED: 'job is system suspended',
                        drmaa.JobState.USER_SUSPENDED: 'job is user suspended',
                        drmaa.JobState.DONE: 'job finished normally',
                        drmaa.JobState.FAILED: 'job finished, but failed'}
    
    status = 0
    
    try:
        jobid = str(jobid)
        status = decodestatus[_session.jobStatus(jobid)]
    except Exception as err:
        if "code 18" in str(err) : #this is an indication of job no longer is running (most probably done)
            return -1
    
    _session.exit()
    return status

# if __name__ == '__main__':
#     print checkStatus(96019)