#!/usr/bin/env python

import common
import os
import subprocess
import sys
import time


resultDir = sys.argv[1]
vid       = sys.argv[2]
cid       = sys.argv[3]
fid       = sys.argv[4]
pluginDir = sys.argv[5]

os.chdir(resultDir)

"""
TODO: Integrate each script as a module in the plugin instead of running them as command. 
"""
common.printtime('Result Dir: '+resultDir)
common.printtime('Variant Caller ID: '+vid)
common.printtime('Coverage Analysis ID: '+cid)
common.printtime('fDiagnostics ID: '+fid)
common.printtime('Plugin Directory: '+pluginDir)

reportFile = 'AgriSumToolkit_%s_%s.report' % (vid,cid)
common.printtime('Report File: '+reportFile)

# CalculateRunSummary.R
command = "Rscript %s/scripts/CalculateRunSummary.R %s %s %s %s" % ( pluginDir, resultDir, vid, cid, fid )
common.printtime('Running: '+command)
subprocess.call(command, shell=True)
common.printtime('Finished: '+command)

# AgriSeqGBSmatrix.py
command = "python %s/scripts/AgriSeqGBSmatrix.py %s %s" % ( pluginDir, resultDir, vid )
common.printtime('Running: '+command)
subprocess.call(command, shell=True)
common.printtime('Finished: '+command)

#TOP BOTTOM
command = "python %s/scripts/Convert2TopBot.py %s %s" % ( pluginDir, resultDir, vid )
common.printtime('Running: '+command)
subprocess.call(command, shell=True)
common.printtime('Finished: '+command)

#TOP BOTTOM
command = "python %s/scripts/Genotypes2TopBottom.py %s %s" % ( pluginDir, resultDir, vid )
common.printtime('Running: '+command)
subprocess.call(command, shell=True)
common.printtime('Finished: '+command)

f=open(reportFile,"w+")
f.write('1')
f.close()


