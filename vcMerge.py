#!/usr/bin/env python

# Copyright (C) 2019 Thermo Fisher Scientific. All Rights Reserved
# Author: Haktan Suren

# This has been tested only on S5XL, S5 Prime with 5.10 and 5.12. 

import common
import glob
import sys
import subprocess
import json
import os
import re
import shutil
from ion.plugin import *

#from django.conf import settings
#from django.template.loader import render_to_string




class vcMerge(IonPlugin):
  version = '1.0.0.0'
  envDict = dict(os.environ)
  author = "longfei.fu@thermofisher.com"
  runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]

  def merge(self):
    report_dir = self.envDict['ANALYSIS_DIR'] # /results/analysis/output/Home/Auto_xxx
    outdir = '/results/analysis' + self.envDict['TSP_URLPATH_PLUGIN_DIR'] # /results/analysis/output/Home/Auto_xxx/plugin_out/variantCallerMerge.x
    this_dir = self.envDict['DIRNAME']
    cmd = "perl %s/merge_main.pl %s %s" % (this_dir,report_dir,outdir)

    print "report dir is: %s" % (report_dir)
    print "plugin outdir is: %s" % (outdir)
    print "cmd is: %s" % (cmd)

    os.system(cmd)
    print "Finished the variantCallerMerge plugin."

  def launch(self,data=None):
    print "Start running the variantCallerMerge plugin."

    start_plugin_json = '/results/analysis' + os.path.join(os.getenv('TSP_URLPATH_PLUGIN_DIR'),'startplugin.json')
    print(start_plugin_json)
    with open(start_plugin_json, 'r') as fh:
      spj = json.load(fh)
      net_location = spj['runinfo']['net_location']

    # /results/analysis/output/Home/S5yanzheng-20211223-chip2-MeanAccuracy_v2_482
    # here report_number is 482
    report_dir = self.envDict['ANALYSIS_DIR'] # /results/analysis/output/Home/Auto_xxx
    report_number = report_dir.split('_')[-1]
    this_plugin_dir = os.getenv('TSP_URLPATH_PLUGIN_DIR').split('/')[-1]
    url_path = os.getenv('TSP_URLPATH_PLUGIN_DIR')
    # ../../../../../report/{{runinfo_pk}}/metal/plugin_out/{{output_dir}}
    
    '''
    full_link = "../../../../../report/%s/metal/plugin_out/%s" % (report_number,this_plugin_dir)
    #print(full_link)

    with open("status_block.html", "w") as html_fp:
      html_fp.write('<html><body><pre>')
      html_fp.write('<tr><td>Output directory:</td>       <a href="%s" target="_blank">%s</a></td></tr>' 
        % (full_link,this_plugin_dir))
      html_fp.write("</pre></body></html>")
    self.merge()
    '''
    self.merge()
    with open("variantCallerMerge_block.html","w") as f:
      f.write('<html><body>To download: "Right Click" -> "Save Link As..."<br>\n')
      for merge_vcf in glob.glob('*.merged.vcf.xls'):
        print(merge_vcf)
        f.write('<a href="%s">%s</a><br>\n'
                % (os.path.join(net_location,url_path,merge_vcf),merge_vcf))
        f.write('</body></html>\n')
    
    return True


if __name__ == "__main__":
    PluginCLI()
