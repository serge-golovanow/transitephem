# -*- coding: utf-8 -*-
"""
Created on Wed May 14 13:51:43 2014

@author: bmorris
"""

#import os
#import webbrowser
import argparse
from calculateEphemerides import calculateEphemerides
#--latitude=45.7778123 --longitude=4.9226306 --elevation=170 --horizon=30 --start=2018,2,1 --end=2018,3,1 --mag=13 --depth=0.01 --twilight=-12 --duration=2

parser = argparse.ArgumentParser(description='Previsions de transits exoplanetes')
parser.add_argument('--latitude', action="store", dest="latitude")
parser.add_argument('--longitude', action="store", dest="longitude")
parser.add_argument('--elevation', action="store", dest="elevation")
parser.add_argument('--horizon', action="store", dest="horizon")
parser.add_argument('--start', action="store", dest="start")
parser.add_argument('--end', action="store", dest="end")
parser.add_argument('--mag', action="store", dest="mag")
parser.add_argument('--depth', action="store", dest="depth")
parser.add_argument('--twilight', action="store", dest="twilight")
parser.add_argument('--duration', action="store", dest="duration")
conf = vars(parser.parse_args())


# Path to ".par" file, with the observatory parameters
parfile = 'mro.par'

# Path to the output HTML file -- this script will assume that the
# "html_out" parameter in the .par file is set to True
#outputPath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'outputs','eventReport.html')

# Calculate the ephemeris for all visible planets within stated limits
calculateEphemerides(parfile,conf)

# Open the HTML output in a new tab of the default browser
# webbrowser.open_new_tab("file:"+2*os.sep+outputPath)