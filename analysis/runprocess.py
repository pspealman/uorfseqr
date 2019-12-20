#Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman
import simplejson as json
SETTINGS = json.load(open('analysis.json','r'))

import os
if not os.path.exists(SETTINGS['processed_dir']):
	os.mkdir(SETTINGS['processed_dir'])

import subprocess
for z in [x['RPF'] for x in SETTINGS['SAMPLES']]:
	subprocess.Popen(["python", "../analysis/newprocess.py", z])

