# Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman
import os
import simplejson as json
SETTINGS = json.load(open('analysis.json','r'))

import subprocess
for x in SETTINGS['SAMPLES']:
	for ch in SETTINGS['chromosomes']:
		for kindof in ['']+['p'+str(i)+'-' for i in range(len(SETTINGS['permutation_seeds']))]:
			fout = None #open(x['name']+ch+'.stdout','w')
			fout2 = None #open(x['name']+ch+'.stderr','w')
			if not os.path.exists('predictions.'+kindof+x['name']):
				os.mkdir('predictions.'+kindof+x['name'])
			subprocess.Popen(["python", "../analysis/newattempt.py", x['name'], ch, kindof], stdout=fout, stderr=fout2)
		# fout.close()
		# fout2.close()

