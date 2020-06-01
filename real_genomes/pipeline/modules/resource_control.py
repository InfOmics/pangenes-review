#!/usr/bin/env python3
import os
import signal
from subprocess import Popen, PIPE, TimeoutExpired, call, DEVNULL, STDOUT
import re
	
def call_program(parameters, software_name):
	
	max_execution_time = 7200 #timeout is in seconds , 2h = 7200s
	
	cmd_list = ['bash','modules/runners.sh', software_name] + parameters
	#print(cmd_list)

	with Popen(cmd_list, shell=False, preexec_fn=os.setsid, stderr=PIPE) as process: #, stdout=PIPE
		try:
			
			err = process.communicate(timeout=max_execution_time)[1] #out, 
			stat = re.findall(r'.*#panreview#.*',err.decode("utf-8"))
			if software_name == 'panx':
				print(err)
			print(stat)
			return stat

		except TimeoutExpired:
			print('Timout (>2h) for '+software_name)
			os.killpg(process.pid, signal.SIGINT)
			
			return None

	

