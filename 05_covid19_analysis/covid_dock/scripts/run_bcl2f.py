import os
from pathlib import Path
import argparse
import datetime
import time
import sys
import logging

def read_run_params():

	import yaml

	with open('/storage/scratch/covid/container/bcl2fastq/config.yaml') as f:
		params = yaml.safe_load(f)
	
	return params


def bashCommunicator(command,logging,output_expected=False):

	import subprocess

	process = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
	stdout, stderr = process.communicate()
	
	if process.returncode != 0:
		print("Process failed %s %d %s %s" % (command,process.returncode, stdout, stderr))
	else:
		if output_expected:
			return [x for x in stdout.split("\n")]


def check_run_completion(run_dir):
	flag_file = Path(run_dir)
	if flag_file.is_file():
		return True
	else:
		return False


def log_msg(msg):
	tt = datetime.datetime.now()
	return ("[%s] %s\n"%(tt.strftime("%Y-%m-%d %H:%M:%S"), msg))


def extract_job_status(job_info):
	from bs4 import BeautifulSoup
	soup = BeautifulSoup(job_info[0], 'html.parser')
	return soup.job_state.contents[0]



if __name__ == '__main__':


	params = read_run_params()

	script = params["script_input"]

	run_dirs = params["run_folders"]

	LOG_FILENAME = params["container"]+"bcl2fastq/log/log_"+str(datetime.datetime.now()).replace(" ","_")+".log"

	logging.basicConfig(filename=LOG_FILENAME,
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')

	logging.info(os.environ)

	os.environ['SGE_CELL']="default"
	os.environ['SGE_ARCH']="lx-amd64"
	os.environ['SGE_ROOT']="/bio/gridengine"
	os.environ['SGE_CLUSTER_NAME']="default"

	logging.info("Added SGE variables..")

	logging.info(os.environ)

	for run in run_dirs:

		run_dir = run_dirs[run]

		logging.info("processing ..."+run_dir)

		if check_run_completion(run_dir+"CopyComplete.txt") and \
		   check_run_completion(run_dir+"SampleSheet.csv") and \
		   not check_run_completion(run_dir+"bcl2fastq_processing.txt") :

			logging.info("file exists--"+run_dir+"CopyComplete.txt")
			logging.info("file exists--"+run_dir+"SampleSheet.txt")
			logging.info("file does NOT exist--"+run_dir+"bcl2fastq_processing.txt")

			cmd = "/bin/touch  "+ run_dir + "bcl2fastq_processing.txt"

			bashCommunicator(cmd,logging,output_expected=False) 
			
			cmd = "/bio/gridengine/bin/lx-amd64/qsub -V -m abe -M ssubedi@houstonmethodist.org /storage/apps/pipelines/bcl2fastq/runBcl2fastq_qsubResearch.sh  -d " + run_dir + " -i novaseq -b 1 "
			 
			logging.info(cmd)
			
			cmd_ret = bashCommunicator(cmd,logging,output_expected=True) 
			logging.info(cmd_ret)

		
			job_id = ""
			for message in cmd_ret:
				if "runBcl2fastq_qsubResearch.sh" in message:
					job_id = message.split( )[2]
					break

			cmd = "/usr/bin/echo  job_id:"+ job_id + " > " + run_dir + "bcl2fastq_job_id.txt"

			bashCommunicator(cmd,logging,output_expected=False) 

		else:
			if not check_run_completion(run_dir+"CopyComplete.txt"):
				logging.info("run is not completed!")
			elif check_run_completion(run_dir+"bcl2fastq_processing.txt"):
				logging.info("bcl2fastq is currently running or completed running!")

