import os
from pathlib import Path
import argparse
import datetime
import time
import sys
import logging

def check_run_completion(run_dir):
	flag_file = Path(run_dir)
	if flag_file.is_file():
		return True
	else:
		return False

def read_run_params():

    import yaml

    with open('/storage/scratch/covid/container/bcl2fastq/config.yaml') as f:
        params = yaml.safe_load(f)

    return params


def bashCommunicator(command,output_expected=False):

    import subprocess

    process = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Process failed %s %d %s %s" % (command,process.returncode, stdout, stderr))
    else:
        if output_expected:
            return [x for x in stdout.split("\n")]

def extract_job_exit_status(job_id):
	try:
		exit_info = bashCommunicator("/bio/gridengine/bin/lx-amd64/qacct -j "+job_id ,output_expected=True)
		exit_status = ""
		for st in exit_info:
			if "exit_status" in st:
				exit_status = st
		return exit_status.replace(" ","").replace("exit_status","")

	except Exception as e:
		return "ERROR"

if __name__ == '__main__':

	params = read_run_params()

	script = params["container"] + params["script_input"]

	run_dirs = params["run_folders"]

	LOG_FILENAME = params["container"]+"bcl2fastq/log/log_rclone_"+str(datetime.datetime.now()).replace(" ","_")+".log"

	logging.basicConfig(filename=LOG_FILENAME,
						format='%(asctime)s %(levelname)-8s %(message)s',
    					level=logging.INFO,
    					datefmt='%Y-%m-%d %H:%M:%S')

	os.environ['SGE_CELL']="default"
	os.environ['SGE_ARCH']="lx-amd64"
	os.environ['SGE_ROOT']="/bio/gridengine"
	os.environ['SGE_CLUSTER_NAME']="default"

	for run in run_dirs:

		run_dir = run_dirs[run]
		print("checking rclone for ... "+run_dir)
		logging.info("checking rclone for ... "+run_dir)

		if check_run_completion(run_dir+"bcl2fastq_processing.txt") and \
		   check_run_completion(run_dir+"bcl2fastq_job_id.txt") and \
		   not check_run_completion(run_dir+"rclone_processing.txt"):

			logging.info("file exists--"+run_dir+"bcl2fastq_processing.txt")
			logging.info("file exists--"+run_dir+"bcl2fastq_job_id.txt")
			logging.info("file does not exist--"+run_dir+"rclone_processing.txt")

			job_id = ""
			with open(run_dir+"bcl2fastq_job_id.txt") as f:
				line = f.readline()
				job_id = line.split(':')[1]
				job_id = job_id.rstrip()

			logging.info("job_id is... "+ job_id)

			job_exit_status = extract_job_exit_status(job_id)

			logging.info("job exit status is... "+ job_exit_status)

			if job_exit_status == "0":

				cmd = "/bin/touch  "+ run_dir + "rclone_processing.txt"

				bashCommunicator(cmd,output_expected=False)

				cmd = "/usr/bin/nohup /storage/apps/opt/rclone/rclone-v1.55.0-linux-amd64/rclone copy "+ run_dir+"out2/  box:Houston-SARSCov2/"+run+"-fastq/ --include 'VCoV-*.fastq.gz'  & "

				logging.info(cmd)

				cmd_out = bashCommunicator(cmd,output_expected=True)

				logging.info("process id is... "+str(cmd_out))

			else:
				logging.info("rlcone not submitted for... "+job_id)
				logging.info("waiting for bcl2fastq to complete... "+ run)
				logging.info("rclone did not process for... "+run_dir)

		else:
			if check_run_completion(run_dir+"rclone_processing.txt"):
				logging.info("rclone is already processing or completed... "+ run)
			else:
				logging.info("run is not completed...rclone will run after bcl2fastq is completed... "+ run)
				logging.info("rclone did not process for... "+run_dir)
