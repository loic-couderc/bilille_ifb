#!/bin/bash
set -e #Exit immediately if a command exits with a non-zero status.
set -u #Treat unset variables as an error when substituting.
set -x #Print commands and their arguments as they are executed.

curl -q https://raw.githubusercontent.com/lccouderc/CRISPRDetect/master/CRISPRFinder_beta2_11.py  > /usr/local/bin/run_crispr_detect.py
chmod a+x /usr/local/bin/run_crispr_detect.py