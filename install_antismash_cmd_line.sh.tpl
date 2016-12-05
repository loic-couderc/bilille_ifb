#!/bin/bash
set -e #Exit immediately if a command exits with a non-zero status.
set -u #Treat unset variables as an error when substituting.
set -x #Print commands and their arguments as they are executed.

curl -q https://bitbucket.org/antismash/docker/raw/HEAD/standalone-lite/run_antismash > /usr/local/bin/run_antismash
chmod a+x /usr/local/bin/run_antismash
sed -i 's|readonly DATABASE_DIR="/data/databases"|readonly DATABASE_DIR="{{ bilille.databases_dir }}"|' /usr/local/bin/run_antismash