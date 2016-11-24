#!/bin/bash
set -e #Exit immediately if a command exits with a non-zero status.
set -u #Treat unset variables as an error when substituting.
set -x #Print commands and their arguments as they are executed.

mkdir -p "{{ bilille.output_dir }}"