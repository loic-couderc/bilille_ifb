#!/bin/bash
set -e #Exit immediately if a command exits with a non-zero status.
set -u #Treat unset variables as an error when substituting.
set -x #Print commands and their arguments as they are executed.

if [ ! -d "{{ bilille.install_dir }}/antismash/example" ]; then
	cd "{{ bilille.install_dir }}/antismash/"
	wget http://antismash.secondarymetabolites.org/upload/example/NC_003888.3.zip && unzip NC_003888.3.zip
	mv 2acd7e9e-4872-48d4-bae9-cac30ec52622 example
	rm NC_003888.3.zip
fi
