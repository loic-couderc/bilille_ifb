#!/bin/bash
set -e #Exit immediately if a command exits with a non-zero status.
set -u #Treat unset variables as an error when substituting.
set -x #Print commands and their arguments as they are executed.

if [ ! -d "{{ bilille.databases_dir }}" ]; then

	## 1_ set up the pfam database
	mkdir -p "{{ bilille.databases_dir }}/pfam" && cd "{{ bilille.databases_dir }}/pfam"
	
	curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz > Pfam-A.hmm.gz
	gunzip Pfam-A.hmm.gz
	apt-get install hmmer
	hmmpress Pfam-A.hmm #prepare an HMM database for hmmscan
	
	## 2_ set up the ClusterBlast database
	mkdir -p "{{ bilille.databases_dir }}/clusterblast" && cd "{{ bilille.databases_dir }}/clusterblast"
	wget https://bitbucket.org/antismash/antismash/downloads/clusterblast_dmnd07.tar.xz
	tar -xJf clusterblast_dmnd07.tar.xz
	rm clusterblast_dmnd07.tar.xz
fi
	
