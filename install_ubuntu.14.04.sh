#!/bin/bash
set -e #Exit immediately if a command exits with a non-zero status.
set -u #Treat unset variables as an error when substituting.
set -x #Print commands and their arguments as they are executed.

bilille_install_dir=$(pwd)

if [ ! -d "$bilille_install_dir" ]; then
 echo 'You have to provide path to the current install dir'
 exit 0;
fi

# configure
pip install j2cli
confs=( "${bilille_install_dir}/antismash/config/settings.py.tpl"
        "${bilille_install_dir}/crispr_detect/config/settings.py.tpl"
        "${bilille_install_dir}/daemon/docker-bilille.conf.tpl"
        "${bilille_install_dir}/docker-compose.yml.tpl"
        "${bilille_install_dir}/download_antismash_databases.sh.tpl"
        "${bilille_install_dir}/download_antismash_example_data.sh.tpl"
        "${bilille_install_dir}/create_output_dir.sh.tpl"
)

for cfg in ${confs[@]}; do
	j2 ${cfg} "${bilille_install_dir}/config.ini" > "${cfg%.*}"
done


# install

## 1_ setup databases
sh "${bilille_install_dir}/download_antismash_databases.sh"

## 2_ setup outpudir
sh "${bilille_install_dir}/create_output_dir.sh"

## 3_ setup antismash example
sh "${bilille_install_dir}/download_antismash_example_data.sh"

## 4_ setup the daemon
ln -f -s "${bilille_install_dir}/daemon/docker-bilille.conf" "/etc/init/docker-bilille.conf"
initctl reload-configuration
#systemctl enable

## 5_ reboot and check the application
#docker-compose up --build

#rsync -av . root@192.54.201.34:/root/mydisk/bilille

#docker stop $(docker ps -a -q) && docker rm $(docker ps -a -q) && docker-compose up



initctl start docker-bilille



#sudo apt-get install python3-pip
#pip3 install virtualenv
#virtualenv flask
#source flask/bin/activate
#pip3 install -r requirements.txt


