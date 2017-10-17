#!/bin/bash
set -e #Exit immediately if a pipeline returns a non-zero status

#0_ configure
echo "configure"
confs=(
	"/websmash/config/settings.py.tpl"
	"/webannot/config/settings.py.tpl"
)

for cfg in ${confs[@]}; do
	j2 ${cfg} "/config/config.ini" > "${cfg%.*}"
done

#1_ run redis_server
echo "run redis server"
redis-server&

#2_ run runsmash
echo "run the antiSMASH job runner(s)"
/runsmash/runSMASH --queue "redis://localhost:6379/0" --workdir /websmash/upload --statusdir /websmash/upload/status --name runner01 --long-running --cpus $(grep -c ^processor /proc/cpuinfo)&

#3_ run the antiSMASH job status monitor
echo "the antiSMASH job runner(s)"
ln -s /websmash/example/ /websmash/upload/example
/runsmash/watchStatus --queue "redis://localhost:6379/0" --statusdir /websmash/upload/status&

#4_ run nginx
echo "run nginx"
nginx -g "daemon off;"&

#5_ run websmash
echo "run the antiSMASH web interface"
cd /websmash && /usr/local/bin/gunicorn --env WEBSMASH_CONFIG=/websmash/config/settings.py -b 0.0.0.0:8000 websmash:app&

#6_ run crispr
echo "run the crispr web interface"
cd /webannot && /webannot/flask/bin/gunicorn --env WEBANNOT_CONFIG=/webannot/config/settings.py -b 0.0.0.0:8001 webannot:app&

#supervisord -n
while true; do sleep 1000; done #to keep docker container alive
