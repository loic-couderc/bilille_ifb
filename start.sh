#!/bin/bash

#0_ configure
echo "configure"
confs=(
	"/antismash/config/settings.py.tpl"
	"/crispr_detect/config/settings.py.tpl"
)

for cfg in ${confs[@]}; do
	j2 ${cfg} "/config/config.ini" > "${cfg%.*}"
done

#1_ run redis_server
echo "run redis server"
redis-server&

#2_ run runsmash
echo "run the antiSMASH job runner(s)"
/runsmash/runSMASH --queue "redis://localhost:6379/0" --workdir /upload --statusdir /upload/status --name runner01 --long-running&

#3_ run nginx
echo "run nginx"
nginx -g "daemon off;"&

#4_ run websmash
echo "run the antiSMASH web interface"
cd /websmash && /usr/local/bin/gunicorn --env WEBSMASH_CONFIG=/websmash/config/settings.py -b 0.0.0.0:8000 websmash:app&

#5_ run crispr
echo "run the crispr web interface"
cd /crispr_detect && /crispr_detect/crispr_detect/flask/bin/gunicorn --env CRISPR_DETECT_CONFIG=/crispr_detect/crispr_detect/config/settings.py -b 0.0.0.0:8001 crispr_detect:app&

#supervisord -n
while true; do sleep 1000; done #to keep docker container alive
