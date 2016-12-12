#Building of docker image

0. Make sur to install databases as mentionned here: http://docs.antismash.secondarymetabolites.org/install/#antismash-standalone-lite

1. Build docker image:
	```bash
	docker build -t bilille .
	```

2. run the container (update path accordingly to your environment)
	```bash
	docker run --rm --name bilille \
	--publish 80:80 \
	--volume /home/lcouderc/workspace/bilille_docker_only:/config:ro \
	--volume /data/databases:/databases:ro \
	--volume /tmp/crispr:/crispr_detect/upload:rw \
	--volume /tmp/websmash:/websmash/upload:rw \
	bilille
	```

3. test
	```bash
	docker exec -i -t bilille /bin/bash -c "cd /crispr_detect && flask/bin/python -m unittest discover tests"
	docker exec -i -t bilille /bin/bash -c "cd /websmash && pip install -r test_requirements.txt && python -m unittest discover tests"
	docker exec -i -t bilille /bin/bash -c "run_crispr_detect -h"
	docker exec -i -t bilille /bin/bash -c "run_antismash -h"
	curl -I localhost
	curl -I http://localhost/upload/example/index.html
	```

4. Use your favorite navigator and connect to localhost