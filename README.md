#Building of docker image

0. Make sur to install databases as mentionned here: http://docs.antismash.secondarymetabolites.org/install/#antismash-standalone-lite

1. Build docker image:
	```bash
	docker build -t microbannot .
	```
	you can save some space by "squashing" the image (https://github.com/goldmann/docker-squash/):
	```bash
	pip install https://github.com/goldmann/docker-squash/archive/master.zip
	docker-squash microbannot -t microbannot
	```

2. run the container (update path accordingly to your environment)
	```bash
	docker run --rm --name microbannot \
	--publish 80:80 \
	--volume /home/lcouderc/workspace/microbannot:/config:ro \
	--volume /data/databases:/databases:ro \
	--volume /tmp/webannot:/webannot/upload:rw \
	--volume /tmp/websmash:/websmash/upload:rw \
	microbannot
	```

3. test
	```bash
	docker exec -i -t microbannot /bin/bash -c "cd /webannot && flask/bin/python -m unittest discover tests"
	docker exec -i -t microbannot /bin/bash -c "cd /websmash && pip install -r test_requirements.txt && python -m unittest discover tests"
	docker exec -i -t microbannot /bin/bash -c "run_crispr_detect -h"
	docker exec -i -t microbannot /bin/bash -c "run_antismash -h"
	curl -I localhost
	curl -I http://localhost/upload/example/index.html
	```

4. Use your favorite navigator and connect to localhost

5. Optional (push image)
	```bash
	docker tag <id> docker-registry.genouest.org/bilille/microbannot
	docker login docker-registry.genouest.org
	docker push docker-registry.genouest.org/bilille/microbannot
	```
The following images are identical:
* docker pull lcouderc/microbannot
* docker pull docker-registry.genouest.org/bilille/microbannot 
