# CRISPR web server component
# VERSION 0.1.0
FROM debian:jessie
MAINTAINER Loïc Couderc <loic.couderc@univ-lille1.fr>

RUN apt-get update && apt-get install -y \
    python3-pip

ADD . /crispr_detect
WORKDIR  /crispr_detect
RUN pip3 install -r requirements.txt && \
    pip3 install gunicorn

VOLUME ['/upload', '/config']

ENTRYPOINT ["/usr/local/bin/gunicorn"]
CMD ["--env", "CRISPR_DETECT_CONFIG=/config/settings.py",  "-b", "0.0.0.0:8000", "crispr_detect:app" ]
