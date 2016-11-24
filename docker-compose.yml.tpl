version: '2'
services:
  antismash-redis:
    restart: always
    image: redis

  antismash-run:
    restart: always
    volumes:     
      - "{{ bilille.output_dir }}/antismash:/upload:rw"
      - "{{ bilille.databases_dir }}:/databases:ro"
    links:
      - antismash-redis:redis
    image: antismash/runsmash-lite
    environment:
      - ANTISMASH_EMAIL_FROM={{ antismash.default_mail_sender }}
      - ANTISMASH_EMAIL_ERROR={{ antismash.default_mail_sender }}
      - ANTISMASH_EMAIL_HOST={{ antismash.mail_server }}
      - ANTISMASH_EMAIL_ENCRYPT=TLS
      - ANTISMASH_EMAIL_USER={{ antismash.mail_username }}
      - ANTISMASH_EMAIL_PASSWORD={{ antismash.mail_password }}
    command: /runsmash/runSMASH --queue {{ antismash.redis_url }} --workdir /upload --statusdir /upload/status --name runner01 --long-running

  antismash-job-monitor:
    restart: always
    volumes:
      - "{{ bilille.output_dir }}/antismash:/upload"
      - "{{ bilille.install_dir }}/antismash/config:/config"
    links:
      - antismash-redis:redis

    image: antismash/runsmash-lite
    command: /runsmash/watchStatus --queue {{ antismash.redis_url }} --statusdir /upload/status


  nginx-proxy:  
    image: nginx
    ports:
      - "80:80"
    volumes:
      - "{{ bilille.install_dir }}/nginx:/etc/nginx/"
    restart: always
    links:
      - site1
      - site2


  site1:
    restart: always
    volumes:
      - "{{ bilille.output_dir }}/antismash:/upload"
      - "{{ bilille.install_dir }}/antismash/config:/config"
    links:
      - antismash-redis:redis
    image: antismash/websmash
#    container_name: site1


  site2:
    restart: always
    volumes:
      - {{ bilille.output_dir }}/crispr_detect:/upload
      - {{ bilille.install_dir }}/crispr_detect/config:/config
    build:
      context: ./crispr_detect
#    container_name: site2

#  site2:
#    restart: always
#    image: tutum/hello-world
#    container_name: site2





