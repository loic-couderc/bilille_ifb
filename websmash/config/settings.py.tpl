############# Configuration #############
DEBUG = {{ bilille.debug }}
SECRET_KEY = 'replace this by a long random awesome string'
RESULTS_PATH = "/upload"
RESULTS_URL = "/upload"

# Flask-Mail settings
MAIL_SERVER = {{ bilille.mail_server }}
MAIL_PORT = {{ bilille.mail_port }}
MAIL_USERNAME = {{ bilille.mail_username }}
MAIL_PASSWORD = {{ bilille.mail_password }}
MAIL_USE_TLS = True
MAIL_USE_SSL = False
DEFAULT_MAIL_SENDER = {{ bilille.default_mail_sender }}
DEFAULT_RECIPIENTS = {{ bilille.default_recipients }}


# Flask-Redis settings
REDIS_URL = "redis://localhost/0"
#########################################
