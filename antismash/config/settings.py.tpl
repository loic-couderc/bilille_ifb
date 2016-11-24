############# Configuration #############
DEBUG = {{ antismash.debug }}
SECRET_KEY = {{ antismash.secret_key }}
RESULTS_PATH = "/upload"
RESULTS_URL = "/upload"

# Flask-Mail settings
MAIL_SERVER = {{ antismash.mail_server }}
MAIL_PORT = {{ antismash.mail_port }}
MAIL_USERNAME = {{ antismash.mail_username }}
MAIL_PASSWORD = {{ antismash.mail_password }}
MAIL_USE_TLS = True
MAIL_USE_SSL = False
DEFAULT_MAIL_SENDER = {{ antismash.default_mail_sender }}
DEFAULT_RECIPIENTS = {{ antismash.default_recipients }}

# Flask-Redis settings
REDIS_URL = {{ antismash.redis_url }}
#########################################'
