############# Configuration #############
DEBUG = {{ microbannot.debug }}
SECRET_KEY = 'replace this by a long random awesome string'
RESULTS_PATH = "/websmash/upload"
RESULTS_URL = "/upload"

# Flask-Mail settings
MAIL_SERVER = "mail"
DEFAULT_MAIL_SENDER = "antismash@example.org"
DEFAULT_RECIPIENTS = ["antismash@example.org"]

# Flask-Redis settings
REDIS_URL = "redis://localhost/0"
#########################################
