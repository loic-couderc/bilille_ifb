############# Configuration #############
DEBUG = False
SECRET_KEY = "replace this by a long random string"
RESULTS_PATH = "/upload"
RESULTS_URL = "/upload"

# Flask-Mail settings
MAIL_SERVER = 'smtps.univ-lille1.fr'
MAIL_PORT = 587
MAIL_USERNAME = 'lcouderc'
MAIL_PASSWORD = '********'
MAIL_USE_TLS = True
MAIL_USE_SSL = False
DEFAULT_MAIL_SENDER = 'loic.couderc@univ-lille1.fr'
DEFAULT_RECIPIENTS = ['loic.couderc@univ-lille1.fr']

# Flask-Redis settings
REDIS_URL = "redis://redis/0"
#########################################'
