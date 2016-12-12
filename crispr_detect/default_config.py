WTF_CSRF_ENABLED = True
# The WTF_CSRF_ENABLED setting activates the cross-site request forgery
# prevention
SECRET_KEY = "you-will-never-guess-my-cristal-secret-key"
UPLOAD_FOLDER = '/tmp/crispr'
DEBUG = True
THREADED = True
CRISPR_OUPUT_MOUNT= "/tmp/crispr" #is the dir used for --volume /tmp/crispr:/crispr_detect/upload:rw \
WEBSMASH_OUTPUT_MOUNT= "/tmp/websmash" #is the dir used for --volume /tmp/websmash:/websmash/upload:rw \