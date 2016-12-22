WTF_CSRF_ENABLED = True
SECRET_KEY = "you-will-never-guess-this-awesome-secret-key" #The WTF_CSRF_ENABLED setting activates the cross-site request forgery prevention
UPLOAD_FOLDER = '/crispr_detect/upload'
DEBUG = {{ microbannot.debug }}
THREADED = True
CRISPR_OUTPUT_MOUNT={{ microbannot.crispr_output_mount }} #is the dir used for --volume /tmp/crispr:/crispr_detect/upload:rw \
WEBSMASH_OUTPUT_MOUNT={{ microbannot.websmash_output_mount }} #is the dir used for --volume /tmp/websmash:/websmash/upload:rw \
