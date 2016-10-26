from flask import Flask
#from flask_bootstrap import Bootstrap

app = Flask(__name__)
app.config.from_object('config')
#crispr_detect.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
#Bootstrap(crispr_detect)

from crispr_detect import views
