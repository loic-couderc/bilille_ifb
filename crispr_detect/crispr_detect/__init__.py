from flask import Flask
from werkzeug import SharedDataMiddleware
#from flask_reverse_proxy import FlaskReverseProxied

#from flask_bootstrap import Bootstrap

app = Flask(__name__)
app.config.from_object('config')
#crispr_detect.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
#Bootstrap(crispr_detect)
#proxied = FlaskReverseProxied(app)

app.wsgi_app = SharedDataMiddleware(app.wsgi_app,
                                    {'/display': app.config['UPLOAD_FOLDER']})
#print(app.app_context().url_prefix())

from crispr_detect import views
