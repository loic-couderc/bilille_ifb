from flask import Flask
from werkzeug import SharedDataMiddleware
#from flask_reverse_proxy import FlaskReverseProxied

#from flask_bootstrap import Bootstrap

app = Flask(__name__)
app.config.from_object('default_config')
app.config.from_envvar('WEBANNOT_CONFIG', silent=True)

#proxied = FlaskReverseProxied(app)

app.wsgi_app = SharedDataMiddleware(app.wsgi_app,
                                    {'/display': app.config['UPLOAD_FOLDER']})
# print(app.app_context().url_prefix())

from webannot import views