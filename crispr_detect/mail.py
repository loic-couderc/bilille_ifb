from flask import Flask
from flask_mail import Mail, Message

app = Flask(__name__)
mail = Mail(app)

app.config['MAIL_DEBUG'] = True
app.config['MAIL_SERVER'] = 'smtps.univ-lille1.fr'
app.config['MAIL_PORT'] = 587
app.config['MAIL_USERNAME'] = 'lcouderc'
app.config['MAIL_PASSWORD'] = '********'
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USE_SSL'] = False
mail = Mail(app)


@app.route('/')
def index():
    msg = Message('Hello', sender='loic.couderc@univ-lille1.fr',
                  recipients=['loic.couderc@univ-lille1.fr'])
    msg.body = 'Hello Flask message sent from Flask-Mail'
    mail.send(msg)
    return 'Sent'

if __name__ == '__main__':
    app.run(debug=True, port=5005)
