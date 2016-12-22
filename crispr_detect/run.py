#!flask/bin/python
from crispr_detect import app

if __name__ == '__main__':
    app.run(threaded=True)