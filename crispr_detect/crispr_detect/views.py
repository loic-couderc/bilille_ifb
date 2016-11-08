from crispr_detect import app
from flask import render_template, abort, flash, redirect, request, redirect, url_for
from jinja2 import TemplateNotFound
from .forms import CrisprFinderForm
#from marshmallow import Schema, fields
from flask_restful import reqparse
from werkzeug.utils import secure_filename
from webargs.flaskparser import parser
import werkzeug
import threading
import os, sys
import uuid
from .crispr_detect import FindCRISPRs

def parse_crispr_finder_form(request):
    parser = reqparse.RequestParser()
    
    parser.add_argument('sequence', type=werkzeug.datastructures.FileStorage, location='files')#, required=True)
    parser.add_argument('k_mer_size_filter', type=int)
    parser.add_argument('pattern', type=str)
    parser.add_argument('window_size', type=int)
    parser.add_argument('allowed_mismatch', type=int)
    parser.add_argument('spacer_dr_match_limit', type=int)
    parser.add_argument('min_dr', type=int)
    parser.add_argument('max_dr', type=int)
    parser.add_argument('min_spacer_dr_ratio', type=float)
    parser.add_argument('max_spacer_dr_ratio', type=float)
    parser.add_argument('first_pass_limit', type=int)
    parser.add_argument('search_tracrrna', type=bool)
    return parser.parse_args()

# @copy_current_request_context
def crispr_finder_runner(**args):
    try:

        processing = create_flag_file(args['outputpath'], 'processing')
        stdout_file = os.path.join(os.path.join(args['outputpath'], 'stdout'))
        stderr_file = os.path.join(os.path.join(args['outputpath'], 'stderr'))
        #global old_stdout, old_stderr
        with open(stdout_file, 'w') as current_stdout, open(stderr_file, 'w') as current_stderr:
            old_stdout, old_stderr = sys.stdout, sys.stderr
            sys.stdout, sys.stderr = current_stdout, current_stdout
            findCRISPRs = FindCRISPRs(args['inputpath'],
                                      args['outputpath'],
                                      args['k_mer_size_filter'],
                                      args['pattern'],
                                      args['window_size'],
                                      args['allowed_mismatch'],
                                      args['spacer_dr_match_limit'],
                                      args['min_dr'],
                                      args['max_dr'],
                                      args['min_spacer_dr_ratio'],
                                      args['max_spacer_dr_ratio'],
                                      args['first_pass_limit'],
                                      args['search_tracrrna'])
            findCRISPRs.analyze()
            #create_flag_file(args['outputpath'], 'terminated')
    except:
        pass
        #create_flag_file(args['outputpath'], 'terminated')
        #create_flag_file(args['outputpath'], 'aborted')
    finally:
        os.unlink(processing)
        sys.stdout, sys.stderr = old_stdout, old_stderr

def create_flag_file(basedir, name):
    path = os.path.join(basedir, name)
    with open(path, 'w') as flag:
        return path


@app.route('/')
@app.route('/index/')
def index():
    user = "tpt" #TODO : clean this
    return render_template('index.html',
                            title='Home',
                            user=user)


@app.route('/crispr_finder/', methods=['GET', 'POST'])
# @use_args(CmdSchema())
def crispr_finder():
    if request.method == 'POST':
        #form = CrisprFinderForm(request.form) #http://stackoverflow.com/questions/30175285/flask-wtf-upload-file-error
        form = CrisprFinderForm()
        if form.validate():
            job_id = str(uuid.uuid4())

            results_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
            os.mkdir(results_dir)

            args = parse_crispr_finder_form(request)
            sequence_file = args['sequence']
            sequence_filename = secure_filename(sequence_file.filename)
            sequence_filepath = os.path.join(results_dir, sequence_filename)
            sequence_file.save(sequence_filepath)

            args['inputpath'] = sequence_filepath
            args['outputpath'] = results_dir
            args['job_id'] = job_id
            run_crispr = threading.Thread(target=crispr_finder_runner, kwargs=args)
            run_crispr.start()
            #print(url_for("crispr_finder_result", uuid=job_id))
            return redirect(url_for("crispr_finder_result", uuid=job_id))
        else:
            return render_template('crispr_finder.html', form=form, page_title='CRISPR Finder')
    return render_template('crispr_finder.html', form=CrisprFinderForm(), page_title='CRISPR Finder')


@app.route('/crispr_finder/result/<uuid>')
def crispr_finder_result(uuid):
    flag = os.path.join(app.config['UPLOAD_FOLDER'], uuid, 'processing')
    #if not request.script_root:
    #    # this assumes that the 'index' view function handles the path '/'
    #    request.script_root = url_for('/crispr_finder/result/', _external=True)
    return render_template('crispr_finder_result.html', uuid=uuid, dirpath=os.path.join(app.config['UPLOAD_FOLDER'], uuid), processing_flag=os.path.isfile(flag))

@app.route('/antismash')
def antismash():
	return redirect("/antismash")

# TODO : fonction runner qui reccupère les données du formulaire, puis créé le dossier de travail et lance le job dans un thread
#@use_kwargs(CmdSchema())
#def runner():
#    pass

"""
findCRISPRs = FindCRISPRs(args['inputpath'],
                          args['outputpath'],
                          args['k_mer_size_filter'],
                          args['pattern'],
                          args['window_size'],
                          args['allowed_mismatch'],
                          args['spacer_dr_match_limit'],
                          args['min_dr'],
                          args['max_dr'],
                          args['min_spacer_dr_ratio'],
                          args['max_spacer_dr_ratio'],
                          args['first_pass_limit'],
                          args['search_tracrrna'])
"""
