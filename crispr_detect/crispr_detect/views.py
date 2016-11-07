from crispr_detect import app
from flask import render_template, abort, flash, redirect, request, redirect, url_for
from jinja2 import TemplateNotFound
from .forms import CrisprFinderForm
from marshmallow import Schema, fields
from webargs.flaskparser import parser

import threading
import os, sys
import uuid
from .crispr_detect import FindCRISPRs


class CmdSchema(Schema):
    sequence = fields.Str(required=True)

    k_mer_size_filter = fields.Integer(required=True)
    pattern = fields.Str(required=True)
    window_size = fields.Integer(required=True)
    allowed_mismatch = fields.Integer(required=True)
    spacer_dr_match_limit = fields.Integer(required=True)
    min_dr = fields.Integer(required=True)
    max_dr = fields.Integer(required=True)
    min_spacer_dr_ratio = fields.Float(required=True)
    max_spacer_dr_ratio = fields.Float(required=True)
    first_pass_limit = fields.Integer(required=True)
    search_tracrrna = fields.Boolean(required=True)

    class Meta:
        strict = True


# @copy_current_request_context
def crispr_finder_job(**args):
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


def get_dirpath(dirname):
    #TODO: Take into account the user session for remote access
    return (os.path.join(app.config['UPLOAD_FOLDER'], dirname))


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
        form = CrisprFinderForm(request.form)
        if form.validate():
            # print request.form
            args = parser.parse(CmdSchema, request)
            # print args
            # args = {k: v for k, v in args.iteritems() if v != u''}
            job_id = str(uuid.uuid4())
            args['job_id'] = job_id
            outputpath = get_dirpath(job_id)
            os.mkdir(outputpath)
            inputname = 'input.dna'
            inputpath = os.path.join(outputpath, inputname)
            with open(inputpath, 'w') as f:
                f.write(args['sequence'])
            # args.pop('sequence', None)
            args['inputpath'] = inputpath
            args['outputpath'] = outputpath
            run_crispr = threading.Thread(target=crispr_finder_job, kwargs=args)
            run_crispr.start()
            print(url_for("crispr_finder_result", uuid=job_id))
            return redirect(url_for("crispr_finder_result", uuid=job_id))
        else:
            return render_template('crispr_finder.html', form=form, page_title='CRISPR Finder')
    return render_template('crispr_finder.html', form=CrisprFinderForm(), page_title='CRISPR Finder')


@app.route('/crispr_finder/result/<uuid>')
def crispr_finder_result(uuid):
    flag = os.path.join(get_dirpath(uuid), 'processing')
    #if not request.script_root:
    #    # this assumes that the 'index' view function handles the path '/'
    #    request.script_root = url_for('/crispr_finder/result/', _external=True)
    return render_template('crispr_finder_result.html', uuid=uuid, dirpath=get_dirpath(uuid),processing_flag=os.path.isfile(flag))

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
