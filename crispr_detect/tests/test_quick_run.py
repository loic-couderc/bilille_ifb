import unittest
from io import BytesIO
import tempfile
import os
import shutil
import re
import json
import warnings
from flask import url_for
import crispr_detect


class CrisprFormTestCase(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter('ignore', DeprecationWarning)

        crispr_detect.app.config['TESTING'] = True
        crispr_detect.app.config['DEBUG'] = True
        crispr_detect.app.config['WTF_CSRF_ENABLED'] = False
        self.tempdir = tempfile.TemporaryDirectory()
        crispr_detect.app.config['UPLOAD_FOLDER'] = self.tempdir.name
        crispr_detect.app.config['SERVER_NAME'] = 'localhost'
        self.app = crispr_detect.app.test_client()
        self.ctx = crispr_detect.app.app_context()
        self.ctx.push()

        self.data_form = dict(
            sequence=(BytesIO(b'>header|nACTGCCGTCA'), "filename.fasta"),
        )

    def tearDown(self): pass
        #shutil.rmtree(self.tempdir.name)
        
    

    def test_quick_run_form_no_sequence(self):
        del self.data_form['sequence']
        response = self.app.post(url_for(
            'run_them_all'), content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        assert response.status_code == 400
        assert b'This field is required' in response.data

    def test_antismash_not_available(self):
        old_url = crispr_detect.app.config['WEBSMASH_URL']
        crispr_detect.app.config['WEBSMASH_URL'] = "http://fake_url/"
        response = self.app.post(url_for(
            'run_them_all'), content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        crispr_detect.app.config['WEBSMASH_URL'] = old_url
        assert b'503' in response.data

    def test_quick_run_form_sequence_file_extension(self):
        self.data_form['sequence'] = (
            BytesIO(b'my file contents'), "filename.wrong")

        response = self.app.post(url_for(
            'run_them_all'), content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        assert response.status_code == 400
        assert b'Fasta only' in response.data

    def test_quick_run_form_success(self):
        response = self.app.post(url_for(
            'run_them_all'), content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        assert response.status_code == 200
        assert str(response.data).count("SUCCESS") == 2



if __name__ == '__main__':
    unittest.main()
