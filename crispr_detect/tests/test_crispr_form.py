import unittest
from io import BytesIO
import tempfile
import os
import re
import json
import warnings
from flask import url_for
import crispr_detect



class CrisprFormTestCase(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter('ignore', DeprecationWarning)

        crispr_detect.app.config['TESTING'] = True
        crispr_detect.app.config['WTF_CSRF_ENABLED'] = False
        self.tempdir = tempfile.TemporaryDirectory()
        crispr_detect.app.config['UPLOAD_FOLDER'] = self.tempdir.name
        crispr_detect.app.config['SERVER_NAME'] = 'localhost'
        self.app = crispr_detect.app.test_client()
        self.ctx = crispr_detect.app.app_context()
        self.ctx.push()

        self.data_form = dict(
            sequence=(BytesIO(b'my file contents'), "filename.fasta"),
            k_mer_size_filter=3,
            pattern="####_####",
            window_size=200,
            allowed_mismatch=1,
            spacer_dr_match_limit=20,
            min_dr=23,
            max_dr=25,
            min_spacer_dr_ratio=0.6,
            max_spacer_dr_ratio=2.5,
            first_pass_limit=200,
            search_tracrrna='False'
        )

    def tearDown(self):
        self.tempdir.cleanup()

    def test_crispr_form_no_sequence(self):
        del self.data_form['sequence']
        response = self.app.post(url_for(
            'crispr_finder'), content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        #response = self.app.post('/crispr_finder/', content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        assert response.status_code == 400
        assert b'This field is required' in response.data

    def test_crispr_form_sequence_file_extension(self):
        self.data_form['sequence'] = (
            BytesIO(b'my file contents'), "filename.wrong")

        response = self.app.post(url_for(
            'crispr_finder'), content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        #response = self.app.post('/crispr_finder/', content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        assert response.status_code == 400
        assert b'Fasta only' in response.data

    def parse_uuid(self, response):
        # Crispr Finder  : 84771680-68bd-4807-9ded-cf7fd7324364
        match = re.search(
            r'CRISPRDetect\s+:\s+([\-\w)]+)', str(response.data))
        return match.group(1)

    def parse_cmd_parameters(self, string):
        match = re.search(r"(\[.*\])", string)
        return eval(match.group(1))

    def send_crispr_form(self, data):
        response = self.app.post(url_for(
            'crispr_finder'), content_type='multipart/form-data', data=data, follow_redirects=True)
        #response = self.app.post('/crispr_finder/', content_type='multipart/form-data', data=self.data_form, follow_redirects=True)
        assert response.status_code == 200

        uuid = self.parse_uuid(response)

        while True:
            response = self.app.get(url_for('processing'), query_string={
                'uuid': uuid}, content_type='application/json')
            #response = self.app.get('/_processing', query_string={ 'uuid' : uuid }, content_type='application/json')
            processing = json.loads(response.get_data(as_text=True))
            if not processing:
                break
        with open(os.path.join(crispr_detect.app.config['UPLOAD_FOLDER'], uuid, 'stdout'), 'r') as stdout_file:
            params = self.parse_cmd_parameters(str(stdout_file.read()))
            return params

    def test_k_mer_size_filter_crispr_form(self):
        self.data_form['k_mer_size_filter'] = 4
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            4, '####_####', 200, 1, 20, 23, 25, 0.6, 2.5, 200, False]

    def test_pattern_crispr_form(self):
        self.data_form['pattern'] = '#####_#####'
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '#####_#####', 200, 1, 20, 23, 25, 0.6, 2.5, 200, False]

    def test_window_size_filter_crispr_form(self):
        self.data_form['window_size'] = 201
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 201, 1, 20, 23, 25, 0.6, 2.5, 200, False]

    def test_allowed_mismatch_crispr_form(self):
        self.data_form['allowed_mismatch'] = 2
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 2, 20, 23, 25, 0.6, 2.5, 200, False]

    def test_spacer_dr_match_limit_crispr_form(self):
        self.data_form['spacer_dr_match_limit'] = 21
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 1, 21, 23, 25, 0.6, 2.5, 200, False]

    def test_min_dr_crispr_form(self):
        self.data_form['min_dr'] = 22
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 1, 20, 22, 25, 0.6, 2.5, 200, False]

    def test_max_dr_crispr_form(self):
        self.data_form['max_dr'] = 27
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 1, 20, 23, 27, 0.6, 2.5, 200, False]

    def test_min_spacer_dr_ratio_crispr_form(self):
        self.data_form['min_spacer_dr_ratio'] = 0.5
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 1, 20, 23, 25, 0.5, 2.5, 200, False]

    def test_max_spacer_dr_ratio_crispr_form(self):
        self.data_form['max_spacer_dr_ratio'] = 2.7
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 1, 20, 23, 25, 0.6, 2.7, 200, False]

    def test_first_pass_limit_crispr_form(self):
        self.data_form['first_pass_limit'] = 199
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 1, 20, 23, 25, 0.6, 2.5, 199, False]

    def test_search_tracrrna_crispr_form(self):
        self.data_form['search_tracrrna'] = True
        params = self.send_crispr_form(self.data_form)
        assert params[2:] == [
            3, '####_####', 200, 1, 20, 23, 25, 0.6, 2.5, 200, True]

if __name__ == '__main__':
    unittest.main()
