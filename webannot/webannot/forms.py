import io
from flask_wtf import FlaskForm as Form
from wtforms import TextAreaField, StringField, IntegerField, FloatField, SelectField
from wtforms.validators import DataRequired
from flask_wtf.file import FileField, FileAllowed, FileRequired
from Bio import SeqIO
from Bio.Alphabet import generic_dna

class CrisprFinderForm(Form):

    sequence = FileField(validators=[FileAllowed(
        ['fasta', 'fsa', 'fa'], 'Fasta only (.fa/.fsa/.fasta)!')
    ])
    sequence2 = TextAreaField()

    k_mer_size_filter = IntegerField(default=3,
                                     description='is used to compare segments of DRs and spacers')
    pattern = StringField(default="####_####",
                          description='is a combination of \'#\' (character) and \'_\' (space) used as seeds during second pass. First pass used continuous seed of length k_mer_size (number of \'#\')')
    window_size = IntegerField(default=200,
                               description='defines size of window in which k_mers will be searched for')
    allowed_mismatch = IntegerField(default=1,
                                    description='defines a number of mismatch allowed in repeat')
    spacer_dr_match_limit = IntegerField(default=20,
                                         description='is a maximum number of matches of length k_mer_size_filter per spacer between DR and spacer',
                                         label='Spacer DR Match Limit')
    min_dr = IntegerField(default=23,
                          description='defines minimum length in bp of repeat part',
                          label='Min DR')
    max_dr = IntegerField(default=55,
                          description='defines maximum length in bp of repeat part',
                          label='Max DR')
    min_spacer_dr_ratio = FloatField(default=0.6,
                                     description='is minimum quotient allowed length of spacer versus DR',
                                     label='Min Spacer DR Ratio')
    max_spacer_dr_ratio = FloatField(default=2.5,
                                     description='is maximum quotient allowed length of spacer versus DR',
                                     label='Max Spacer DR Ratio')
    first_pass_limit = IntegerField(default=200,
                                    description='is maximum allowed distance between two regions with repeats')
    search_tracrrna = SelectField(choices=[('False', 'no'), ('True', 'yes')], default='False', description='',
                                  label='Search tracrRNA')

    def validate(self):
        rv = Form.validate(self)

        if not self.sequence.data.filename and not self.sequence2.data:
            message = "A sequence have to be provided. Please, upload a file or past a sequence."
            self.sequence.errors.append(message)
            self.sequence2.errors.append(message)
            return False
        elif self.sequence.data.filename and self.sequence2.data:
            message = "Please, choose a file OR past your sequence"
            self.sequence.errors.append(message)
            self.sequence2.errors.append(message)
            return False

        if self.sequence.data.filename:
            stream = io.StringIO(self.sequence.data.stream.read().decode("utf-8"))
            self.sequence.data.stream.seek(0)
            message = self._valid_fasta(stream)
            if message is not None:
                self.sequence.errors.append(message)
                return False
        elif self.sequence2:
            stream = io.StringIO(self.sequence2.data)
            message = self._valid_fasta(stream)
            if message is not None:
                self.sequence2.errors.append(message)
                return False
        return rv

    def _valid_fasta(self, stream):
        try:
            seq_count = 0
            for record in SeqIO.parse(stream, "fasta"):
                seq_count+=1
                if record.seq == "":
                    return "Empty sequence found"
                if seq_count > 1:
                    return "Only one sequence have to be provided"
            if seq_count == 0:
                return "No valid sequence(s) provided"
        except Exception as e:
            return "Can't parse fasta sequence(s): %s" % e
