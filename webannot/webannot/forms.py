from flask_wtf import FlaskForm as Form
from wtforms import TextAreaField, StringField, IntegerField, FloatField, SelectField
from wtforms.validators import DataRequired
from flask_wtf.file import FileField, FileAllowed, FileRequired


class CrisprFinderForm(Form):

    # sequence = TextAreaField('Sequence', validators=[
    #    DataRequired('Please provide a valid DNA sequence')])
    sequence = FileField(validators=[
        FileRequired(), FileAllowed(
            ['fasta', 'fsa', 'fa'], 'Fasta only (.fa/.fsa/.fasta)!')
    ])
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
