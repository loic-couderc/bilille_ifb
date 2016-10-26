from flask_wtf import FlaskForm as Form
from wtforms import FileField, TextAreaField, StringField, IntegerField, FloatField, SelectField
from wtforms.validators import DataRequired

class CrisprFinderForm(Form):

    sequence = TextAreaField('Sequence', validators=[
        DataRequired('Please provide a valid DNA sequence')])
    #sequence = FileField( validators=[
    #    DataRequired()])
    k_mer_size_filter = IntegerField(default=3, validators=[
        DataRequired()])
    pattern = StringField(default="####_####", validators=[
        DataRequired()])
    window_size = IntegerField(default=200)
    allowed_mismatch = IntegerField(default=1)
    spacer_dr_match_limit = IntegerField(default=20)
    min_dr = IntegerField(default=23)
    max_dr = IntegerField(default=55)
    min_spacer_dr_ratio = FloatField(default=0.6)
    max_spacer_dr_ratio = FloatField(default=2.5)
    first_pass_limit = IntegerField(default=200)
    search_tracrrna = SelectField(choices=[('False','no'),('True','yes')], default='False')
