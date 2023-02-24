from flask import Flask, render_template, session, redirect, url_for, flash
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired

app = Flask(__name__)
app.config['SECRET_KEY'] = 'chaveboa'

class NameForm(FlaskForm):
    molecule = StringField( validators=[DataRequired()])
    submit = SubmitField('Submit')

@app.route('/', methods=['GET', 'POST'])
def index():
    form = NameForm()
    if form.validate_on_submit():
        session['name'] = form.molecule.data
        return redirect(url_for('index'))
    return render_template('index.html', form=form, name=session.get('name'))

@app.route('/upload', methods=['post'])
def upload():
    form = NameForm()
    if form.validate_on_submit():
        print(form.molecule)
    return render_template('job.html', molecule=form.molecule.data)


def main():
    app.run()
    return

if __name__ == '__main__':
    main()