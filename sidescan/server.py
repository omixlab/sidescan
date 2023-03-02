import os
from flask import Flask, render_template, session, redirect, url_for, request
from werkzeug.utils import secure_filename
#from flask_wtf import FlaskForm
#from wtforms import StringField, SubmitField
#from wtforms.validators import DataRequired
from flask_bootstrap import Bootstrap

app = Flask(__name__)
bootstrap = Bootstrap(app)
app.config['UPLOAD_PATH'] = '/home/lucasmocellin/project/sidescan/sidescan/static/imgs'

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/', methods=['POST'])
def upload_files():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        #if file_ext not in app.config['UPLOAD_EXTENSIONS']:
          #  abort(400)
        uploaded_file.save(os.path.join(app.config['UPLOAD_PATH'], filename))
        print (filename)
    return render_template('job.html', filename=filename)


#class NameForm(FlaskForm):
  #  molecule = StringField( validators=[DataRequired()])
 #   submit = SubmitField('Submit')



#@app.route('/', methods=['GET', 'POST'])
#def index():
   #form = NameForm()
   #if form.validate_on_submit():
    #    session['name'] = form.molecule.data
    #    return redirect(url_for('index'))
   #return render_template('index.html', form=form, name=session.get('name'))

#@app.route('/upload', methods=['post'])
#def upload():
   # form = NameForm()
  #  if form.validate_on_submit():
  #      print(form.molecule)
  #  return render_template('job.html', molecule=form.molecule.data)


def main():
    app.run(port=5003)
    return

if __name__ == '__main__':
    main()