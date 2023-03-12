import os
from flask import Flask, render_template, session, redirect, url_for, request, jsonify
from werkzeug.utils import secure_filename
#from flask_wtf import FlaskForm
#from wtforms import StringField, SubmitField
#from wtforms.validators import DataRequired
from flask_bootstrap import Bootstrap
from celery.result import AsyncResult
import celery
import uuid
import json


app = Flask(__name__)
broker = celery.Celery(broker='redis://localhost:6379', backend='redis://localhost:6379')
bootstrap = Bootstrap(app)
app.config['UPLOAD_PATH'] = 'sidescan/static/uploads'

from sidescan.worker import run_sidescan

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/', methods=['POST'])
def upload_files():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)

    project_id = str(uuid.uuid4())
    project_directory = os.path.join(app.config['UPLOAD_PATH'], project_id)
    project_output = os.path.join(project_directory, 'predictions.json')

    os.system(f"mkdir -p {project_directory}")
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        
        project_input = os.path.join(project_directory, filename)
        uploaded_file.save(project_input)

        async_result = run_sidescan.apply_async((project_input, project_output))

    return redirect(f'/jobs/{async_result.task_id}')

@app.route('/jobs/<job_id>')
def jobs(job_id):
    result = AsyncResult(job_id, app=broker)
    if result.status == 'SUCCESS':
        return render_template('job.html', result_status=result.status, job_id=job_id, result_data=result.get())
    else:
        return render_template('job.html', result_status=result.status, job_id=job_id)



def main():
    app.run(port=5008)
    return

if __name__ == '__main__':
    main()