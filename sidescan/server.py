import os
from flask import Flask, render_template, session, redirect, url_for, request, jsonify
from werkzeug.utils import secure_filename
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from celery.result import AsyncResult
import celery
import uuid
import json


app = Flask(__name__)
broker = celery.Celery(broker='redis://localhost:6379', backend='redis://localhost:6379')
bootstrap = Bootstrap(app)
app.config['UPLOAD_PATH'] = 'sidescan/static/uploads'

db = SQLAlchemy()
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///:memory:"
db.init_app(app)

class Job(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    email = db.Column(db.String)
    status = db.Column(db.String)
    celery_job_id = db.Column(db.String)

with app.app_context():
    db.create_all()

from sidescan.worker import run_sidescan

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/', methods=['POST'])
def upload_files():

    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)

    job = Job()
    db.session.add(job)
    db.session.commit()

    project_id = str(job.id)
    project_directory = os.path.join(app.config['UPLOAD_PATH'], project_id)
    project_output = os.path.join(project_directory, 'predictions.json')

    print(project_directory)

    os.system(f"mkdir -p {project_directory}")
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        
        project_input = os.path.join(project_directory, filename)
        uploaded_file.save(project_input)

        async_result = run_sidescan.apply_async((project_input, project_output))
        job.celery_job_id = str(async_result.task_id)
        job.status = 'QUEUED'
        db.session.commit()

    return redirect(f'/jobs/{job.id}')

@app.route('/jobs/<int:job_id>')
def jobs(job_id):
    job = db.get_or_404(Job, job_id)
    result = AsyncResult(job.celery_job_id, app=broker)
    job.status = result.status
    db.session.commit()
    if result.status == 'SUCCESS':
        return render_template('job.html', result_status=result.status, job_id=job_id, result_data=result.get())
    else:
        return render_template('job.html', result_status=result.status, job_id=job_id)

def main():
    app.run(port=5009)
    return

if __name__ == '__main__':
    main()