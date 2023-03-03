from .server import broker
import json
import os

@broker.task
def run_sidescan(project_input, project_output):
    exit_code = os.system(f'sidescan-search -i {project_input} -m data/models.pickle -o {project_output}')
    if exit_code == 0:
        return json.loads(open(project_output).read())
    else:
        return {'mgs': 'error'}