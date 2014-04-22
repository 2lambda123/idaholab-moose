from InputParameters import InputParameters
import os, sys, shutil

MOOSE_DIR = os.path.abspath(os.path.join('..', 'moose'))
FRAMEWORK_DIR = os.path.abspath(os.path.join('..', 'moose', 'framework'))
#### See if MOOSE_DIR is already in the environment instead
if os.environ.has_key("MOOSE_DIR"):
  MOOSE_DIR = os.environ['MOOSE_DIR']
  FRAMEWORK_DIR = os.path.join(MOOSE_DIR, 'framework')
if os.environ.has_key("FRAMEWORK_DIR"):
  FRAMEWORK_DIR = os.environ['FRAMEWORK_DIR']

class Job(object):
  def getValidParams():
    params = InputParameters()
    params.addRequiredParam('type', "The type of test of Tester to create for this test.")
    params.addParam('template_script', FRAMEWORK_DIR + '/scripts/ClusterLauncher/pbs_submit.sh', "The template job script to use.")
    return params
  getValidParams = staticmethod(getValidParams)

  def __init__(self, name, params):
    self.specs = params

  # Called from the current directory to copy files (usually from the parent)
  def copyFiles(self, job_file):
    for file in os.listdir('../'):
      if os.path.isfile('../' + file) and file != job_file:
        shutil.copy('../' + file, '.')

  # Called to prepare a job script if necessary
  def prepareJobScript(self):
    return

  # Called to launch the job
  def launch(self):
    return
