import numpy as np
import os
import json

from extrempy.constant import kb, J2eV


def fermi_dirac(E, mu, T):
    return 1/(np.exp((E-mu)/(kb*T*J2eV)) + 1)


class InputGenerator:

    def __init__(self, work_path):

        self.work_path = work_path

    def generate_submit(self, 
                        job_template_path,
                        job_name,
                        platform = 'bh',
                        job_group_id=None):

        self.platform = platform

        if self.platform == 'bh':

            with open(job_template_path, 'r') as f:
                job_param = json.load(f)

            self.job_name = job_name
            job_param["job_name"] = self.job_name

            if job_group_id is not None:
                job_param["job_group_id"] = job_group_id

            with open(os.path.join(self.work_path, 'job.json'), 'w') as f:
                json.dump(job_param, f, indent=4)

    def submit(self):
        cwd = os.getcwd()
        os.chdir(self.work_path)
        if self.platform == 'bh':
            try:
                os.system('mkdir  ../'+self.job_name)
            except:
                pass
            pwd = 'bohr job submit -i job.json -p ./ -r ../'+self.job_name
            os.system(pwd)
        os.chdir(cwd)




