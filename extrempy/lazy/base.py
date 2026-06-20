import numpy as np
import os
import json
import subprocess

from extrempy.constant import kb, J2eV


def fermi_dirac(E, mu, T):
    return 1/(np.exp((E-mu)/(kb*T*J2eV)) + 1)


class InputGenerator:

    def __init__(self, work_path):

        self.work_path = work_path

    def generate_submit(self,
                        job_template_path,
                        job_name,
                        platform='bh',
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

        elif self.platform == 'slurm':
            self._generate_slurm_submit(job_template_path, job_name)

    def _generate_slurm_submit(self, machine_json_path, job_name):
        """Generate sbatch script from dpgen machine.json FP section.

        Supports two formats:
          1. dpdispatcher 'batch' dict: machine.batch.slurm_partition, ...
          2. dpdispatcher 'resources' dict: resources.queue_name, ...
        """
        with open(machine_json_path, 'r') as f:
            machine = json.load(f)
        fp_conf = machine.get('fp', [{}])[0]
        resources = fp_conf.get('resources', {})
        batch = fp_conf.get('machine', {}).get('batch', {})
        command = fp_conf.get('command', 'mpirun vasp_std')
        custom_flags = resources.get('custom_flags', [])

        if batch:
            partition = batch.get('slurm_partition', '')
            nodes = batch.get('slurm_nodes', 1)
            ntasks = batch.get('slurm_ntasks_per_node', 32)
            wall_time = batch.get('slurm_time', '24:00:00')
            gres = batch.get('slurm_gres', '')
            extra = batch.get('slurm_args', '')
        else:
            partition = resources.get('queue_name', '')
            nodes = resources.get('number_node', 1)
            ntasks = resources.get('cpu_per_node', 32)
            wall_time = ''

        lines = ['#!/bin/bash']
        lines.append(f'#SBATCH -J {job_name}')
        if partition:
            lines.append(f'#SBATCH -p {partition}')
        lines.append(f'#SBATCH -N {nodes}')
        lines.append(f'#SBATCH --ntasks-per-node={ntasks}')
        if wall_time:
            lines.append(f'#SBATCH -t {wall_time}')
        for flag in custom_flags:
            flag = flag.strip()
            if flag.startswith('#SBATCH') and '--job-name' not in flag:
                lines.append(flag)
        lines.append('')
        lines.append(f'cd {self.work_path}')
        lines.append('')
        lines.append(command)

        sbatch_path = os.path.join(self.work_path, 'job.sbatch')
        with open(sbatch_path, 'w') as f:
            f.write('\n'.join(lines) + '\n')

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
        elif self.platform == 'slurm':
            try:
                result = subprocess.run(
                    ['sbatch', 'job.sbatch'],
                    capture_output=True, text=True, check=True)
                print(f"  sbatch submitted: {result.stdout.strip()}")
            except subprocess.CalledProcessError as e:
                print(f"  ERROR: sbatch failed: {e.stderr}")
            except FileNotFoundError:
                print("  ERROR: 'sbatch' not found. Is Slurm installed?")
        os.chdir(cwd)




