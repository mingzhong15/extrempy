import numpy as np
import os
import json
import subprocess

from extrempy.constant import kb, J2eV


def fermi_dirac(E, mu, T):
    return 1/(np.exp((E-mu)/(kb*T*J2eV)) + 1)


def parse_machine_json(path, section='model_devi'):
    """Extract Slurm configuration from a dpgen-style ``machine.json``.

    Supports two formats:
      1. **batch** format — ``machine.batch.{slurm_partition, …}`` (dpdispatcher output)
      2. **resources** format — ``resources.{queue_name, number_node, …}``

    Parameters
    ----------
    path : str
        Path to ``machine.json`` (``~`` is expanded).
    section : str
        JSON top-level key to read, e.g. ``'model_devi'`` or ``'fp'``.

    Returns
    -------
    dict
        Keys: ``partition``, ``nodes``, ``ntasks_per_node``, ``wall_time``,
        ``gres``, ``command``, ``source_list``, ``custom_flags``, ``envs``.
        Missing values are ``None``.
    """
    with open(os.path.expanduser(path)) as f:
        machine = json.load(f)

    raw = machine.get(section)
    if isinstance(raw, list):
        raw = raw[0] if raw else {}
    if not raw:
        raw = {}

    res = raw.get('resources', {})
    bat = raw.get('machine', {}).get('batch', {})

    cfg = {}

    if bat:
        # dpdispatcher batch format
        cfg['partition'] = bat.get('slurm_partition')
        cfg['nodes'] = bat.get('slurm_nodes') or res.get('number_node')
        cfg['ntasks_per_node'] = (bat.get('slurm_ntasks_per_node')
                                  or res.get('cpu_per_node'))
        cfg['wall_time'] = bat.get('slurm_time')
        cfg['gres'] = bat.get('slurm_gres')
    else:
        # resources format (common for model_devi on Slurm clusters)
        cfg['partition'] = res.get('queue_name')
        cfg['nodes'] = res.get('number_node')
        cfg['ntasks_per_node'] = res.get('cpu_per_node')
        cfg['wall_time'] = None
        ngpu = res.get('gpu_per_node', 0)
        cfg['gres'] = f'gpu:{ngpu}' if ngpu else None

    cfg['command'] = raw.get('command')
    cfg['source_list'] = res.get('source_list', [])
    cfg['custom_flags'] = res.get('custom_flags', [])
    cfg['envs'] = res.get('envs', {})

    return cfg


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
        """Generate sbatch script from dpgen machine.json FP section."""
        tmpl = parse_machine_json(machine_json_path, section='fp')

        partition = tmpl.get('partition') or ''
        nodes = tmpl.get('nodes') or 1
        ntasks = tmpl.get('ntasks_per_node') or 32
        wall_time = tmpl.get('wall_time') or ''
        command = tmpl.get('command') or 'mpirun vasp_std'
        custom_flags = tmpl.get('custom_flags', [])

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




