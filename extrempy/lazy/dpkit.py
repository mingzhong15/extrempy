
from extrempy.constant import *
import json

class DPKITParamGenerator:

    def __init__(self, work_dir, type_map, is_ele = True):

        self.work_dir = work_dir
        self.type_map = type_map
        self.is_ele = is_ele

    def _get_training_set(self, set_dir, prefix='*'):

        self.set_dir = set_dir
        set_list = glob.glob(os.path.join(set_dir, prefix, '*'))

        # 构建出相对路径
        self.dataset_list = []

        for set_path in set_list:
            # Calculate relative path from PWD to set_path
            relative_path = os.path.relpath(set_path, self.work_dir)
            self.dataset_list.append(relative_path)

    def _generate_training_params(self, stop_batch=1000000):

        training_param = {
            "model": {
                "type_map": self.type_map,
                "descriptor": {
                    "type":         "se_e2_a",
                    "rcut_smth": 	1.8,
                    "rcut": 	    6.0,
                    "neuron":	    [25, 50, 100],
                    "type_one_side":True,
                    "resnet_dt": 	False,
                    "axis_neuron": 	8,
                    "seed": 	    np.random.randint(0, 100000)
                },
                "fitting_net": {
                    "neuron":		[240, 240, 240],
                    "resnet_dt":    True,
                    "numb_fparam":	0,
                    "seed":         np.random.randint(0, 100000)
                }
            },
            "learning_rate": {
                "type": "exp",
                "start_lr": 1e-3,
                "stop_lr":1e-8
            },
            "loss": {
                "start_pref_e": 0.2,
                "limit_pref_e": 1,
                "start_pref_f": 1000,
                "limit_pref_f": 1,
                "start_pref_v": 0.2,
                "limit_pref_v": 1
            },
            "training": {
                "stop_batch": stop_batch,
                "seed":     np.random.randint(0, 100000),
                "disp_file": "lcurve.out",
                "disp_freq": 100,
                "numb_test": 10,
                "save_freq": 1000,
                "save_ckpt": "model.ckpt",
                "disp_training": True,
                "time_training": True,
                "profiling": False,
                "profiling_file": "timeline.json",
                "training_data": {
                    "systems": self.dataset_list,
                    "set_prefix": "set",
                    "batch_size": 1
                }
            },
            "_comment": "that's all"
        }

        
        if self.is_ele:
            training_param["model"]["fitting_net"]["numb_fparam"] = 1
        else:
            training_param["model"]["fitting_net"]["numb_fparam"] = 0
            
        with open( os.path.join(self.work_dir, 'input.json'), 'w') as f:
            json.dump(training_param, f, indent=4)

    def _generate_job_params(self, job_file, job_name):

        with open(job_file, 'r') as f:
            job_param = json.load(f)

        if self.is_ele:
            ele_label = 'f'
        else:
            ele_label = 'no_f'

        self.job_name = job_name + '_' + ele_label
        job_param["job_name"] = self.job_name

        with open(os.path.join(self.work_dir, 'job.json'), 'w') as f:
            json.dump(job_param, f, indent=4)

    def _submit_job(self, platform='bh'):

        if platform == 'bh':
            os.chdir(self.work_dir)

            try:
                os.system('mkdir  ../'+self.job_name)
            except:
                pass

            pwd = 'lbg job submit -i job.json -p ./ -r ../'+self.job_name
            os.system(pwd)
        else:
            print('platform not supported')