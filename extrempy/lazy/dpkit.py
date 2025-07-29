
from extrempy.constant import *
import json

def _generate_training_params(type_map, trainig_data, 
                              is_ele=True, stop_batch=1000000):

    training_param = {
        "model": {
            "type_map": type_map,
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
                "systems": trainig_data,
                "set_prefix": "set",
                "batch_size": 1
            }
        },
        "_comment": "that's all"
    }

    
    if is_ele:
        training_param["model"]["fitting_net"]["numb_fparam"] = 1
    else:
        training_param["model"]["fitting_net"]["numb_fparam"] = 0
        
    return training_param

def _generate_job_params(work_dir, job_file, job_name):

    with open(job_file, 'r') as f:
        job_param = json.load(f)

    job_param["job_name"] = job_name

    with open(os.path.join(work_dir, 'job.json'), 'w') as f:
        json.dump(job_param, f, indent=4)
    