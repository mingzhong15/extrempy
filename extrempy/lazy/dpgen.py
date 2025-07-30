from extrempy.constant import *
import json
from .lib import _get_mass_map

def _generate_dpgen_template():

    training_param = {
        "model": {
            "type_map": ["A"],
            "descriptor": {
                "type": "se_e2_a",
                "rcut_smth": 	1.8,
                "rcut": 	6.0,
                "neuron":	[25, 50, 100],
                "type_one_side": True,
                "resnet_dt": 	False,
                "axis_neuron": 	8,
                "seed": 	0
            },
            "fitting_net": {
                "neuron":		[240, 240, 240],
                "resnet_dt": True,
                "numb_fparam":	0,
                "seed": 0
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
            "stop_batch": 200000,
            "disp_file": "lcurve.out",
            "disp_freq": 100,
            "numb_test": 10,
            "save_freq": 1000,
            "save_ckpt": "model.ckpt",
            "disp_training": True,
            "time_training": True,
            "profiling": False,
            "profiling_file": "timeline.json",
            "_comment": "that's all"
        }
    }

    template_param = {
        "type_map": ["A"],
        "mass_map": [1.0],
        "init_data_prefix": "./",
        "init_data_sys": [],
        "sys_configs_prefix": "./",
        "sys_configs": [
            ["phase-1.vasp"]
        ],
        "_comment": " that's all ",
        "numb_models": 4,
        "train_param": "input.json",
        "default_training_param": training_param,

        "training_init_model": False,
        "training_reuse_iter": 1,
        "training_reuse_numb_steps": 100000,
        "training_reuse_start_lr": 1e-4,
        "use_relative": True,
        "model_devi_f_avg_relative": False,
        "epsilon": 0,
        "model_devi_dt": 0.0005,
        "model_devi_skip":0,
        "model_devi_f_trust_lo": 0.1,
        "model_devi_f_trust_hi": 1.0,
        "model_devi_clean_traj": False,
        "model_devi_jobs": [],
        "fp_style": "vasp",
        "use_ele_temp":1,
        "shuffle_poscar": False,
        "ratio_failed": 0.5,
        "fp_task_max": 120,
        "fp_task_min": 0,
        "fp_pp_path": "./",
        "fp_pp_files": ["POTCAR.A"],
        "fp_incar": "INCAR.A"
    }

    return template_param


class DPGENParamGenerator:

    def __init__(self, type_map, json_file=None):

        if json_file is None:
            self.jparam = _generate_dpgen_template()
        else:
            with open(json_file, 'r') as f:
                self.jparam = json.load(f)

        self.type_map = type_map

        self.jparam["type_map"] = self.type_map
        self.jparam["mass_map"] = _get_mass_map(self.type_map)

        self.jparam["default_training_param"]["model"]["type_map"] = self.type_map

    def _set_init_data(self, SET_DIR, prefix='*'):

        set_list = glob.glob(os.path.join(SET_DIR, prefix))

        SET_LIST = []

        for set_path in set_list:
            # Calculate relative path from PWD to set_path
            relative_path = os.path.relpath(set_path, SET_DIR)
            SET_LIST.append(relative_path)

        self.jparam["init_data_prefix"] = SET_DIR
        self.jparam['init_data_sys'] = SET_LIST

    def _set_sys_configs(self, CONF_DIR, CONF_LIST):

        self.jparam["sys_configs_prefix"] = CONF_DIR
        self.jparam["sys_configs"] = CONF_LIST
    
    def _set_model_devi_settings(self, dt = 0.001, 
                                 f_trust = [0.1, 1.0],
                                 is_train_init=False, 
                                 is_relative=True,
                                 epsilon=1.0):

        self.jparam["model_devi_dt"] = dt
        self.jparam["model_devi_skip"] = 0
        self.jparam["model_devi_f_trust_lo"] = f_trust[0]
        self.jparam["model_devi_f_trust_hi"] = f_trust[1]

        self.jparam["training_init_model"] = is_train_init

        if is_train_init:
            self.jparam["training_reuse_iter"] = 1
            self.jparam["training_reuse_numb_steps"] = 100000
            self.jparam["training_reuse_start_lr"] = 1e-4

        self.jparam["use_relative"] = is_relative
        if is_relative:
            self.jparam["epsilon"] = epsilon
            self.jparam["use_relative"] = True
            
    def _set_model_traninig_settings(self, stop_batch=200000, is_ele_temp=True):
        
        self.jparam["default_training_param"]["training"]["stop_batch"] = stop_batch

        if is_ele_temp:
            self.jparam["use_ele_temp"] = 1
            self.jparam["default_training_param"]["model"]["fitting_net"]["numb_fparam"] = 1
        else:
            self.jparam["use_ele_temp"] = 0
            self.jparam["default_training_param"]["model"]["fitting_net"]["numb_fparam"] = 0

    def _set_fp_settings(self, fp_pp_path, fp_pp_files, fp_incar):

        self.jparam["fp_accurate_threshold"] = 0.95
        self.jparam["fp_accurate_soft_threshold"] = 0.0

        self.jparam["fp_pp_path"] = fp_pp_path
        self.jparam["fp_pp_files"] = fp_pp_files
        self.jparam["fp_incar"] = os.path.join(fp_pp_path, fp_incar)

    def _set_model_devi_jobs(self, sys_idx, Tmin, Tmax, Pmin, Pmax,
                             is_pimd=False,
                             ensemble='npt', 
                             numb_iters=5,
                             delta_T=2000,
                             init_steps=1000,
                             trj_freq=20,
                             numb_frame_per_iter_per_PT = 10,
                             nbeads=0):

        T_list = _generate_temp_list(Tmin, Tmax)
        p_list = (np.array(_generate_pres_list(Pmin, Pmax))*1e4).tolist()

        self.T_list = T_list
        self.p_list = p_list
        self.sys_list = sys_idx

        T_ranges = split_temperature_range(T_list, delta_T)

        init_iters_numb = len(self.jparam["model_devi_jobs"])
        print('Existing %.d iterations '%init_iters_numb)
        
        real_idx = init_iters_numb

        for range_idx, current_T_list in enumerate(T_ranges):

            print(f"\ntemperatures are divided into {range_idx + 1}/{len(T_ranges)} ranges")
            print(f"for temperatures ranges: {min(current_T_list):.2f}K - {max(current_T_list):.2f}K")
            
            for iter_idx in range(numb_iters):

                if real_idx < len(init_steps):
                    nsteps = init_steps[real_idx]
                else:
                    nsteps = init_steps[-1] * pow(2, iter_idx - len(init_steps)+1)

                # 创建新的字典

                if ensemble == 'nvt' or ensemble == 'npt':
             
                    new_job = {
                        "sys_idx": sys_idx,
                        "temps": current_T_list,
                    }

                    strs_1 = f" Iter. {real_idx}: nsteps = {nsteps},\n"
                    strs_2 = f" temp ({len(current_T_list)}) = {current_T_list},\n "

                    if ensemble == 'nvt':
                        pass
                    elif ensemble == 'npt':
                        new_job['press'] = p_list
                        strs_2 += f" press ({len(p_list)}) = {p_list} \n"
                    
                    if not is_pimd:
                        pass
                    else:
                        new_job['nbeads'] = nbeads
                        strs_1 += f" nbeads = {nbeads} (PIMD) \n"

                    print(strs_1 + strs_2)

                    new_job["ensemble"] = ensemble


                new_job["trj_freq"] = trj_freq
                new_job["nsteps"] = nsteps
                new_job["_idx"] = real_idx

                self.jparam["model_devi_jobs"].append(new_job)
                real_idx += 1

                #print(self.jparam["model_devi_jobs"][-1]['temps'], self.jparam["model_devi_jobs"][-1]['press'])

        print('total number of tasks per iteration: ', len(T_list) * len(p_list) * numb_frame_per_iter_per_PT ) 

        self.jparam["fp_task_max"] = len(T_list) * len(p_list) * numb_frame_per_iter_per_PT
        self.jparam["fp_task_min"] = len(T_list) * len(p_list) * 1

    
def _generate_dpgen_machine_from_file(machine_file, prefix, is_pimd=False, nbeads=8):

    with open(machine_file, 'r') as f:
        machine_param = json.load(f)

    machine_param['train'][0]['machine']['remote_profile']['input_data']['job_name'] = prefix + '_dpgen_dp'
    machine_param['model_devi'][0]['machine']['remote_profile']['input_data']['job_name'] = prefix + '_dpgen_md'
    machine_param['fp'][0]['machine']['remote_profile']['input_data']['job_name'] = prefix + '_dpgen_fp'
    
    if not is_pimd:
        machine_param['model_devi'][0]['command'] = 'lmp -i input.lammps -v restart 0'
    else:
        machine_param['model_devi'][0]['command'] = 'mpirun --allow-run-as-root -np %.d lmp'%(nbeads)

        if nbeads == 8:
            gpu_type = 'c8_m32_1 * NVIDIA V100'
        elif nbeads == 16:
            gpu_type = 'c16_m62_1 * NVIDIA T4'
        elif nbeads == 32:
            gpu_type = 'c32_m64_cpu'
        else:
            raise ValueError('nbeads must be 8, 16, or 32')

        machine_param['model_devi'][0]['machine']['remote_profile']['input_data']['scass_type'] = gpu_type
        machine_param['model_devi'][0]['machine']['remote_profile']['input_data']['image_address'] = 'registry.dp.tech/dptech/dpmd:2.2.8-cuda12.0'
    
    return machine_param


def split_temperature_range(T_list, delta_T):
    """
    将温度列表按照给定的温度间隔划分成多个子列表
    
    参数:
    T_list: 原始温度列表
    delta_T: 温度区间大小
    
    返回:
    list of lists: 划分后的温度区间列表
    """
    if not T_list:
        return []
    
    # 确保温度列表是排序的
    T_list = sorted(T_list)
    T_ranges = []
    current_range = []
    
    for temp in T_list:
        if not current_range:
            current_range.append(temp)
        elif temp - current_range[0] <= delta_T:
            current_range.append(temp)
        else:
            T_ranges.append(current_range)
            current_range = [temp]
    
    if current_range:
        T_ranges.append(current_range)
    
    return T_ranges

def _generate_temp_list(xmin, xmax, scale=0.2, delta_x=50):

    x_list = [xmin]
    current_x = xmin
    while(True):
        current_x = x_list[-1] +  max(delta_x, int(x_list[-1] * scale /10)*10 )
        x_list.append( current_x )
        if current_x > xmax:
            break
    return x_list

def _generate_pres_list(xmin, xmax, delta_min=1, Nx_max=100):

    scale = xmax / xmin
    x_list = [xmin]

    if scale >= 100:
        print('generating from Log Mode ...')
        
        current_x = xmin
        while(True):
            current_x = x_list[-1] * 10
            if current_x > delta_min or current_x > xmax:
                break
            x_list.append( current_x )
            
    if x_list[-1] < xmax:
        print('generating from Linear Mode ...')
        
        delta = xmax - x_list[-1]
        delta_x = 1/2*x_list[-1]

        Nx = delta/delta_x

        if Nx > Nx_max:
            while(True):
                delta_x *= 2
                Nx = delta/delta_x
                if Nx < Nx_max:
                    break
                delta_x = delta/Nx

        current_x = x_list[-1]
        while(True):
            
            if delta_x/current_x <= 0.2:
                delta_x *= 2
            current_x = x_list[-1] +  delta_x
            x_list.append( current_x )
            
            if current_x > xmax:
                break
            
    return x_list