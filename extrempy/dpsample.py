from extrempy.constant import *
from extrempy.dataset import SetSys

from monty.serialization import loadfn,dumpfn

class SampleSys():
    
    # ========================================= #
    # read basic information for DPGEN run
    # dpgen_dir : the working directory of DPGEN
    # ========================================= #
    def __init__(self, dpgen_dir, printf = True):
        
        self.dir = dpgen_dir
        self.param_file = os.path.join(dpgen_dir, 'param.json')
        self.jdata =loadfn(self.param_file)
        
        self.N_iter = len(glob.glob( os.path.join( dpgen_dir ,'iter.*' )))
        #self.N_iter = len(self.jdata['model_devi_jobs'])
        
        self.confs_list = self.jdata['sys_configs']
        self.label_list = []
        
        self.printf = printf
        
        print("DPGenerator System contains %.d Iterations"%(self.N_iter) )
        print("There are %.d initial configurations for exploration: \t "%(len(self.confs_list)) )
        for cc in self.confs_list:
            label = cc[0].split('.')[0]
            print(label,'\t\t')
            self.label_list.append(label)

    # =============================
    # part 1 DP training process check 
    # =============================
    def _get_loss(self, iter_idx, model_idx=0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        model_idx = '%.3d'%model_idx
    
        path = os.path.join(self.dir,iter_idx, '00.train', model_idx, 'lcurve.out')

        return np.loadtxt(path)
      
    def _plot_loss(self, ax, iter_idx, model_idx=0):
        #fig, ax = plt.subplots(1,3, figsize=(8,2),dpi=200)

        data = self._get_loss(iter_idx, model_idx)

        ll_list = ['Energy (eV)','Force (eV/$\\rm{\mathring A}$)','Virial (eV)']
        for idx in range(3):
            if idx == 0:
                ax[idx].plot(data[:,0],data[:,2+idx],'o',ms=1,mew=0.5, label='Iter.%.3d (DP%.3d)'%(iter_idx, model_idx))
            else:
                ax[idx].plot(data[:,0],data[:,2+idx],'o',ms=1,mew=0.5)

            ax[idx].set_xscale('log')
            ax[idx].set_yscale('log')
            ax[idx].set_title(ll_list[idx])
            ax[idx].set_xlabel('stopbatch')
          
    # =============================
    # part 2 model deviation during DPMD check 
    # ============================= 
  
    def _get_model_devi(self, iter_idx, sys_idx = 0, case_idx = 0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        task_idx = 'task.%.3d'%sys_idx+'.%.6d'%case_idx  
        
        path = os.path.join(self.dir,iter_idx, '01.model_devi', task_idx, 'model_devi.out')

        return np.loadtxt(path)
     
    def _get_thermo_from_md(self, iter_idx, sys_idx = 0, case_idx = 0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        task_idx = 'task.%.3d'%sys_idx+'.%.6d'%case_idx  
        
        path = os.path.join(self.dir,iter_idx, '01.model_devi', task_idx, 'input.lammps')

        temp = float(os.popen("grep ' TEMP ' " + path).readlines()[0].split()[-1])
        ele_temp = float(os.popen("grep ' ELE_TEMP ' " + path).readlines()[0].split()[-1])
        press = float(os.popen("grep ' PRES ' " + path).readlines()[0].split()[-1]) * bar2Pa /1e9
      
        return temp, ele_temp, press
      
    def _plot_model_devi(self, ax, iter_idx, sys_idx = 0, case_idx = 0, show_ele_temp = False):
        
        lo= self.jdata['model_devi_f_trust_lo']
        hi= self.jdata['model_devi_f_trust_hi']

        data = self._get_model_devi(iter_idx, sys_idx , case_idx)
        
        ax.plot(data[:,0],data[:,4],'o', mfc='none',mew=0.5,ms=2)
        
        ax.axhline(lo,linestyle='--',lw=0.5,color='k')
        ax.axhline(hi,linestyle='--',lw=0.5,color='k')  
        
        try:
            t,te,p = self._get_thermo_from_md(iter_idx, sys_idx , case_idx)

            output = 'T_i = %.d K'%t
            if show_ele_temp:
                output += 'T_e = %.d K'%te
            output += ' p = %.2f GPa, '%p
            output += self.label_list[sys_idx]
            ax.set_title(output, fontsize=7)
        except:
            pass

        ax.set_ylabel('model_devi ($\\rm{eV/\mathring A}$)')
        ax.set_xlabel('timestep')
        ax.set_xlim(data[0,0],data[-1,0])

          
    # =============================
    # part 3 configuration sampled check 
    # ============================= 

    def _get_thermo_from_data(self, iter_idx, sys_idx = 0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        data_idx = 'data.%.3d'%sys_idx
        
        path = os.path.join(self.dir, iter_idx, '02.fp', data_idx)
        
        # return (temps, press, vol, energy)
        ss = SetSys(path, is_printf=self.printf)
        ss._read_thermo()

        return ss.temp, ss.pres, ss.vol, ss.natoms
    
    def _plot_sampling(self, ax, iter_idx,  sys_idx = 0, color='dodgerblue', is_label=False, label=''):
        
        self.temps, self.press, self.vol, self.natoms = self._get_thermo_from_data(iter_idx, sys_idx)

        #print(self.temps)
        
        if is_label:
            ax.plot(self.press, self.temps, 
                'o', color=color, alpha=0.4, mew=0.5, mfc='none',ms=3, label=label)   
        else:
            ax.plot(self.press, self.temps, 
                'o', color=color, alpha=0.4, mew=0.5, mfc='none',ms=3)   
        
        if self.printf:
            output = 'Iter %.6d '%iter_idx
            output += ' add %.d frames '%self.temps.shape[0]
            output += '-- sys.%.3d '%sys_idx 
            output += '(' + self.label_list[sys_idx] +')'

            print(output)

        ax.set_xlabel('$p$ (GPa)')
        ax.set_ylabel('$T$ (K)')


    def _plot_all_sampling(self, ax, color):

        count = np.zeros(len(self.confs_list))
        
        self.frame_sys = np.zeros(len(self.confs_list))
        
        for i in range(self.N_iter):
            
            for sys_idx in range(len(self.confs_list)):
                
                try:
                  #label only plots in the first time for single system
                    if count[sys_idx] == 0:
                        self._plot_sampling(ax, iter_idx = i, sys_idx = sys_idx, 
                                        color=color[sys_idx], is_label=True, 
                                            label='DPGEN ('+self.label_list[sys_idx]+')')
                        count[sys_idx] +=1 
                    else:
                        self._plot_sampling(ax, iter_idx = i, sys_idx = sys_idx, 
                                        color=color[sys_idx])
                
                    self.frame_sys[sys_idx] += self.temps.shape[0]
                
                    if self.printf:
                        output = 'Iter %.6d --- sys. %.3d'%(i, sys_idx)
                        output += ' %.d frames in total '%self.frame_sys[sys_idx]
                        
                        print(output)
                
                except:

                    
                    pass        

        
    # in case of electron-temperature sampling condition
    def _obtain_tt_te_sampling(self):
        
        input_list = glob.glob( os.path.join(self.dir, 'iter.*', '01.model_devi', 'task*', 'input.lammps') )

        tt = []
        te = []

        for path in input_list:
            temp = float(os.popen("grep ' TEMP ' " + path).readlines()[0].split()[-1])
            ele_temp = float(os.popen("grep ' ELE_TEMP ' " + path).readlines()[0].split()[-1])
            #press = float(os.popen("grep ' PRES ' " + path).readlines()[0].split()[-1]) * bar2Pa /1e9

            tt = np.append(tt, temp)
            te = np.append(te, ele_temp)
            
        return tt,te          
      
    # =============================
    # part 4 extract internal energy & entropy  
    # ============================= 
  
    def _get_extra_raw(self, iter_idx, sys_idx = 0, str_s='energy  without entropy', save_raw = 'internal.raw'):
        
        iter_idx = 'iter.%.6d'%iter_idx

        path = os.path.join(self.dir, iter_idx, '02.fp', 'task.%.3d.*'%sys_idx, 'OUTCAR')

        fp_list = glob.glob( path )

        tmp = []
        
        for task_idx in range(len(fp_list)):
            
            outcar = os.path.join(self.dir, iter_idx, '02.fp', 'task.%.3d.%.6d'%(sys_idx, task_idx), 'OUTCAR')
        
            #print(outcar)
        
            cmd = "grep '"+str_s+"' " + outcar
           
            tmp = np.append(tmp, float(os.popen(cmd).readlines()[0].split()[3]) )
            
        save_dir = os.path.join(self.dir, iter_idx, '02.fp', 'data.%.3d'%sys_idx)
        np.savetxt( os.path.join(save_dir, save_raw), tmp.reshape(-1,1) )
    
        print(os.path.join(save_dir, save_raw),' is generated')

        
    # =============================
    # part 5 collect all configurations (including fparam.raw)  
    # =============================                 

    # 注意，该命令是提取已有的data.00*，而不是从OUTCAR重新提取
    def _collect_data(self, out_dir, set_numb = 20000, prefix='', exe_path = '~/raw_to_set.sh'):
        
        for sys_idx in range(len(self.label_list)):

            file_list = glob.glob(os.path.join(self.dir, 'iter.*', '02.fp', prefix+'data.%.3d'%sys_idx))
            file_list.sort()
            
            print('Sys.%.3d is working, there are %.d sub-datasets '%(sys_idx, len(file_list)))
            
            if len(file_list) != 0:

                fparam = []
                internal = []
                ms = dpdata.MultiSystems()

                for file in file_list:

                    sys = dpdata.LabeledSystem( file, fmt='deepmd/raw' )
                    ms.append(sys)

                    try:
                        tmp = np.loadtxt( os.path.join(file,'fparam.raw') )
                        fparam = np.append(fparam, tmp)
                        
                        is_fparam = True
                        
                    except:
                        is_fparam = False
                        
                    try:
                        tmp    = np.loadtxt( os.path.join(file,'internal.raw') )
                        internal = np.append(internal, tmp)
                        
                        is_internal = True
                    except:
                        is_internal = False                        

                outdir_ss = os.path.join(out_dir, self.label_list[sys_idx])

                try:
                    print('create output directories : '+outdir_ss)
                    os.mkdir( outdir_ss)
                except:
                    pass

                natom = sys.get_natoms()
                aparam = np.expand_dims(fparam, 0).repeat(natom,axis=0).T
                    
                ms.to_deepmd_raw(outdir_ss)
                
                if is_fparam:
                    np.savetxt( os.path.join(outdir_ss, 'fparam.raw'), fparam)
                    np.savetxt( os.path.join(outdir_ss, 'aparam.raw'), aparam)
                
                    print('NOTE: fparam.raw is generated')
                if is_internal:
                    np.savetxt( os.path.join(outdir_ss, 'internal.raw'), internal)
                    print('NOTE: internal.raw is generated')
                    
                os.chdir(outdir_ss)
                os.system('mv ./*/*.raw ./')
                
                os.system(exe_path + ' %.d'%(set_numb))

            print('Sys.%.3d is done'%(sys_idx))


import json


def _get_mass_map(type_map):

    mass_map_ref = {
        'H': 1.00794, 'He': 4.0026,

        'Li': 6.941,  'Be': 9.0122,  'B': 10.811,  'C': 12.0107, 
        'N':  14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797, 

        'Na': 22.9897, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 
        'P':  30.9738, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,

        'K': 39.0983, 'Ca': 40.078, 'Sc': 44.9559, 'Ti': 47.867, 
        'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938, 'Fe': 55.845,
        'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 
        'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216, 'Se': 78.96,
        'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 
        'Y': 88.9059, 'Zr': 91.224, 'Nb': 92.9064, 'Mo': 95.94,
        'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 
        'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.710, 
        'Sb': 121.760, 'Te': 127.60, 'I': 126.904, 'Xe': 131.293, 
        
        'Cs': 132.9055, 'Ba': 137.327, 'La': 138.9055, 'Ce': 140.116, 
        'Pr': 140.9077, 'Nd': 144.242, 'Pm': 145.00, 'Sm': 150.36, 
        'Eu': 151.964,  'Gd': 157.25, 'Tb': 158.9253, 'Dy': 162.50, 
        'Ho': 164.9303, 'Er': 167.259, 'Tm': 168.9342, 'Yb': 173.04,
        'Lu': 174.967,  'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 
        'Re': 186.207,  'Os': 190.23, 'Ir': 192.217,  'Pt': 195.084, 
        'Au': 196.9665, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 
        'Bi': 208.9804, 'Po': 209.0,  'At': 210.0, 'Rn': 222.0, 
        'Fr': 223.0,    'Ra': 226.0,  'Ac': 227.0, 'Th': 232.0381, 
        'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 
        'Am': 243.06138, 'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 
        'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0984, 'No': 259.1010, 
        'Lr': 262.1096, 'Rf': 261.1086, 'Db': 262.1138, 'Sg': 266.1218, 
        'Bh': 264.1254, 'Hs': 265.1388, 'Mt': 266.1516, 'Ds': 269.1649, 
        'Rg': 272.174, 'Cn': 277.1926, 'Nh': 278.1884, 'Fl': 281.1972, 
        'Mc': 282.1935, 'Lv': 285.2078, 'Ts': 286.2102, 'Og': 289.2184
    }

    mass_map = []
    for key in type_map:
        mass_map.append(mass_map_ref[key])

    return mass_map


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


    def _set_init_data(self, SET_DIR, SET_LIST):

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