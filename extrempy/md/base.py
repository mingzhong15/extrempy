from extrempy.constant import *
import pandas as pd

import re

def read_and_sort_files(directory, pattern='dump.*'):
    # 使用glob.glob读取文件列表
    file_list = glob.glob(os.path.join(directory, pattern))
    
    # 提取文件名中的数字并排序
    file_list.sort(key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
    
    return file_list

def calculate_dump_freq(file_list):
    # 提取文件名中的数字
    numbers = [int(re.findall(r'\d+', os.path.basename(f))[0]) for f in file_list]
    
    # 计算相邻数字的差值
    intervals = [j - i for i, j in zip(numbers[:-1], numbers[1:])]
    
    # 返回最小间隔
    return min(intervals)

class MDSys():
    
    def __init__(self, root_dir, dt = 1.0, traj_dir='traj', is_printf=True):
    
        self.root = root_dir
        self.traj_dir = os.path.join(self.root, traj_dir)

        dump_list = read_and_sort_files( self.traj_dir, pattern='dump.*')
        self.dump_freq = calculate_dump_freq(dump_list)

        self.dt = dt
        
        # =============================== 
        # read basic parameter
        # ===============================
        self.numb_frames = len(dump_list)
        self.numb_atoms = np.loadtxt(dump_list[0],skiprows=9)[:,0].shape[0]
        
        self.type_list = np.unique(np.loadtxt(dump_list[0],skiprows=9)[:,1])
        self.numb_types = len(self.type_list)
        
        keys = np.array(pd.read_csv(dump_list[0]).values[7][0].split()[2:])

        if 'type' in keys:
            self.type_idx = np.where(keys=='type')[0][0]
            
        if 'xu' in keys:
            self.x_idx = np.where(keys=='xu')[0][0]
            print('----xu, yu, zu contained----')
            
        if 'vx' in keys:
            self.vx_idx = np.where(keys=='vx')[0][0]
            print('----vx, vy, vz contained----')
            
        # usually valid for isochoric ensemble
        tmp = dpdata.System(dump_list[0], fmt='lammps/dump')
        self.cells = tmp.data['cells']
        
        # output basic parameters
        
        if is_printf:        
            print('%.d dump files in total'%self.numb_frames)
            print('%.d atoms in single frame'%self.numb_atoms )
    
            print('%.d types of atom in the frame'%self.numb_types)
    
    # -------------------------------------- #
    # read dump file
    # -------------------------------------- #
    def _read_dump(self, idx):
        
        tmp = np.loadtxt( os.path.join( self.traj_dir, 'dump.%.d'%idx), skiprows=9)

        self.type = tmp[:,self.type_idx]        

        try:
            self.position = tmp[:,self.x_idx:self.x_idx+3]
        except:
            pass

        try:
            self.velocity = tmp[:,self.vx_idx:self.vx_idx+3]
        except:
            pass
        
