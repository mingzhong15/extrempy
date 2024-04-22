from extrempy.constant import *
import pandas as pd

class MDSys():
    
    def __init__(self, root_dir, traj_dir='traj', is_printf=True):
    
        self.root = root_dir
        self.traj_dir = os.path.join(self.root, traj_dir)
        
        dump_list = glob.glob( os.path.join( self.traj_dir, 'dump.*') )
        dump_list.sort()
        
        # =============================== 
        # read basic parameter
        # ===============================
        self.numb_frames = len(dump_list)
        self.numb_atoms = np.loadtxt(dump_list[0],skiprows=9)[:,0].shape[0]
        
        self.type_list = np.unique(np.loadtxt(dump_list[0],skiprows=9)[:,1])
        self.numb_elements = len(self.type_list)
        
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
    
            print('%.d types of atom in the frame'%self.numb_elements)
    
    # -------------------------------------- #
    # 
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
        
