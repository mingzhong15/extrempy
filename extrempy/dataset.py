from extrempy.constant import *

from dscribe.descriptors import SOAP
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

class SetSys():

    def __init__(self, set_dir, is_printf=True):
        
        self.SET_DIR = set_dir
        self.is_printf = is_printf

        if self.is_printf:
            print('read dataset from '+self.SET_DIR + ' ...')
        
        try:
            file_list = ['energy', 'force', 'virial', 'box']
            
            for file in file_list:
                os.path.exists( os.path.join(self.SET_DIR, 'set.000', file+'.npy') )

            if self.is_printf:
                print('basic dataset checked')
        except:

            print(file+ ' not found')
            
         # 通过type.raw获得原子数信息
        self.natoms = np.loadtxt(os.path.join(self.SET_DIR, 'type.raw')).shape[0]
        
        self.energy = np.load(os.path.join(self.SET_DIR,'set.000', 'energy.npy'))/self.natoms
        
        self.frames = self.energy.shape[0]

        if self.is_printf:
            print('%.d frames (%.d-atom-system) contained'%(self.frames, self.natoms))

        
    # 读取数据包含的力信息
    def _read_force(self,):

        force = np.load(os.path.join(self.SET_DIR,'set.000', 'force.npy')).reshape(-1, 3)

        # 计算每个原子受力的幅值，对单个frame的受力进行平均
        f_norm = np.linalg.norm(force, axis=1)
        self.f_ave = np.average( f_norm.reshape(-1, self.natoms), axis=1)
        self.f_std = np.std( f_norm.reshape(-1, self.natoms), axis=1)

        # 对不同温度的受力进行平均
        try:
            self.temp_list = np.unique(self.temp)
        except:
            self._read_thermo()
            
            self.temp_list = np.unique(self.temp)
        
        self.f_tave = []
        self.f_tstd = []
        for tt in self.temp_list:
            
            cri = self.temp == tt
            
            self.f_tave = np.append(self.f_tave, np.average(f_norm.reshape(-1, self.natoms)[cri]))
            self.f_tstd = np.append(self.f_tstd, np.std(f_norm.reshape(-1, self.natoms)[cri]))
            
    # 读取数据包含的热力学信息（温度、压强）    
    def _read_thermo(self, ):

        # 读取virial tensor的对角元，并做平均
        virials = np.load(os.path.join(self.SET_DIR, 'set.000', 'virial.npy'))
        alg = (virials[:,0] + virials[:,4] + virials[:,-1])/3

        # 读取box.npy来获得体积信息，单位是 Angstrom^3
        box = np.load(os.path.join(self.SET_DIR,'set.000',  'box.npy'))
        self.vol = np.linalg.det(box.reshape(-1,3,3))

        # 计算virial贡献的压强，并做单位转换，从eV/Ang^3转换为GPa
        self.pres =  alg/self.vol  /J2eV * (m2A**3) * Pa2GPa

        try:
            fpar = np.load(os.path.join(self.SET_DIR,'set.000', 'fparam.npy')).reshape(-1,)

            if self.is_printf:
                print('fparam is readed')

            # 读取fparam.npy获得温度信息，计算得到热动压，做单位转换
            therm_term = self.natoms/self.vol * kb * fpar * (m2A**3) *1e-9

            self.pres += therm_term
        except:

            print('fparam.npy is not found in ' + os.path.join(self.SET_DIR,'set.000', 'fparam.npy') )
            fpar = np.ones(self.energy.shape) * 0

        self.temp = fpar
        
        
    def _read_soap(self, INTERVAL = 10, target_species = None,  is_frame_average = True):

        print('soap descriptor is estimated at ', self.SET_DIR)
        
        confs = dpdata.LabeledSystem( self.SET_DIR, fmt="deepmd/npy")

        if target_species is None:
            self.species = list(np.unique(confs.get_atom_names()))
        else:
            self.species = target_species

        soap = SOAP(
                species=self.species,
                r_cut=6.0,            # 重要参数：截断半径
                periodic=True,        # 假设数据为周期性体系（根据实际情况调整）
                n_max=8,              # 径向基函数数量
                l_max=6,              # 角动量量子数
                sigma=1.0,            # 高斯宽度
                sparse=False,
                
            )
        
        all_descriptors = []
        # 遍历每个构型
        for idx in range(0, len(confs), INTERVAL):   
            # 转换为ASE Atoms对象（dpdata -> ASE转换）
            atoms = confs[idx].to_ase_structure()
            atom_symbols = atoms[0].get_chemical_symbols()
            mask = [symbol in self.species for symbol in atom_symbols]
            # 生成当前构型所有原子的SOAP描述符
            descriptors = soap.create(atoms[0][mask])

            if is_frame_average:
                # 对构型内所有原子的SOAP进行平均 (不推荐）
                descrip_ave = np.mean(descriptors, axis=0)
                all_descriptors.append(descrip_ave)

            else:
                all_descriptors.append(descriptors)


        if is_frame_average:
            X = np.concatenate(all_descriptors, axis=0)
            X = X.reshape(-1, descrip_ave.shape[0])      

            self.soap_frame = X
            print('vector of SOAP is generate : ', self.soap_frame.shape)
        else:
            X = np.concatenate(all_descriptors, axis=0)
            X = X.reshape(-1, descriptors.shape[1])

            self.soap_atomic = X

            print('vector of SOAP is generate : ', self.soap_atomic.shape)

def _get_case_list(SET_DIR, file='type.raw', index=-2):
    
    p_list = []
    tmp = glob.glob( os.path.join(SET_DIR, file))
    
    for tmp0 in tmp:
    
        p_list = np.append(p_list, tmp0.split('/')[index])

    return np.unique(p_list)

def _read_soap_from_pos_file( REF_DIR, target_species=None, fmt='vasp/poscar' ):

    label_ref = []
    all_soap_ref = []

    pos_list = glob.glob( os.path.join( REF_DIR, '*.POSCAR') )

    all_descriptors = []
    
    for idx, pos_file in enumerate(pos_list):

        print('read '+  pos_file )

        confs = dpdata.System(pos_file, fmt=fmt)
        if target_species is None:
            species = list(np.unique(confs.get_atom_names()))
        else:
            species = target_species
    
        soap = SOAP(
            species=species,
            r_cut=6.0,            # 重要参数：截断半径
            periodic=True,        # 假设数据为周期性体系（根据实际情况调整）
            n_max=8,              # 径向基函数数量
            l_max=6,              # 角动量量子数
            sigma=1.0,            # 高斯宽度
            sparse=False,
        )
        
        atoms = confs.to_ase_structure()
        atom_symbols = atoms[0].get_chemical_symbols()
        mask = [symbol in species for symbol in atom_symbols]
        descriptors = soap.create(atoms[0][mask])
        
        all_descriptors.append(descriptors) 

    X = np.concatenate(all_descriptors, axis=0)
    X = X.reshape(-1, descriptors.shape[1])

    return X

class MultiSetSys():

    def __init__(self, SET_DIR_LIST, is_printf=True):

        self.SET_DIR_LIST = SET_DIR_LIST
        self.is_printf = is_printf

        self.systems = []

        self._read_all_systems()
        
        # local atomic environment
        self.soap_a_all = None
        self.label_a_all = None
        self.soap_f_all = None
        self.label_f_all = None

        # thermodynamic information
        self.pres_a_all = None
        self.temp_a_all = None
        self.pres_f_all = None
        self.temp_f_all = None

        # force information
        self.f_fave_all = None
        self.f_fstd_all = None

    def _read_all_systems(self,):

        for SET_DIR in self.SET_DIR_LIST:

            sys = SetSys(SET_DIR, is_printf=self.is_printf)
            self.systems.append(sys)

            if self.is_printf:
                print('system '+ SET_DIR + ' is readed')

    def _read_all_soap(self, INTERVAL = 10, target_species = None,  is_frame_average = True):

        for s_idx, sys in enumerate(self.systems):

            sys._read_soap(INTERVAL = INTERVAL, target_species = target_species, is_frame_average = is_frame_average)

            if is_frame_average:
                sys.label_f = np.ones(sys.soap_frame.shape[0])*s_idx
            else:
                sys.label_a = np.ones(sys.soap_atomic.shape[0])*s_idx

        if is_frame_average:
            self.soap_f_all = np.concatenate([sys.soap_frame for sys in self.systems], axis=0)  
            self.label_f_all = np.concatenate([sys.label_f for sys in self.systems], axis=0)

        else:
            self.soap_a_all = np.concatenate([sys.soap_atomic for sys in self.systems], axis=0)
            self.label_a_all = np.concatenate([sys.label_a for sys in self.systems], axis=0)

    def _read_all_thermo(self, INTERVAL = 10, is_frame_average = True):

        for s_idx, sys in enumerate(self.systems):

            sys._read_thermo()

            # (frame, 1) -> (frame, natoms)
            if not is_frame_average:

                #print('SOAP shape: ', sys.soap_atomic.shape, 'pres shape: ', sys.pres.shape)
                natom = int(sys.soap_atomic.shape[0]/ sys.pres[::INTERVAL].shape[0])

                sys.pres_a = np.repeat(sys.pres[::INTERVAL], natom, axis=0)
                sys.temp_a = np.repeat(sys.temp[::INTERVAL], natom, axis=0)

        if is_frame_average:
            self.pres_f_all = np.concatenate([sys.pres[::INTERVAL] for sys in self.systems], axis=0)
            self.temp_f_all = np.concatenate([sys.temp[::INTERVAL] for sys in self.systems], axis=0)
        else:
            self.pres_a_all = np.concatenate([sys.pres_a for sys in self.systems], axis=0)
            self.temp_a_all = np.concatenate([sys.temp_a for sys in self.systems], axis=0)

    def _read_all_force(self, INTERVAL = 10, is_frame_average = True):

        if not is_frame_average:
            raise KeyError('force is not frame averaged')

        for s_idx, sys in enumerate(self.systems):

            sys._read_force()

        self.f_fave_all = np.concatenate([sys.f_ave[::INTERVAL] for sys in self.systems], axis=0)
        self.f_fstd_all = np.concatenate([sys.f_std[::INTERVAL] for sys in self.systems], axis=0)

