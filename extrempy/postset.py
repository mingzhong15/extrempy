from extrempy.constant import *
from extrempy.dataset import SetSys

from dscribe.descriptors import SOAP
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


class PostSetSys(SetSys):

    def __init__(self, set_dir, is_printf=True):
        super().__init__(set_dir, is_printf)

    def _read_force(self,):
        super()._read_force()

    def _read_thermo(self, ):
        super()._read_thermo()

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
