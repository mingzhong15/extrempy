from extrempy.constant import *

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
        
        
