from extrempy.constant import *
from extrempy.dataset import SetSys

class TestSys():
    
    def __init__(self, test_dir = './'):
        
        # path to DP prediction
        self.TEST_DIR = test_dir

        
        self.e_dft = []
        self.e_dp  = []
        
        self.f_dft_fave = []
        self.f_dft_fstd = []
        self.f_dp_fave  = []
        
        self.fc_rmse_fave = []
        self.fc_mae_fave  = []
        self.fcn_mae_fave = []
        self.fcn_rmse_fave = []
        
        self.fmag_mae_fave = []
        self.fmag_rmse_fave = []
        
        self.v_dft_fave = []
        self.v_dp_fave  = []
        self.v_rmse_fave = []
        self.v_mae_fave  = []
        
        self.p_dft  = []
        self.p_dp   = []
        
        self.vols = []
        self.temps = []
        self.press = []
        
    # 读取set里的温度、压强
    def _read_set(self, set_dir):
        
        self.set_sys = SetSys(set_dir)
        
        try:
            self.set_sys._read_thermo( )
        
            self.temps = np.append(self.temps, self.set_sys.temp)
            self.press = np.append(self.press, self.set_sys.pres)
            self.vols  = np.append(self.vols,  self.set_sys.vol)
        except:
            
            pass
        
    # 读取dp test的各类输出文件
    def _read_dp_test(self, prefix,  INTERVAL=1,  is_print=True):
        
        self.INTERVAL = INTERVAL
        self.prefix = prefix
        
        self._update_natoms(is_print)
        
        self._read_energy_test( is_print)
        self._read_force_test( is_print)
        self._read_virial_test( is_print)            
 
    def _update_natoms(self, is_print):
        
        energy = np.loadtxt( os.path.join(self.TEST_DIR, self.prefix+'.e.out') )
        force = np.loadtxt( os.path.join(self.TEST_DIR, self.prefix+'.f.out') )
        
        nframes = energy.shape[0]
        self.natoms = int(force.shape[0]/nframes)
        
        if is_print:
            print('Total atoms in this phase is %.d '%self.natoms)

    # 读取*.e.out，获得能量测试统计数据
    def _read_energy_test(self, is_print):
        energy = np.loadtxt( os.path.join(self.TEST_DIR, self.prefix+'.e.out') )/self.natoms
        
        self.e_dft = np.append(self.e_dft, energy[:,0])
        self.e_dp  = np.append(self.e_dp,  energy[:,1])
        
        # 计算能量的偏差值
        self.e_delta  = self.e_dft - self.e_dp
        
        # 对能量的RMSE/MAE进行所有frame的统计平均
        e_rmse = np.sqrt( np.average(self.e_delta**2) )
        e_mae = np.average( np.abs(self.e_delta) )
        
        if is_print:
            print('MAE of energy is %.5f meV/atom'%(e_mae*1e3))
            print('RMSE of energy is %.5f meV/atom'%(e_rmse*1e3))
            
    # (nframes, N*natom) -> (nframes, )
    # (nframes, N)       -> (nframes, )
    def frame_stats(self, data, is_atomic = True):
        
        lens = int(data.shape[1])
        
        if is_atomic:
            f_ave = np.average( data.reshape(-1, lens*self.natoms), axis=1)
            f_std = np.std( data.reshape(-1, lens*self.natoms), axis=1) 
            #如果除以数据量，是评估误差，不是样本误差 / np.sqrt( lens*self.natoms )
            
        else:
            f_ave = np.average( data.reshape(-1, lens), axis=1)
            f_std = np.std( data.reshape(-1, lens), axis=1) #/ np.sqrt( lens*self.natoms )
            
        return f_ave, f_std    
    
    # ==================================
    # f_rmse = sqrt ( average( (f_{x, DFT} - f_{x, DP}) ** 2 ) )
    # 读取*.f.out，获得受力测试统计数据
    def _read_force_test(self, is_print):

        force = np.loadtxt( os.path.join(self.TEST_DIR, self.prefix+'.f.out'))
        
        force_dft = force[:,:3]
        self.force_dp  = force[:,3:]
        force_delta = force_dft - self.force_dp
        
        # 0 力的幅值
        f_dft_fave, f_dft_fstd = self.frame_stats(np.linalg.norm(force_dft, axis=1).reshape(-1,1))
        f_dp_fave, f_dp_fstd = self.frame_stats(np.linalg.norm(self.force_dp, axis=1).reshape(-1,1))
        
        self.f_dft_fave = np.append(self.f_dft_fave, f_dft_fave)
        self.f_dft_fstd = np.append(self.f_dft_fstd, f_dft_fstd)
        self.f_dp_fave = np.append(self.f_dp_fave, f_dp_fave)

        
        # 计算力的统计平均值 （注意是分量的统计平均值，而不是力大小，f_dft-ave = 1/3 f_norm）
        #f_dft_fave, f_dft_fstd = np.sqrt(self.frame_stats(force_dft**2))
        #f_dp_fave, f_dp_fstd = np.sqrt(self.frame_stats(self.force_dp**2))
       
        

        # 1 力分量的统计偏差
        # 对每个atom的力进行平均，得到frame的平均RMSE/MAE 
        # （注意是分量的统计平均值，而不是力大小）
        fc_rmse_fave, fc_rmse_fstd = np.sqrt(self.frame_stats(force_delta**2))
        fc_mae_fave, fc_mae_fstd = self.frame_stats( np.abs(force_delta) )
        
        self.fc_rmse_fave = np.append(self.fc_rmse_fave, fc_rmse_fave)
        self.fc_mae_fave = np.append(self.fc_mae_fave, fc_mae_fave)

        # output statistics for whole dataset
        self.fc_rmse = np.sqrt(np.average(self.fc_rmse_fave**2))
        self.fc_mae = np.average(self.fc_mae_fave)
        if is_print:
            print('MAE of force component is %.5f eV/Angstrom'%self.fc_mae)
            print('RMSE of force component is %.5f eV/Angstrom'%self.fc_rmse)
            
        # 2 力偏差的幅值（考虑了受力方向偏差引入的误差）
        # force_delta是一个矢量，可以取模，转化为标量
        fc_norm = np.linalg.norm(force_delta, axis=1).reshape(-1, 1)
        
        fcn_mae_fave, fcn_mae_fstd = self.frame_stats( fc_norm )
        fcn_rmse_fave, fcn_rmse_fstd = np.sqrt(self.frame_stats( fc_norm**2 ))
        
        self.fcn_mae_fave = np.append(self.fcn_mae_fave, fcn_mae_fave)
        self.fcn_rmse_fave = np.append(self.fcn_rmse_fave, fcn_rmse_fave)
        
        # output statistics for whole dataset
        fcn_mae = np.average(self.fcn_mae_fave)
        fcn_rmse = np.average(self.fcn_rmse_fave)
        if is_print:
            print('MAE of the force error magnitude is %.5f eV/Angstrom'%fcn_mae)
            print('RMSE of the force error magnitude is %.5f eV/Angstrom'%fcn_rmse)
            
        # 3. 力幅值的偏差 （忽略了方向偏差带来的误差）
        mag_dft = np.linalg.norm(force_dft, axis=1).reshape(-1, 1)
        mag_dp  = np.linalg.norm(self.force_dp,  axis=1).reshape(-1, 1)
        
        fmag_mae_fave, fmag_mae_fstd = self.frame_stats( np.abs(mag_dft - mag_dp) )
        fmag_rmse_fave, fmag_rmse_fstd = np.sqrt(self.frame_stats( (mag_dft-mag_dp)**2))
        
        self.fmag_mae_fave = np.append(self.fmag_mae_fave, fmag_mae_fave)
        self.fmag_rmse_fave = np.append(self.fmag_rmse_fave, fmag_rmse_fave)
        
        # output statistics for whole dataset
        fmag_mae = np.average(self.fmag_mae_fave)
        fmag_rmse = np.average(self.fmag_rmse_fave)   
        if is_print:
            print('MAE of the force magnitude is %.5f eV/Angstrom'%fmag_mae)
            print('RMSE of the force magnitude is %.5f eV/Angstrom'%fmag_rmse)
        
        
    def _read_virial_test(self, is_print):
        virial = np.loadtxt(os.path.join(self.TEST_DIR, self.prefix+'.v.out')) 
        
        is_atomic = False
        
        v_dft = virial[:, :9]
        v_dp  = virial[:, 9:]
        v_delta = v_dft - v_dp

        # magnitude of virials for a single frame
        v_dft_fave, v_dft_fstd = np.sqrt(self.frame_stats(v_dft**2, is_atomic))
        v_dp_fave, v_dp_fstd = np.sqrt(self.frame_stats(v_dp**2, is_atomic))
        
        self.v_dft_fave = np.append(self.v_dft_fave, v_dft_fave)
        self.v_dp_fave = np.append(self.v_dp_fave, v_dp_fave)
        
        # rmse/mae of force for a single frame
        v_rmse_fave, v_rmse_fstd = np.sqrt(self.frame_stats(v_delta**2, is_atomic))
        v_mae_fave, v_mae_fstd = self.frame_stats( np.abs(v_delta), is_atomic )
        
        self.v_rmse_fave = np.append(self.v_rmse_fave, v_rmse_fave)
        self.v_mae_fave = np.append(self.v_mae_fave, v_mae_fave)

        # output statistics for whole dataset
        self.v_rmse = np.sqrt(np.average(self.v_rmse_fave**2))
        self.v_mae = np.average(self.v_mae_fave)
        if is_print:
            print('MAE of virial is %.5f eV'%self.v_mae)
            print('RMSE of virial is %.5f eV'%self.v_rmse)        


        # diagonal component of virial tensors [eV]
        alg_dft = (virial[:,0] + virial[:,4] + virial[:,8])/3 
        alg_dp  = (virial[:,9] + virial[:,13] + virial[:,17])/3
        
        try:
            # pressure [GPa]
            p_dft = alg_dft/(self.set_sys.vol)/J2eV * (m2A**3) *1e-9 
            p_dp  = alg_dp/(self.set_sys.vol)/J2eV * (m2A**3) *1e-9

            self.p_dft = np.append(self.p_dft, p_dft)
            self.p_dp = np.append(self.p_dp, p_dp)
            p_delta = p_dft - p_dp

            self.p_rmse = np.sqrt(np.average(p_delta**2))
            self.p_mae = np.average( np.abs(p_delta) )

            if is_print:
                print('MAE of pressure is %.5f GPa'%self.p_mae)
                print('RMSE of pressure is %.5f GPa'%self.p_rmse)
        except:  
            
            print('pressure esitimation failed, please check the path of SET_DIR')
            
            pass
            
    def _read_model_devi(self, MD_DIR, case ):
        
        data = np.loadtxt( os.path.join(MD_DIR, 'md.'+case+'.out'), comments='#')
        
        self.max_devi_f = data[:,4]
        self.min_devi_f = data[:,5]
        self.ave_devi_f = data[:,6]
        self.devi_e     = data[:,7]
        
        print('model deviation is estimated for %.d frames'%(self.devi_e.shape[0]))

        
# 多个DP模型用于评估模型偏差(model deviation)
class MultiTestSys():
    
    def __init__(self, test_dir, model_list):
        
        self.test_dir = test_dir
        self.model_list = model_list
        
        self.e_list = []
        self.f_list = []
        
        self.numb_model = len(model_list)
    
    def _read_model_devi(self, prefix):
        
        for model_dir in self.model_list:
            
            ss = TestSys( os.path.join(self.test_dir, model_dir) )
            ss._read_dp_test( prefix, is_print = False )
            ss._update_natoms( is_print = False)
            
            self.e_list = np.append(self.e_list, ss.e_dp)
            self.f_list = np.append(self.f_list, ss.force_dp)
        
            print(model_dir, ' is done \n')
            
        # energy matrix [n_model x n_frame x 1]
        self.e_list = self.e_list.reshape(self.numb_model, -1, 1)
        # force matrix [n_model x n_frame x n_atom x 3]
        self.f_list = self.f_list.reshape(self.numb_model, -1, ss.natoms, 3)
        
        self.devi_e = np.std(self.e_list, axis=0)
        self.devi_f = np.std(self.f_list, axis=0)
        
        # 对三个分量进行平均
        self.devi_f_ave = np.average( self.devi_f, axis = -1 )
        self.devi_f_sum = np.sum( self.devi_f, axis = -1 )
        # 对三个分量取模
        self.devi_f_norm= np.linalg.norm( self.devi_f, axis= -1)
        
        # frame-average -> np.average(data, axis=-1)
        
        #devi_f_ave = np.average( np.average(ms.devi_f, axis=-1), axis= -1)
        #devi_f_max = np.max( np.average(ms.devi_f, axis=-1), axis= -1)

