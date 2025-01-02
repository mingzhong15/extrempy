# N_toal -> numb_frames
# N_atom -> numb_atoms
# N_type -> numb_element

from extrempy.constant import *
from extrempy.md.base import MDSys

class SEDCalc(MDSys):

    def __init__(self, *args, SKIP=0, INTERVAL=1, **kwargs):
        super().__init__(*args, **kwargs)

        self.SKIP = SKIP
        self.INTERVAL = INTERVAL

        self.dump_idx = np.arange(self.SKIP, self.numb_frames, self.INTERVAL)
        self.numb_frames_used = self.dump_idx.shape[0]

        self.delta_k = 2 * np.pi / np.diag(self.cells[0])

    def _read_dump(self, tdx):
        return super()._read_dump(tdx)  
    
    # --------------------------------------------------- #
    # initialization of complex arrays (Real and Imaginary part)
    # --------------------------------------------------- #
    
    # 1 component 
    def _init_rkt(self):

        self.rkt_re = np.zeros((self.numb_frames_used, 1))
        self.rkt_im = np.zeros((self.numb_frames_used, 1))

    # 3 components : vx, vy, vz
    def _init_jkt(self):

        if self.numb_types == 1:

            self.jkt_re = np.zeros((self.numb_frames_used, 3))
            self.jkt_im = np.zeros((self.numb_frames_used, 3))

        elif self.numb_types > 1:
        
            self.jkt_re = np.zeros((self.numb_frames_used, 3 * ( self.numb_types + 1) ))
            self.jkt_im = np.zeros((self.numb_frames_used, 3 * ( self.numb_types + 1) ))    

    # --------------------------------------------------- #
    # calculation the correlation fuction for single atomic configuration
    # --------------------------------------------------- #
    def _calc_theta(self, k_vec):

        self.theta = np.dot(self.position, k_vec)

    # particle current correlation
    # j(k,t) = \sum_i v_i(t) * exp( i * k * r_i (t) )
    
    # real part : cos (k * r_i)
    # imaginary part : sin (k * r_i)

    def _calc_jkt(self, tdx):

        # jkt (all atoms)
        for idx in range(3):

            self.jkt_re[tdx, idx] = np.sum(  self.velocity[:,idx] * np.cos(self.theta) )
            self.jkt_im[tdx, idx] = np.sum( -self.velocity[:,idx] * np.sin(self.theta) )
        
        # jkt (each elements)
        if self.numb_types > 1:
  
            for edx in range(self.numb_types):
                
                cri = self.type == edx+1
                for idx in range(3):
                    self.jkt_re[tdx, idx+3*(edx+1)] = np.sum(  self.velocity[:,idx][cri] * np.cos(self.theta[cri]) )
                    self.jkt_im[tdx, idx+3*(edx+1)] = np.sum( -self.velocity[:,idx][cri] * np.sin(self.theta[cri]) )


    # calculate the correlation spectra for each component (jdx) of particle current 
    def _calc_ckw(self, dt, jdx):

        freq = np.fft.fftfreq(self.numb_frames_used, d=dt)

        tmp = np.transpose( np.vstack((self.jkt_re[:, jdx],self.jkt_im[:,jdx])) )
        j_kt = np.zeros(self.numb_frames_used, dtype=np.complex128)
        j_kt.real = tmp[:,0]
        j_kt.imag = tmp[:,1]

        j_kw = np.fft.fft(j_kt)
        j_kw_conj = np.fft.fft( j_kt.conjugate() )

        ckw = np.abs(  (j_kw+j_kw_conj)/2   )**2 / ( self.numb_frames_used * self.numb_atoms)

        out = np.transpose(np.vstack((freq, ckw)))
        return out[np.lexsort(out[:,::-1].T)]
        

    def _calc_sed_from_traj( self, save_dir, k_vec_tmp, nk = 1, suffix=None, is_save = True):

        self._init_jkt()

        unit_k_vec = k_vec_tmp / np.linalg.norm(k_vec_tmp)

        if suffix is None:
            suffix = '%.3f%.3f%.3f'%(unit_k_vec[0], unit_k_vec[1], unit_k_vec[2])

        k_vec = self.delta_k * unit_k_vec * nk

        for tdx in range( self.dump_idx.shape[0]):

            self._read_dump( self.dump_idx[tdx] * self.dump_freq)

            self._calc_theta(k_vec)
            self._calc_jkt(tdx)

        if is_save:
            np.savetxt( os.path.join(save_dir, 'jkt_%.dk_'%(nk)+suffix+'.txt'), np.hstack([self.jkt_re, self.jkt_im]) )

        labels = ['vx','vy','vz']

        tmp = np.zeros((self.numb_frames_used, 2))
        for jdx in range(3):

            ckw = self._calc_ckw(self.dt*self.INTERVAL*self.dump_freq, jdx)

            if jdx == 0:
                tmp = np.array(ckw)
            else:
                tmp[:,1] += ckw[:,1]

            if is_save:
                np.savetxt( os.path.join(save_dir, 'ckw_'+labels[jdx]+'_%.dk_'%(nk)+suffix+'.txt'), ckw, header='component-'+labels[jdx]+' nk=%.d frequency (1/fs) spectrum (arb. unit.)'%nk, comments='# ')

        if is_save:
            np.savetxt( os.path.join(save_dir, 'ckw_tot_%.dk_'%(nk)+suffix+'.txt'), tmp, header='nk=%.d frequency (1/fs) spectrum (arb. unit.)'%nk, comments='# ')
    
