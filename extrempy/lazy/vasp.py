import numpy as np
import os
import dpdata
import json

from extrempy.constant import kb, J2eV
from .base import InputGenerator,fermi_dirac

class VASPGenerator(InputGenerator):

    def __init__(self, *arg, potcar_lib_path,
                poscar_file = None, **kwargs):

        super().__init__(*arg, **kwargs)

        self.potcar_lib_path = potcar_lib_path

        self.poscar_file = poscar_file
        self._generate_poscar()
        
        self._generate_potcar()


    def _generate_poscar(self):

        try:
            sys = dpdata.System(self.poscar_file, fmt='vasp/poscar')
            self.numb_atom = sys.get_atom_numbs()
            self.atom_names = sys.get_atom_names()

            os.system('cp '+self.poscar_file+' '+os.path.join(self.work_path, 'POSCAR'))

            for i in range(len(self.atom_names)):
                print('POSCAR contains %d '%(self.numb_atom[i]) + self.atom_names[i] +' atoms')

        except:
            print('POSCAR file is not found')
            pass
        
    def _generate_potcar(self):

        cmd = 'cat '

        for element in self.atom_names:

            potcar_file = os.path.join(self.potcar_lib_path, element, 'POTCAR')

            suffix = ['_pv', '_sv']

            count = 0
            while not os.path.exists(potcar_file):
                

                potcar_file = os.path.join(self.potcar_lib_path, element+suffix[count], 'POTCAR')

                count += 1

                if count >= len(suffix):
                    break

            cmd += potcar_file

            with open(potcar_file, 'r') as f:
                f.readline()
                val = f.readline().split()
                self.numb_valence = int(float(val[0]))

        cmd += ' >> ' +os.path.join(self.work_path, 'POTCAR.tmp') + ' \n'

        cmd += 'mv '+os.path.join(self.work_path, 'POTCAR.tmp') + ' '+os.path.join(self.work_path, 'POTCAR') + ' \n'

        os.system(cmd)


    def set_params(self, encut=600, ele_temp=300, nbands=None, scale = 1.0, nband_min=5):

        self.encut = encut
        self.ele_temp = ele_temp

        nbands_0 = int(self.numb_valence * self.numb_atom[0] / 2 * scale) + nband_min * self.numb_atom[0]

        if nbands is None:
            self.nbands = nbands_0
        else:
            self.nbands = nbands

    def generate_incar(self, md_steps=100, dt = 1, 
                             latt_temp=None, mode='scf'):

        # basic settings
        if latt_temp is None:
            latt_temp = self.ele_temp

        cmd = '# CONTROL \n'
        cmd += 'ISTART = 0 \n'
        cmd += 'ICHARG = 2 \n'
        cmd += 'LWAVE = .FALSE. \n'
        cmd += 'LCHARG = .FALSE. \n'
        cmd += ' \n'

        cmd += '# ELECTRON \n'
        cmd += 'ENCUT = %.d \n'%(self.encut)
        cmd += 'NELM = 100 \n'
        cmd += 'ALGO = Normal \n'
        cmd += 'PREC = High \n'
        cmd += 'ISMEAR = -1 \n'
        cmd += 'SIGMA = %.16f \n'%(self.ele_temp * kb * J2eV)
        cmd += 'EDIFF = 1E-6 \n'
        cmd += 'NBANDS = %.d \n'%(self.nbands)
        cmd += ' \n'

        cmd += '# XC FUNCTIONAL \n'
        cmd += 'GGA = PS \n'
        cmd += ' \n'

        if mode == 'scf':
            cmd += '# ION \n'
            cmd += 'IBRION = -1 \n'
            cmd += 'NSW = 0 \n'
        elif mode == 'md':
            cmd += '# ION \n'
            cmd += 'IBRION = 0 \n'
            cmd += 'ISIF = 2 \n'
            cmd += 'NSW = %.d \n'%(md_steps)
            cmd += 'POTIM = %.2f \n'%(dt)
            cmd += 'TEBEG = %.d \n'%(latt_temp)
            cmd += 'TEEND = %.d \n'%(latt_temp)
            cmd += 'SMASS = 0 \n'
            cmd += 'MDALGO = 2 \n'
            cmd += ' \n'

        cmd += ' \n'
        cmd += '# PARALLIZATION \n'
        cmd += 'LREAL = Auto \n'
        cmd += 'KPAR = 4 \n'
        cmd += 'NPAR = 4 \n'
        cmd += ' \n'
        cmd += '# K-POINTS \n'
        cmd += 'KSPACING = 0.5 \n'
        cmd += 'KGAMMA = .TRUE. \n'
        cmd += ' \n'
        cmd += ' \n'

        with open(os.path.join(self.work_path, 'INCAR'), 'w') as f:
            f.write(cmd)


class VASPParamGenerator:

    def __init__(self, work_path, poscar_path, potcar_path, job_path):
        
        self._read_poscar(poscar_path)
        self._read_potcar(potcar_path)

        self.work_path = work_path
        self.poscar_path = poscar_path
        self.potcar_path = potcar_path
        self.job_path = job_path
        

        for i in range(len(self.atom_names)):

            print('POSCAR contains %d '%(self.numb_atom[i]) + self.atom_names[i] +' atoms')

        print('POTCAR contains %d valence electrons'%(self.numb_valence))

        
    def _read_poscar(self, poscar_path):

        sys = dpdata.System(poscar_path, fmt='vasp/poscar')
        self.numb_atom = sys.get_atom_numbs()
        self.atom_names = sys.get_atom_names()

    def _read_potcar(self, potcar_path):

        with open(potcar_path, 'r') as f:

            f.readline()

            val = f.readline().split()

            self.numb_valence = int(float(val[0]))

    def _generate_job_params(self, job_name):

        with open(self.job_path, 'r') as f:
            job_param = json.load(f)

        self.job_name = job_name
        job_param["job_name"] = self.job_name

        with open(os.path.join(self.work_path, 'job.json'), 'w') as f:
            json.dump(job_param, f, indent=4)


    def _set_params(self, encut=600, ele_temp=300, nbands=None, nband_min=5):

        self.encut = encut
        self.ele_temp = ele_temp

        nbands_0 = int(self.numb_valence * self.numb_atom[0] / 2) + nband_min * self.numb_atom[0]

        if nbands is None:
            self.nbands = nbands_0
        else:
            self.nbands = nbands

    def _generate_param_input(self, md_steps=100, dt = 1, 
                              latt_temp=None, mode='scf'):

        # basic settings
        if latt_temp is None:
            latt_temp = self.ele_temp

        cmd = '# CONTROL \n'
        cmd += 'ISTART = 0 \n'
        cmd += 'ICHARG = 2 \n'
        cmd += 'LWAVE = .FALSE. \n'
        cmd += 'LCHARG = .FALSE. \n'
        cmd += ' \n'

        cmd += '# ELECTRON \n'
        cmd += 'ENCUT = %.d \n'%(self.encut)
        cmd += 'NELM = 100 \n'
        cmd += 'ALGO = Normal \n'
        cmd += 'PREC = High \n'
        cmd += 'ISMEAR = -1 \n'
        cmd += 'SIGMA = %.16f \n'%(self.ele_temp * kb * J2eV)
        cmd += 'EDIFF = 1E-6 \n'
        cmd += 'NBANDS = %.d \n'%(self.nbands)
        cmd += ' \n'

        if mode == 'scf':
            cmd += '# ION \n'
            cmd += 'IBRION = -1 \n'
            cmd += 'NSW = 0 \n'
        elif mode == 'md':
            cmd += '# ION \n'
            cmd += 'IBRION = 0 \n'
            cmd += 'ISIF = 2 \n'
            cmd += 'NSW = %.d \n'%(md_steps)
            cmd += 'POTIM = %.2f \n'%(dt)
            cmd += 'TEBEG = %.d \n'%(latt_temp)
            cmd += 'TEEND = %.d \n'%(latt_temp)
            cmd += 'SMASS = 0 \n'
            cmd += 'MDALGO = 2 \n'
            cmd += ' \n'

        cmd += ' \n'
        cmd += '# PARALLIZATION \n'
        cmd += 'LREAL = Auto \n'
        cmd += 'KPAR = 4 \n'
        cmd += 'NPAR = 4 \n'
        cmd += ' \n'
        cmd += '# K-POINTS \n'
        cmd += 'KSPACING = 0.5 \n'
        cmd += 'KGAMMA = .TRUE. \n'
        cmd += ' \n'
        cmd += ' \n'

        with open(os.path.join(self.work_path, 'INCAR'), 'w') as f:
            f.write(cmd)

        os.system('cp '+self.poscar_path+' '+os.path.join(self.work_path, 'POSCAR'))
        os.system('cp '+self.potcar_path+' '+os.path.join(self.work_path, 'POTCAR'))


    def _submit_jobs(self):

        os.chdir(self.work_path)

        try:
            os.system('mkdir  ../'+self.job_name)
        except:
            pass

        pwd = 'lbg job submit -i job.json -p ./ -r ../'+self.job_name
        os.system(pwd)


class VASPReader:

    def __init__(self, work_path, is_prinft=False):
        self.work_path = work_path

        self.is_prinft = is_prinft

    def _read_outcar(self):

        with open(os.path.join(self.work_path, 'OUTCAR'), 'r') as f:
            
            count = 0

            while True:

                line = f.readline()

                if 'band No.  band energies     occupation' in line and count == 0:

                    band = []
                    band_energy = []
                    band_occupation = []

                    while True:

                        line = f.readline()

                        if len(line.split()) == 3:

                            tmp = line.split()

                            band.append(int(tmp[0]))
                            band_energy.append(float(tmp[1]))
                            band_occupation.append(float(tmp[2]))

                        else:
                            break

                        # if '---------' in line:
                        #     break
                    count += 1
                    self.band = np.vstack((band, band_energy, band_occupation)).T

                    if self.is_prinft:
                        print('%.d Band data read'%(ss.band.shape[0]))

                if 'SIGMA = ' in line:
                    self.sigma = float(line.split()[-1])
                    self.ele_temp = self.sigma/kb/J2eV

                    if self.is_prinft:
                        print('Fermi-Dirac smearing width: %.12f eV (equal to %.2f K)'%(self.sigma, self.ele_temp))

                if 'E-fermi' in line:
                    self.efermi = float(line.split()[2])
                    
                    if self.is_prinft:
                        print('Fermi energy: %.12f eV'%self.efermi)

                if 'free  energy' in line:

                    self.free_energy = float(line.split()[-2])
                    
                    if self.is_prinft:
                        print('Free energy: %.12f eV'%self.free_energy)

                if 'energy  without entropy' in line:

                    self.internal_energy = float(line.split()[-4])

                    if self.is_prinft:
                        print('Internal energy: %.12f eV'%self.internal_energy)

                if 'external pressure' in line:
                    self.press = float(line.split()[3])
                    
                    if self.is_prinft:
                        print('Pressure: %.12f kBar'%self.press)

                if 'Total CPU time used' in line:

                    self.cpu_time = float(line.split()[-1])

                    if self.is_prinft:
                        print('CPU time: %.12f s'%self.cpu_time)

                    break

    def _plot_band(self, ax, cc='dimgray',ele_temp=None):

        ax.plot(self.band[:,1], self.band[:,2]/2, 'o', ms=4, mew=0.5, color=cc,mfc='none' )

        ax.axvline(self.efermi, color=cc, ls=':', lw=1.0)

        xx = np.linspace(self.band[:,1].min(), self.band[:,1].max(), 1000)

        if ele_temp is None:
            ele_temp = self.ele_temp

        yy = fermi_dirac(xx, self.efermi, ele_temp)

        ax.plot(xx, yy, '-', color=cc, lw=1.0)

        y0 = np.zeros_like(xx)

        ax.fill_between(xx, y0, yy, facecolor=cc, alpha=0.1)

        ax.set_xlabel('$E$ (eV)')
        ax.set_ylabel('$g(E)$')

        ax.set_xlim(xx.min(), )

        cri = self.band[:,2] > 0.0

        self.min_nband = (self.band[:,1][cri]).shape[0]

        print( self.min_nband, 'bands are needed to be occupied with > 1e-5 fraction' )
