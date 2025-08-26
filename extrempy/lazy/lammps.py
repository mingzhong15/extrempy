from extrempy.lazy.lib import _get_mass_map, _get_lattice_str
from extrempy.constant import *
import json

class LAMMPSInputGenerator:

    def __init__(self, work_dir, dt = 0.001, type_map = None, 
                 is_pimd = False, nbeads=0):

        self.work_dir = work_dir

        self.strs =  'units                 metal \n'
        self.strs += 'boundary              p p p \n'
        self.strs += 'atom_style            atomic \n'
        self.strs += 'atom_modify           map yes \n\n'

        self.dt = dt
        self.type_map = type_map
        self.mass_map = _get_mass_map(type_map)

        self.is_pimd = is_pimd
        self.nbeads  = nbeads
        self.is_ipi  = False
        self.is_rdf  = False
        self.is_thermo = False

        self.strs += 'variable              ibead uloop 99 pad\n'

        self.strs += 'variable              nbeads equal %.d\n\n'%self.nbeads

        self.strs += 'timestep              %.6f \n'%(self.dt)

    def _set_configuration(self, mode = 'from_file', 
                           data_file = None,
                           lattice_type = None, 
                           lattice_param = None,
                           scale = 1.0, 
                           numb_cell = [1,1,1]):

        self.strs += '\n # -------------------- configurations -------------------- #\n'


        self.data_file = data_file
        
        if mode == 'from_file':

            cmd = 'cp {} {}'.format(self.data_file, os.path.join(self.work_dir, 'init.data'))
            os.system(cmd)

            self.strs += 'read_data             {}\n'.format('init.data')

            self.strs += 'variable              LX equal lx \n'
            self.strs += 'variable              LY equal ly \n'
            self.strs += 'variable              LZ equal lz \n\n'

            self.scale = scale
            self.strs += 'variable              scale equal %.6f \n'%(self.scale)
            self.strs += 'variable              LXnew equal ${LX}*${scale} \n'
            self.strs += 'variable              LYnew equal ${LY}*${scale} \n'
            self.strs += 'variable              LZnew equal ${LZ}*${scale} \n\n'
            
            if scale != 1.0:
                self.strs += 'change_box                 all x final 0 ${LXnew} y final 0 ${LYnew} z final 0 ${LZnew} boundary p p p remap units box\n'

            self.strs += 'replicate             %.d %.d %.d \n\n'%(numb_cell[0], numb_cell[1], numb_cell[2])
            
        elif mode == 'lattice':

            self.strs += 'variable              Nx equal %.d \n'%(numb_cell[0])
            self.strs += 'variable              Ny equal %.d \n'%(numb_cell[1])
            self.strs += 'variable              Nz equal %.d \n\n'%(numb_cell[2])

            self.strs += _get_lattice_str(lattice_type, lattice_param)

            self.strs += 'variable              LX equal lx \n'
            self.strs += 'variable              LY equal ly \n'
            self.strs += 'variable              LZ equal lz \n\n'

        for idx in range(len(self.type_map)):
            self.strs += 'mass                  {} {} \n'.format(idx+1, self.mass_map[idx])
        self.strs += '\n'


    def _set_pes(self, mode ='dpmd', pes_file=None, fparam=None):
        self.strs += '\n # -------------------- PES -------------------- #\n'

        if mode == 'dpmd':

            cmd = 'cp {} {}'.format(pes_file, os.path.join(self.work_dir, 'cp.pb'))
            os.system(cmd)

            if fparam is None:
                self.strs += 'pair_style            deepmd cp.pb \n'
            else:
                self.strs += 'pair_style            deepmd cp.pb fparam {}\n'.format(fparam)
            self.strs += 'pair_coeff            * * \n'
            
        self.strs += 'neighbor              1.0 bin \n'
        self.strs += 'neigh_modify          every 10 delay 0 check no \n\n'

    def _set_thermo(self, thermo_freq = 100):

        self.strs += '\n # -------------------- thermo & dump -------------------- #\n'


        self.thermo_freq = thermo_freq

        self.strs += 'thermo_style          custom step temp press vol density etotal \n'
        self.strs += 'thermo_modify         format float %.4f \n'
        self.strs += 'thermo                {} \n\n'.format(thermo_freq)


        self.strs += 'variable              TEMP equal temp \n'
        self.strs += 'variable              PRESS equal press \n'
        self.strs += 'variable              ETOTAL equal etotal \n'
        self.strs += 'variable              VOL equal vol \n'
        self.strs += 'variable              RHO equal density \n'
        self.strs += 'variable              STEP equal step \n'
        self.strs += 'variable              ENTHALPY equal enthalpy \n\n'


    def _set_thermalization(self, ensemble = 'npt', 
                            run_steps = 10000, 
                            T = 300, p = 1.0,
                            T2 = None):

        self.strs += '\n # -------------------- thermalization -------------------- #\n'

        self.strs += 'variable              p equal %.6f \n'%(p*10000)
        self.strs += 'variable              T equal %.6f \n'%(T)
        if T2 is not None:
            self.strs += 'variable              T2 equal %.6f \n'%(T2)

        self.strs += 'variable              tdamp equal %.6f \n'%(100*self.dt)
        self.strs += 'variable              pdamp equal %.6f \n\n'%(1000*self.dt)

        self.strs += 'velocity              all create $T %.d dist gaussian rot yes \n'%(np.random.randint(0, 10000))

        if ensemble == 'npt':
            self.strs += 'fix                   1 all npt temp $T $T ${tdamp} aniso $p $p ${pdamp} \n'
            self.strs += 'run                   {} \n'.format(run_steps)
            self.strs += 'unfix                 1 \n\n'

        elif ensemble == 'heat-until-melt':
            self.strs += 'fix                   1 all npt temp $T $T ${tdamp} aniso $p $p ${pdamp} \n'
            self.strs += 'run                   {} \n'.format(10000)
            self.strs += 'unfix                 1 \n\n'

            self.strs += 'fix                   1 all npt temp ${T2} ${T2} ${tdamp} aniso $p $p ${pdamp} \n'
            self.strs += 'run                   {} \n'.format(run_steps)
            self.strs += 'unfix                 1 \n\n'            

            self.strs += 'fix                   1 all npt temp $T $T ${tdamp} aniso $p $p ${pdamp} \n'
            self.strs += 'run                   {} \n'.format(run_steps)
            self.strs += 'unfix                 1 \n\n'        

        #self.strs += 'reset_timestep 0 \n\n'

    def _set_ipi(self, run_steps = 1000, ):
        self.strs += '\n # -------------------- i-PI run -------------------- #\n'

        self.is_ipi  = True

        self.strs += 'fix                   1 all ipi IPI 2312 unix \n\n'

        self.strs += 'run                   {} \n\n'.format( (run_steps+2)*self.nbeads)
        

    def _set_dump(self, dump_freq = 100, dump_file = 'dump.lammps'):

        self.dump_file = dump_file

        if self.is_pimd:
            dump_file = dump_file + '${ibead}'
            
        self.strs += 'dump                 1 all custom %.d '%(dump_freq) + dump_file + ' id type xu yu zu vx vy vz \n'
        self.strs += 'dump_modify          1 sort id \n'

        tmp = ''
        for type0 in self.type_map:
            tmp += type0 + ' '
        self.strs += 'dump_modify           1 element ' + tmp + '\n\n'

    def _set_calculations(self, is_rdf = True, is_thermo = True):

        self.strs += '\n # -------------------- compute -------------------- #\n'


        self.is_rdf = is_rdf
        self.is_thermo = is_thermo

        if self.is_rdf:
            tmp = _generate_type_list(self.type_map)

            self.strs += 'compute               RDF all rdf 200 '+tmp+'\n'

            outfile = 'rdf.txt'
            if self.is_pimd:
                outfile = outfile + '${ibead}'
            
            self.strs += 'fix                   rdf all ave/time 100 100 10000 c_RDF[*] file '+outfile+' mode vector\n'

        if self.is_thermo:

            outfile = 'thermo.dat'
            if self.is_pimd:
                outfile = outfile + '${ibead}'

            self.strs += 'fix                   thermoprint all print %.d "${STEP} ${TEMP} ${PRESS} ${VOL} ${LX} ${LY} ${LZ} ${RHO} ${ETOTAL} ${ENTHALPY}" &\n'%(self.thermo_freq)
            self.strs += '                      title "# step temp[K] press[bars] vol[A^3] Lx[A] Ly[A] Lz[A] density[gcc] etotal[eV] enthalpy[eV]" &\n'
            self.strs += '                      file '+outfile+' screen no\n\n'

    def _set_trajectory(self, ensemble = 'npt', 
                        run_steps = 10000, T = 300, p = 1.0):

        self.strs += '\n # -------------------- run -------------------- #\n'


        if ensemble == 'npt':
            if not self.is_pimd: 
                self.strs += 'fix                   1 all npt temp $T $T ${tdamp} aniso $p $p ${pdamp} \n'

            else:
                self.strs += 'fix                   1 all pimd/langevin fmmode physical ensemble npt integrator obabo thermostat PILE_L ${ibead} temp $T tau ${tdamp} scale 1.0 barostat BZP aniso $p taup ${pdamp} \n'

        self.strs += 'run                   {} \n'.format(run_steps)
        self.strs += 'unfix                 1 \n'

    def _write_input(self, file_name = 'input.lammps'):


        if self.is_rdf:
            self.strs += 'unfix                 rdf \n'

        if self.is_thermo:
            self.strs += 'unfix                 thermoprint \n\n'

        if not self.is_ipi:
            self.strs += 'write_data            final.data \n\n'

        with open( os.path.join(self.work_dir, file_name), 'w') as f:
            f.write(self.strs)

    def _generate_job_params(self, job_name, platform= None, ncores=8):

        self.job_name = job_name
        if self.is_pimd:
            self.job_name = job_name + '_pimd'

        job_param = {
            "job_name": self.job_name,
            "command": "mkdir traj && lmp -in input.lammps > log.run",
            "log_file": "log.run",
            "backward_files": ["rdf.txt*", "thermo.dat*", "*log.lammps*", "final.data", self.dump_file+'*'],
            "project_id": platform,
            "platform": "ali",
            "machine_type": "1 * NVIDIA V100_16g",
            "job_type": "container",
            "image_address": "registry.dp.tech/dptech/dpmd:2.2.8-cuda12.0"
        }

        self.ncores = ncores
        
        if self.is_pimd:

            if self.ncores == 8:
                job_param['machine_type'] = 'c8_m32_1 * NVIDIA V100'

            elif self.ncores == 16:
                job_param['machine_type'] = 'c16_m62_1 * NVIDIA T4'
            
            elif self.ncores == 32:
                job_param['machine_type'] = 'c32_m64_cpu'

            else:
                raise ValueError('nbeads must be 8, 16, or 32')

            job_param['command'] = 'mpirun --use-hwthread-cpus --allow-run-as-root -np %.d lmp -in input.lammps -p %.dx%.d -log log'%(self.ncores, self.nbeads, int(self.ncores/self.nbeads))
    
        with open(os.path.join(self.work_dir, 'job.json'), 'w') as f:
            json.dump(job_param, f, indent=4)

    def _submit_job(self, platform='bh'):

        if platform == 'bh':
            os.chdir(self.work_dir)

            try:
                os.system('mkdir  ../'+self.job_name)
            except:
                pass

            pwd = 'lbg job submit -i job.json -p ./ -r ../'+self.job_name
            os.system(pwd)
        else:
            print('platform not supported')


def _generate_type_list(type_map):

    tmp = ''
    for idx in range(len(type_map)):

        for jdx in range(len(type_map)):

            if idx <= jdx:
                str = '%.d %.d'%(idx+1, jdx+1)
                tmp += str + ' '

    tmp = tmp.strip()
    return tmp
