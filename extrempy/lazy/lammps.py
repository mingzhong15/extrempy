from extrempy.lazy.lib import _get_mass_map

def _get_lattice_str(lattice_type, lattice_param):

    if lattice_type == 'fcc':
        strs = 'lattice             fcc %.6f \n'%(lattice_param)
        
    elif lattice_type == 'bcc':
        strs = 'lattice             bcc %.6f \n'%(lattice_param)
    elif lattice_type == 'sc':
        strs = 'lattice             sc %.6f \n'%(lattice_param)
    elif lattice_type == 'hcp':
        strs = 'lattice             hcp %.6f %.6f \n'%(lattice_param[0], lattice_param[1])
    elif lattice_type == 'diamond':
        strs = 'lattice             diamond %.6f \n'%(lattice_param)
    
    elif lattice_type == 'Mg2SiO4-fo1':

        strs =  'variable              a equal 4.75430    \n'
        strs += 'variable              b equal 10.20100   \n'
        strs += 'variable              c equal 5.98190    \n\n'

        strs += 'lattice               custom        1.0                                &  \n'
        strs += '                      a1            $a         0.0        0.0          &  \n'
        strs += '                      a2            0.0        $b         0.0          &  \n'
        strs += '                      a3            0.0        0.0        $c           &  \n'
        strs += '                      basis         0.00000    0.00000    0.00000      &  \n'   
        strs += '                      basis         0.50000    0.50000    0.00000      &  \n'  
        strs += '                      basis         0.00000    0.00000    0.50000      &  \n'
        strs += '                      basis         0.50000    0.50000    0.50000      &  \n'
        strs += '                      basis         0.98450    0.27420    0.25000      &  \n'
        strs += '                      basis         0.01550    0.72580    0.75000      &  \n'
        strs += '                      basis         0.48450    0.22580    0.75000      &  \n'
        strs += '                      basis         0.51550    0.77420    0.25000      &  \n'
        strs += '                      basis         0.42500    0.09640    0.25000      &  \n'
        strs += '                      basis         0.57500    0.90360    0.75000      &  \n'
        strs += '                      basis         0.92500    0.40360    0.75000      &  \n'
        strs += '                      basis         0.07500    0.59640    0.25000      &  \n'
        strs += '                      basis         0.76750    0.08990    0.25000      &  \n'
        strs += '                      basis         0.23250    0.91010    0.75000      &  \n'
        strs += '                      basis         0.26750    0.41010    0.75000      &  \n'
        strs += '                      basis         0.73250    0.58990    0.25000      &  \n'
        strs += '                      basis         0.22710    0.44070    0.25000      &  \n'
        strs += '                      basis         0.77290    0.55930    0.75000      &  \n'
        strs += '                      basis         0.72710    0.05930    0.75000      &  \n'
        strs += '                      basis         0.27290    0.94070    0.25000      &  \n'
        strs += '                      basis         0.26380    0.16920    0.02330      &  \n'
        strs += '                      basis         0.73620    0.83080    0.97670      &  \n'
        strs += '                      basis         0.76380    0.33080    0.97670      &  \n'
        strs += '                      basis         0.23620    0.66920    0.02330      &  \n'
        strs += '                      basis         0.73620    0.83080    0.52330      &  \n'
        strs += '                      basis         0.26380    0.16920    0.47670      &  \n'
        strs += '                      basis         0.23620    0.66920    0.47670      &  \n'
        strs += '                      basis         0.76380    0.33080    0.52330         \n\n'

        strs += 'region                box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice \n\n'
        strs += 'create_box            3 box \n'
        strs += 'create_atoms          3 box  & \n'
        strs += '                     basis   1   1   basis   2  1   basis   3  1   basis   4  1   & \n'
        strs += '                     basis   5   1   basis   6  1   basis   7  1   basis   8  1   & \n'
        strs += '                     basis   9   2   basis  10  2   basis  11  2   basis  12  2   & \n'
        strs += '                     basis  13   3   basis  14  3   basis  15  3   basis  16  3   & \n'
        strs += '                     basis  17   3   basis  18  3   basis  19  3   basis  20  3   & \n'
        strs += '                     basis  21   3   basis  22  3   basis  23  3   basis  24  3   & \n'
        strs += '                     basis  25   3   basis  26  3   basis  27  3   basis  28  3   \n\n'


    elif lattice_type == 'Mg2SiO4-fo2':

        strs =  'lattice         custom    1.0     & \n'
        strs += '                a1    4.6830000877      0.0000000000         0.0000000000  & \n'
        strs += '                a2   -1.2499409133      9.1247875820         0.0000000000  & \n'
        strs += '                a3   -1.5492150373     -0.4930850674         5.0623401655  & \n'
        strs += '                basis   0.9714  0.0966  0.6415  & \n'
        strs += '                basis   0.0419  0.9354  0.0567  & \n'
        strs += '                basis   0.0026  0.6656  0.7362  & \n'
        strs += '                basis   0.5065  0.516   0.8488  & \n'
        strs += '                basis   0.5386  0.7905  0.1388  & \n'
        strs += '                basis   0.0105  0.3666  0.9614  & \n'
        strs += '                basis   0.4746  0.2414  0.5588  & \n'
        strs += '                basis   0.5065  0.5161  0.3489  & \n'
        strs += '                basis   0.0048  0.3875  0.4564  & \n'
        strs += '                basis   0.0083  0.6445  0.2412  & \n'
        strs += '                basis   0.6136  0.1239  0.0951  & \n'
        strs += '                basis   0.3995  0.9081  0.6027  & \n'
        strs += '                basis   0.2702  0.843   0.8193  & \n'
        strs += '                basis   0.7865  0.6995  0.968   & \n'
        strs += '                basis   0.7628  0.2287  0.3717  & \n'
        strs += '                basis   0.2756  0.0474  0.4928  & \n'
        strs += '                basis   0.7429  0.1891  0.8786  & \n'
        strs += '                basis   0.2266  0.3325  0.7297  & \n'
        strs += '                basis   0.2506  0.8032  0.3262  & \n'
        strs += '                basis   0.7833  0.6947  0.4029  & \n'
        strs += '                basis   0.7924  0.4622  0.1735  & \n'
        strs += '                basis   0.7372  0.9845  0.2048  & \n'
        strs += '                basis   0.2204  0.5698  0.5243  & \n'
        strs += '                basis   0.776   0.4475  0.645   & \n'
        strs += '                basis   0.2685  0.1109  0.0031  & \n'
        strs += '                basis   0.237   0.5847  0.0528  & \n'
        strs += '                basis   0.23    0.3374  0.2947  & \n'
        strs += '                basis   0.7446  0.9211  0.695 \n\n'

        strs += 'region          box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice \n\n'

        strs += 'create_box      3 box \n'
        strs += 'create_atoms    3 box  & \n'
        strs += '                basis   1   1   basis   2  1   basis   3  1   basis   4  1   & \n'
        strs += '                basis   5   1   basis   6  1   basis   7  1   basis   8  1   & \n'
        strs += '                basis   9   2   basis  10  2   basis  11  2   basis  12  2   & \n'
        strs += '                basis  13   3   basis  14  3   basis  15  3   basis  16  3   & \n'
        strs += '                basis  17   3   basis  18  3   basis  19  3   basis  20  3   & \n'
        strs += '                basis  21   3   basis  22  3   basis  23  3   basis  24  3   & \n'
        strs += '                basis  25   3   basis  26  3   basis  27  3   basis  28  3   \n\n'

    elif lattice_type == 'Mg2SiO4-fo3':

        strs =  'variable        a equal 2.6400 \n'
        strs += 'variable        b equal 8.5960 \n'
        strs += 'variable        c equal 9.0400 \n\n'

        strs += 'lattice         custom    1.0     & \n'
        strs += '                a1      $a       0.0     0.0     & \n'
        strs += '                a2      0.0      $b      0.0     & \n'
        strs += '                a3      0.0      0.0     $c      & \n'
        strs += '                basis         0.50000    0.37300    0.35000  &  \n'
        strs += '                basis         0.50000    0.62700    0.85000  & \n'
        strs += '                basis         0.00000    0.87300    0.35000  & \n'
        strs += '                basis         0.00000    0.12700    0.85000  & \n'
        strs += '                basis         0.50000    0.88900    0.66200  &  \n'
        strs += '                basis         0.50000    0.11100    0.16200  &  \n'
        strs += '                basis         0.00000    0.38900    0.66200  &  \n'
        strs += '                basis         0.00000    0.61100    0.16200  &  \n'
        strs += '                basis         0.00000    0.85700    0.00000  &   \n'
        strs += '                basis         0.00000    0.14300    0.50000  &  \n'
        strs += '                basis         0.50000    0.35700    0.00000  &  \n'
        strs += '                basis         0.50000    0.64300    0.50000  &  \n'
        strs += '                basis         0.50000    0.98300    0.46000  &   \n'
        strs += '                basis         0.50000    0.01700    0.96000  &  \n'
        strs += '                basis         0.00000    0.48300    0.46000  &  \n'
        strs += '                basis         0.00000    0.51700    0.96000  &  \n'
        strs += '                basis         0.50000    0.23100    0.57200  &  \n'
        strs += '                basis         0.50000    0.76900    0.07200  &  \n'
        strs += '                basis         0.00000    0.73100    0.57200  &  \n'
        strs += '                basis         0.00000    0.26900    0.07200  &  \n'
        strs += '                basis         0.50000    0.70600    0.27100  &  \n'
        strs += '                basis         0.50000    0.29400    0.77100  &  \n'
        strs += '                basis         0.00000    0.20600    0.27100  &  \n'
        strs += '                basis         0.00000    0.79400    0.77100  &  \n'
        strs += '                basis         0.50000    0.46600    0.18100  & \n'
        strs += '                basis         0.50000    0.53400    0.68100  & \n'
        strs += '                basis         0.00000    0.96600    0.18100  &  \n'
        strs += '                basis         0.00000    0.03400    0.68100    \n\n'


        strs += 'region          box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice \n\n'

        strs += 'create_box      3 box \n'
        strs += 'create_atoms    3 box  & \n'
        strs += '                basis   1   1   basis   2  1   basis   3  1   basis   4  1   & \n'
        strs += '                basis   5   1   basis   6  1   basis   7  1   basis   8  1   & \n'
        strs += '                basis   9   2   basis  10  2   basis  11  2   basis  12  2   & \n'
        strs += '                basis  13   3   basis  14  3   basis  15  3   basis  16  3   & \n'
        strs += '                basis  17   3   basis  18  3   basis  19  3   basis  20  3   & \n'
        strs += '                basis  21   3   basis  22  3   basis  23  3   basis  24  3   & \n'
        strs += '                basis  25   3   basis  26  3   basis  27  3   basis  28  3   \n\n'

    else:
        raise ValueError('Lattice type %s not supported'%(lattice_type))
    
    return strs

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
