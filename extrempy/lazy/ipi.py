import xml.etree.ElementTree as ET
from xml.dom import minidom
from extrempy.constant import *

class IPIInputGenerator:

    def __init__(self, work_dir, prefix='', latency_value='1e-3'):

        self.work_dir = work_dir
        self.prefix = prefix
        
        # 创建根元素
        self.simulation = ET.Element("simulation", verbosity="high")

        # 添加ffsocket元素
        ffsocket = ET.SubElement(self.simulation, "ffsocket", 
                                mode="unix", 
                                name="lammps", 
                                pbc="True")
        
        address = ET.SubElement(ffsocket, "address")
        address.text = " IPI "
        latency = ET.SubElement(ffsocket, "latency")
        latency.text = latency_value
        
        # 添加prng元素
        prng = ET.SubElement(self.simulation, "prng")
        seed = ET.SubElement(prng, "seed")
        seed.text = "%.d"%(np.random.randint(0, 10000))

    
    def _set_output(self, is_traj = True, traj_freq = 10,
                          is_thermo = True, thermo_freq = 100):
        
        # 创建output部分
        output = ET.SubElement(self.simulation, "output", prefix=self.prefix)

        if is_thermo:
            # 添加properties元素
            properties = ET.SubElement(output, "properties", 
                                      filename="thermo", 
                                      stride="%.d"%thermo_freq)
            properties.text = " [step,time{picosecond},temperature{kelvin},pressure_md{bar},pressure_cv{bar},kinetic_md{electronvolt},kinetic_cv{electronvolt},potential{electronvolt},conserved{electronvolt},density{g/cm3},volume{angstrom3},cell_h{angstrom}] "

        if is_traj:
                
            # 添加trajectory元素
            trajectory1 = ET.SubElement(output, "trajectory", 
                                       filename="xc", 
                                       stride="%.d"%traj_freq, 
                                       format="xyz", 
                                       cell_units="angstrom")
            trajectory1.text = " x_centroid{angstrom} "
            
            trajectory2 = ET.SubElement(output, "trajectory", 
                                       filename="vc", 
                                       stride="%.d"%traj_freq, 
                                       format="xyz", 
                                       cell_units="angstrom")
            trajectory2.text = " v_centroid{m/s} "
        
        # 添加checkpoint元素
        ET.SubElement(output, "checkpoint", 
                     filename="chk", 
                     stride="100", 
                     overwrite="True")

    def _set_configuration(self, nbeads = 32,
                          is_restart = False, input_file = 'data.xyz'):

        self.is_restart = is_restart
        self.nbeads = nbeads
        
        # 添加system元素
        self.system = ET.SubElement(self.simulation, "system")
        
        # 初始化部分
        initialize = ET.SubElement(self.system, "initialize", nbeads="%.d"%nbeads)

        if self.is_restart:
            file_elem = ET.SubElement(initialize, "file", mode="chk")
        else:
            file_elem = ET.SubElement(initialize, "file", mode="xyz")
            self.velocities = ET.SubElement(initialize, "velocities", mode="thermal", units="kelvin")
            
            #velocities.text = " %.2f "%T
            
        file_elem.text = input_file
        
        # 力场部分
        forces = ET.SubElement(self.system, "forces")
        force = ET.SubElement(forces, "force", forcefield="lammps")


    def _set_md(self, run_steps = 1000, dt = 0.5, 
                      ensemble = 'nvt',
                      T = 300, p = 1.0):  

        self.ensemble = ensemble
        
        # 添加total_steps元素
        total_steps = ET.SubElement(self.simulation, "total_steps")
        total_steps.text = " %.d "%run_steps

        if not self.is_restart:
            self.velocities.text = " %.2f "%T

        # 系综部分
        ensemble_str = ET.SubElement(self.system, "ensemble")
        temperature = ET.SubElement(ensemble_str, "temperature", units="kelvin")
        temperature.text = " %.2f "%T

        if 'p' in ensemble:
            pres = ET.SubElement(ensemble_str, "pressure", units="bar")
            pres.text = " %.6f "%(p*10000) # p in GPa units
    
        # 运动部分
        motion = ET.SubElement(self.system, "motion", mode="dynamics")
        fixcom = ET.SubElement(motion, "fixcom")
        fixcom.text = " True "
        
        dynamics = ET.SubElement(motion, "dynamics", mode=ensemble)
        timestep = ET.SubElement(dynamics, "timestep", units="femtosecond")
        timestep.text = " %.3f "%dt

        if 't' in ensemble:
            thermostat = ET.SubElement(dynamics, "thermostat", mode="pile_g")
            tau = ET.SubElement(thermostat, "tau", units="femtosecond")
            tau.text = " %.2f "%(dt*100)
            pile_lambda = ET.SubElement(thermostat, "pile_lambda")
            pile_lambda.text = " 0.2 "

        if 'p' in ensemble:
            barostat = ET.SubElement(dynamics, "barostat", mode="flexible")
            tau = ET.SubElement(barostat, "tau", units="femtosecond")
            tau.text = " %.2f "%(dt*1000)
        
    def _write_input(self, outfile='in.xml'):

        # 转换为XML字符串并美化输出
        rough_string = ET.tostring(self.simulation, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")

        # 保存到文件
        with open( os.path.join(self.work_dir, outfile), "w") as f:
            f.write(pretty_xml)


    def _generate_job_params(self, job_name, prj_id= None):

        self.job_name = job_name

        job_param = {
            "job_name": self.job_name,
            "command": "nohup i-pi in.xml &>log.ipi & mpirun --allow-run-as-root --use-hwthread-cpus -np %.d lmp -in input.lammps >log.run"%(self.nbeads),
            "log_file": "log.run",
            "backward_files": [self.prefix+"*", "log.*"],
            "project_id": prj_id,
            "platform": "ali",
            "machine_type": "1 * NVIDIA V100_16g",
            "job_type": "container",
            "image_address": "registry.dp.tech/dptech/dp/native/prod-14432/pimd-ipi:v1"
        }

        if self.nbeads <= 8:
            job_param['machine_type'] = 'c8_m32_1 * NVIDIA V100'

        elif self.nbeads <= 16:
            job_param['machine_type'] = 'c20_m76_2 * NVIDIA V100'
        
        elif self.nbeads <= 32:
            job_param['machine_type'] = 'c40_m152_4 * NVIDIA V100'

        else:
            raise ValueError('nbeads must be 8, 16, or 32')


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