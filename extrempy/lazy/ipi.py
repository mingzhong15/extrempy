import xml.etree.ElementTree as ET
from xml.dom import minidom
from extrempy.constant import *
import json
import re
import time
import sys

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


def trans_data2ipixyz(work_dir, 
                       data_file, 
                       output_file, 
                       type_map=None):
    """
    将LAMMPS data文件转换为i-PI兼容的xyz格式
    
    参数:
    data_file: LAMMPS data文件路径
    output_file: 输出的xyz文件路径
    type_map: 原子类型到元素符号的映射字典(例如: {1: 'O', 2: 'H'})
    """
    # 默认原子类型映射
    if type_map is None:
        type_map = {1: 'O', 2: 'H'}  # 根据您的系统调整
    
    # 读取LAMMPS data文件
    with open( os.path.join(work_dir, data_file), 'r') as f:
        lines = f.readlines()
    
    # 解析基本参数
    natoms = 0
    xlo, xhi, ylo, yhi, zlo, zhi = 0, 0, 0, 0, 0, 0
    xy, xz, yz = 0, 0, 0
    atoms_section = False
    masses_section = False
    atoms = []
    
    for line in lines:
        if 'atoms' in line and not natoms:
            natoms = int(line.split()[0])
        elif 'xlo xhi' in line:
            xlo, xhi = map(float, line.split()[:2])
        elif 'ylo yhi' in line:
            ylo, yhi = map(float, line.split()[:2])
        elif 'zlo zhi' in line:
            zlo, zhi = map(float, line.split()[:2])
        elif 'xy xz yz' in line:
            xy, xz, yz = map(float, line.split()[:3])
        elif 'Masses' in line:
            masses_section = True
            atoms_section = False
            continue
        elif 'Atoms' in line:
            atoms_section = True
            masses_section = False
            continue
        elif masses_section:
            # 跳过质量部分
            continue
        elif atoms_section and line.strip():
            # 解析原子行
            parts = line.split()
            if len(parts) >= 5:
                atom_type = int(parts[1])
                x, y, z = map(float, parts[2:5])
                atoms.append((atom_type, x, y, z))
    
    # 计算盒子向量和参数
    a = np.array([xhi - xlo, 0, 0])
    b = np.array([xy, yhi - ylo, 0])
    c = np.array([xz, yz, zhi - zlo])
    
    # 计算盒子长度
    a_len = np.linalg.norm(a)
    b_len = np.linalg.norm(b)
    c_len = np.linalg.norm(c)
    
    # 计算盒子角度
    alpha = np.degrees(np.arccos(np.dot(b, c) / (b_len * c_len)))
    beta = np.degrees(np.arccos(np.dot(a, c) / (a_len * c_len)))
    gamma = np.degrees(np.arccos(np.dot(a, b) / (a_len * b_len)))
    
    # 写入xyz文件
    with open(os.path.join(work_dir, output_file), 'w') as f:
        f.write(f"{natoms}\n")
        f.write(f"# CELL(abcABC): {a_len:.6f} {b_len:.6f} {c_len:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f} cell{{angstrom}} positions{{angstrom}}\n")
        
        for atom in atoms:
            atom_type, x, y, z = atom
            element = type_map.get(atom_type, 'X')  # 使用X表示未知元素
            f.write(f"{element} {x:.8f} {y:.8f} {z:.8f}\n")


def convert_velocity_ms_to_angps(velocity_ms):
    """
    将速度从m/s转换为Å/ps
    1 m/s = 1e10 Å/s = 1e-2 Å/ps (因为1 ps = 1e-12 s)
    """
    return float(velocity_ms) * 1e-2

def extract_cell_info(comment_line):
    """
    从注释行提取CELL信息
    
    注释行格式示例:
    # CELL(abcABC):   22.72307    22.78546    21.31459    90.00000    90.00000   120.03664  Step:           0
    
    返回: (a, b, c, alpha, beta, gamma)
    """
    # 使用正则表达式提取数字
    numbers = re.findall(r"[-+]?\d*\.\d+|\d+", comment_line)
    
    if len(numbers) >= 6:
        try:
            a = float(numbers[0])
            b = float(numbers[1])
            c = float(numbers[2])
            alpha = float(numbers[3])
            beta = float(numbers[4])
            gamma = float(numbers[5])
            return (a, b, c, alpha, beta, gamma)
        except ValueError:
            pass
    
    # 如果无法提取，返回默认值（需要根据您的系统调整）
    print(f"警告: 无法从注释行提取CELL信息，使用默认值: {comment_line}")
    return (20.0, 20.0, 20.0, 90.0, 90.0, 90.0)


def convert_cell_to_lammps(a, b, c, alpha, beta, gamma):
    """
    将晶胞参数转换为LAMMPS的倾斜盒子表示
    
    返回: 三个列表，每个列表包含 [lo, hi, tilt]
    """
    # 将角度转换为弧度
    alpha_rad = alpha * np.pi / 180.0
    beta_rad = beta * np.pi / 180.0
    gamma_rad = gamma * np.pi / 180.0
    
    # 计算倾斜因子
    xy = b * np.cos(gamma_rad)
    xz = c * np.cos(beta_rad)
    yz = (b * c * np.cos(alpha_rad) - xy * xz) / np.sqrt(b*b - xy*xy)
    
    # 计算盒子边界
    # 对于斜盒子，需要确保所有原子都在盒子内
    # 计算最小和最大可能的坐标值
    x_min = min(0.0, xy, xz, xy + xz)
    x_max = max(a, a + xy, a + xz, a + xy + xz)
    
    y_min = min(0.0, yz)
    y_max = max(np.sqrt(b*b - xy*xy), np.sqrt(b*b - xy*xy) + yz)
    
    z_min = 0.0
    z_max = np.sqrt(c*c - xz*xz - yz*yz)
    
    # 调整边界以确保所有原子都在盒子内
    # 对于x方向，如果xy或xz为负，xlo需要调整
    xlo = x_min
    xhi = x_max
    ylo = y_min
    yhi = y_max
    zlo = z_min
    zhi = z_max
    
    return [
        [xlo, xhi, xy],
        [ylo, yhi, xz],
        [zlo, zhi, yz]
    ]
    
    
def trans_xcvc2dump(work_dir, 
                    cut = 1e10,
                    prefix = None, 
                    output_file='dump.lammps',
                    atom_type_map=None):
    """
    合并xc.xyz和vc.xyz文件，创建包含位置和速度信息的文件
    并将速度单位从m/s转换为Å/ps
    """
    
    xc_file = 'xc.xyz'
    vc_file = 'vc.xyz'
    
    if prefix is not None:
        xc_file = prefix+'.'+xc_file
        vc_file = prefix+'.'+vc_file
        
    # 如果没有提供原子类型映射，使用默认映射
    if atom_type_map is None:
        atom_type_map = {'H': 1, 'O': 2, 'C': 3, 'N': 4}  # 根据您的系统修改    
 
    # 检查文件是否存在
    if not os.path.exists( os.path.join(work_dir,xc_file) ) :
        print(f"错误: 文件 {xc_file} 不存在")
        return False
        
    if not os.path.exists( os.path.join(work_dir,vc_file) ):
        print(f"错误: 文件 {vc_file} 不存在")
        return False
    
    # 打开输入和输出文件
    with open(os.path.join(work_dir,xc_file), 'r') as f_xc, open( os.path.join(work_dir,vc_file), 'r') as f_vc, open( os.path.join(work_dir, output_file) , 'w') as f_out:
        frame_count = 0
        start_time = time.time()
        last_update_time = start_time

        while True:

            if frame_count >= cut:
                break

            # 读取xc文件的原子数
            xc_line1 = f_xc.readline().strip()
            if not xc_line1:  # 文件结束
                break
                
            # 读取vc文件的原子数
            vc_line1 = f_vc.readline().strip()
            if not vc_line1:  # 文件结束
                break
            
            # 检查原子数是否一致
            if xc_line1 != vc_line1:
                print(f"错误: 第 {frame_count+1} 帧原子数不一致: xc={xc_line1}, vc={vc_line1}")
                break
                
            # 读取xc文件的注释行
            xc_line2 = f_xc.readline().strip()
            if not xc_line2:
                break
                
            # 读取vc文件的注释行
            vc_line2 = f_vc.readline().strip()
            if not vc_line2:
                break
            
            # 从注释行提取CELL信息（使用xc文件的注释行）
            cell_info = extract_cell_info(xc_line2)
            
            # 写入dump文件头
            f_out.write("ITEM: TIMESTEP\n")
            f_out.write(f"{frame_count}\n")
            f_out.write("ITEM: NUMBER OF ATOMS\n")
            f_out.write(f"{xc_line1}\n")
            f_out.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            
            # 写入盒子边界信息
            a, b, c, alpha, beta, gamma = cell_info
            box_bounds = convert_cell_to_lammps(a, b, c, alpha, beta, gamma)
            for dim in box_bounds:
                f_out.write(f"{dim[0]:.8e} {dim[1]:.8e} {dim[2]:.8e}\n")
            
            f_out.write("ITEM: ATOMS id type xu yu zu vx vy vz\n")
            
            # 读取并合并原子数据
            num_atoms = int(xc_line1)
            for atom_id in range(1, num_atoms + 1):
                # 读取坐标数据
                xc_atom = f_xc.readline().split()
                if len(xc_atom) < 4:
                    print(f"错误: 第 {frame_count+1} 帧第 {atom_id} 个原子坐标数据不完整")
                    break
                    
                # 读取速度数据
                vc_atom = f_vc.readline().split()
                if len(vc_atom) < 4:
                    print(f"错误: 第 {frame_count+1} 帧第 {atom_id} 个原子速度数据不完整")
                    break
                
                # 检查原子类型是否一致
                if xc_atom[0] != vc_atom[0]:
                    print(f"警告: 第 {frame_count+1} 帧第 {atom_id} 个原子类型不一致")
                
                # 获取原子类型
                element = xc_atom[0]
                atom_type = atom_type_map.get(element, 1)  # 默认为1
                
                # 提取坐标
                x = float(xc_atom[1])
                y = float(xc_atom[2])
                z = float(xc_atom[3])
                
                # 提取并转换速度单位: m/s -> Å/ps
                vx = convert_velocity_ms_to_angps(vc_atom[1])
                vy = convert_velocity_ms_to_angps(vc_atom[2])
                vz = convert_velocity_ms_to_angps(vc_atom[3])
                
                # 写入dump格式的原子数据
                f_out.write(f"{atom_id} {atom_type} {x:.8f} {y:.8f} {z:.8f} {vx:.8f} {vy:.8f} {vz:.8f}\n")
            
            frame_count += 1

            # 每处理一定数量的帧或一定时间后更新进度
            current_time = time.time()
            if current_time - last_update_time > 2.0 or frame_count % 1000 == 0:  
                elapsed = current_time - start_time
                frames_per_second = frame_count / elapsed if elapsed > 0 else 0
                # 使用回车符\r回到行首，实现单行更新
                sys.stdout.write(f"\r已处理 {frame_count} 帧，耗时 {elapsed:.1f} 秒，平均 {frames_per_second:.1f} 帧/秒")
                sys.stdout.flush()
                last_update_time = current_time

        # 处理完成后显示最终结果
        total_time = time.time() - start_time
        print(f"\n转换完成! 共处理 {frame_count} 帧，总耗时 {total_time:.1f} 秒")
    
    print(f"转换完成! 共处理 {frame_count} 帧")
    return True