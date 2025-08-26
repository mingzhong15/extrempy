def _get_mass_map(type_map):

    mass_map_ref = {
        'H': 1.00794, 'He': 4.0026,

        'Li': 6.941,  'Be': 9.0122,  'B': 10.811,  'C': 12.0107, 
        'N':  14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797, 

        'Na': 22.9897, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 
        'P':  30.9738, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,

        'K': 39.0983, 'Ca': 40.078, 'Sc': 44.9559, 'Ti': 47.867, 
        'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938, 'Fe': 55.845,
        'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 
        'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216, 'Se': 78.96,
        'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 
        'Y': 88.9059, 'Zr': 91.224, 'Nb': 92.9064, 'Mo': 95.94,
        'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 
        'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.710, 
        'Sb': 121.760, 'Te': 127.60, 'I': 126.904, 'Xe': 131.293, 
        
        'Cs': 132.9055, 'Ba': 137.327, 'La': 138.9055, 'Ce': 140.116, 
        'Pr': 140.9077, 'Nd': 144.242, 'Pm': 145.00, 'Sm': 150.36, 
        'Eu': 151.964,  'Gd': 157.25, 'Tb': 158.9253, 'Dy': 162.50, 
        'Ho': 164.9303, 'Er': 167.259, 'Tm': 168.9342, 'Yb': 173.04,
        'Lu': 174.967,  'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 
        'Re': 186.207,  'Os': 190.23, 'Ir': 192.217,  'Pt': 195.084, 
        'Au': 196.9665, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 
        'Bi': 208.9804, 'Po': 209.0,  'At': 210.0, 'Rn': 222.0, 
        'Fr': 223.0,    'Ra': 226.0,  'Ac': 227.0, 'Th': 232.0381, 
        'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 
        'Am': 243.06138, 'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 
        'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0984, 'No': 259.1010, 
        'Lr': 262.1096, 'Rf': 261.1086, 'Db': 262.1138, 'Sg': 266.1218, 
        'Bh': 264.1254, 'Hs': 265.1388, 'Mt': 266.1516, 'Ds': 269.1649, 
        'Rg': 272.174, 'Cn': 277.1926, 'Nh': 278.1884, 'Fl': 281.1972, 
        'Mc': 282.1935, 'Lv': 285.2078, 'Ts': 286.2102, 'Og': 289.2184
    }

    mass_map = []
    for key in type_map:
        mass_map.append(mass_map_ref[key])

    return mass_map

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

        strs =  'variable        a equal 4.75430  \n'
        strs += 'variable        b equal 10.20100   \n'
        strs += 'variable        c equal 5.98190 \n'

        strs += 'lattice               custom        1.0                                &  \n'
        strs += '                      a1            $a         0.0        0.0          &  \n'
        strs += '                      a2            0.0        $b         0.0          &  \n'
        strs += '                      a3            0.0        0.0        $c           &  \n'
        strs += '                      basis         0.00000    0.00000    0.00000      &  \n'
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
        strs += 'create_atoms          3 box  & \n\n'

        strs += '                     basis   1   1   basis   2  1   basis   3  1   basis   4  1   & \n'
        strs += '                     basis   5   1   basis   6  1   basis   7  1   basis   8  1   & \n\n'
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
        strs += 'variable        b equal 8.5960   \n'
        strs += 'variable        c equal 9.0400 \n'

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