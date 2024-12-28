import os  
import numpy as np
import matplotlib.pyplot as plt
import time
import glob
import dpdata

from monty.serialization import loadfn,dumpfn
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'

# ==================================== #
#  Constants
# ==================================== #

e = 1.602e-19                           #charge of a electron, C / electron volt, J
epsilon0 = 8.854e-12                    #vaccum dielectric constant, F/m
me = 9.109e-31                          #mass of a electron , kg
h = 6.626e-34                           #Planck constant, J*s 
hbar = h/(2*np.pi)                      #reduced Planck constant, J*s
NA = 6.022e23                           #Avogadro constant, 1/mol
kb = 1.38e-23                            # Boltzmann constant, J/K
aB = 4*np.pi*epsilon0*hbar**2 / (me*e**2)# Bohr radius, m
Ry = 13.6                                # Ryberg, eV

# ==================================== #
#  Real Units <---> S.I. Units
# ==================================== #
cal2J = 4.184 
J2cal = 1/cal2J

kg2g = 1e3
g2kg = 1/kg2g

m2cm = 1e2
cm2m = 1/m2cm

mJ2J = 1e-3        
J2mJ = 1/mJ2J

nm2A = 10
A2nm = 1/nm2A
# ==================================== #
#  Metal Units <---> S.I. Units
# ==================================== #

J2eV = 6.2442e18

m2A = 1e10

gmol2g = 1/NA

kg2gmol = 1e3/gmol2g

s2ps = 1e12
s2fs = 1e15

gcc2kgm3 = 1e-3/(cm2m**3)

bar2Pa = 1e5
kbar2GPa = 1e3*bar2Pa/1e9
Pa2GPa = 1e-9

# ==================================== #
#  Atomic Units <---> Metal/S.I. Units
# ==================================== #

m2bohr = 1/aB
A2bohr = m2bohr/m2A
eV2Ry = 1/Ry

# ==================================== #
#  PHONON
# ==================================== #

c = 3*10**8                               # speed of light, m/s
omega2k = s2ps*1/(c) *cm2m                # 1/ps -> 1/cm
meV2THz = 1e-3 / hbar /J2eV / 2/np.pi / s2ps

# ==================================== #
#  Basic Function 
# ==================================== #

def grep_file(file, keyword):

    f = open(file, 'r')

    lines = f.readlines()

    for ll in lines:

        if keyword in ll:

            return ll
        
    return None

def _estimate_sel(rcut, num_rho):
    return 4/3 * np.pi * rcut**3 * num_rho

def radial_average_2d(x,y):
    
    x_uniq = np.unique(x)
    y_rave = np.zeros(x_uniq.shape)
    
    for i in range(x_uniq.shape[0]):
        xi = x_uniq[i]
        y_rave[i] = np.sum(y[x == xi])
    
    return x_uniq, y_rave

def chunk_average(x, y, N):
    
    x_chunk = np.linspace(x.min(), x.max(), N)
    y_chunk = np.zeros(x_chunk.shape)
    for i in range(N-1):
        ll = x_chunk[i]
        rr = x_chunk[i+1]
        
        cri1 = x > ll
        cri2 = x < rr
        y_chunk[i] = np.average(y[cri1 * cri2])
    
    return x_chunk, y_chunk

def chunk_sum(x, y, N):
    
    x_chunk = np.linspace(x.min(), x.max(), N)
    y_chunk = np.zeros(x_chunk.shape)
    for i in range(N-1):
        ll = x_chunk[i]
        rr = x_chunk[i+1]
        
        cri1 = x > ll
        cri2 = x < rr
        y_chunk[i] = np.sum(y[cri1 * cri2])
    
    return x_chunk, y_chunk

def trapezoidal(x,y):

    delta = x[1:]-x[:-1]
    n = y.shape[0] - 1
    
    height = (y[1:]+y[:-1])/2
    
    return np.sum(delta*height)
# ==================================== #
#  Basic For Plot
# ==================================== #

sci_color = np.array(['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e'])

cc_matrix = np.array([
    ['#FD6D5A', '#FEB40B', '#6DC354', '#994487', '#518CD8', '#443295'],
    ['#264653', '#2A9D8F', '#E9C46A', '#F4A261', '#E76F51', '#253777'],
    ['#4C87D6', '#F38562', '#F2B825', '#D4C114', '#88B421', '#199FE0'],
    ['#037CD2', '#00AAAA', '#927FD3', '#E54E5D', '#EAA700', '#F57F4B'],
]
)

def generate_colormap(N, colormap=plt.cm.rainbow):

    return [colormap(int(x*colormap.N/N)) for x in range(N)]   

def parity_plot(ax, data_dft, data_dp, INTERVAL, cl, ll=''):

    ax.plot(data_dft[::INTERVAL], data_dp[::INTERVAL],'o',color=cl,
            ms=2, mew=0.5, label=ll)

    vmin = np.min(data_dft)
    vmax = np.max(data_dft)

    x1 = np.linspace(vmin,vmax)
    ax.plot(x1,x1,'--k', lw=1)
    ax.set_xlim(vmin,vmax)
    ax.set_ylim(vmin,vmax)

def phonon_plot(FILE, ax, color, ll, lw, mk, INTERVAL, label):
    
    data = np.loadtxt(FILE)
    nbranch = data.shape[1] - 1
    
    k_max = data[-1,0]

    for i in range(1,nbranch+1):

        inver_cm2THz = 1/omega2k
        inver_cm2meV = 1/omega2k*s2ps*hbar*J2eV*1e3*2*np.pi
        #INTERVAL = 5
        if i == 1:
            ax.plot(data[::INTERVAL,0]/k_max,data[::INTERVAL,i]*inver_cm2THz,ls=ll,marker=mk, 
                    markersize=1.5, mew=0.5, color=color, linewidth=lw,  label=label)

        else:
            ax.plot(data[::INTERVAL,0]/k_max,data[::INTERVAL,i]*inver_cm2THz,ls=ll,marker=mk, 
                    markersize=1.5, mew=0.5, color=color, linewidth=lw, )        
              
    k_path = data[::100,0]/k_max
    
    #ax.set_xlim(data[0,0],data[-1,0])
        
    return k_path 


from scipy import interpolate
# =====================================================================
# return `y` value for `x` array, according to reference data `data_ref`
# =====================================================================
def _obtain_uniform(x, x_ref, y_ref):
 
    f = interpolate.interp1d(y_ref, x_ref,kind='slinear')
    return f(x)


def dump_to_xyz(FILE, OUTFILE):

    file = open(FILE, 'r')
    outfile = open(OUTFILE, 'w')

    line = file.readline()

    while line:

        if 'NUMBER OF ATOMS' in line:
            line = file.readline()
            natoms = int(line)
            outfile.write(line)

        if 'BOX BOUNDS' in line:
             
            for kk in range(3):           
                line = file.readline()

                xl, xh = line.split()

                ret = '%.11f %.11f '%(float(xl),float(xh))

                outfile.write(ret)

            outfile.write('\n')    
            
        if 'ATOMS' in line:

            for i in range(natoms):
                line = file.readline()

                output = line.split()[1:]

                #idx = int(output[0])
                x = float(output[1])
                y = float(output[2])
                z = float(output[3])

                outfile.write('%.11f %.11f %.11f \n'%(x,y,z))

        line = file.readline()
    
    file.close()
    outfile.close()