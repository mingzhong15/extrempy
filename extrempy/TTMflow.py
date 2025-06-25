from extrempy.constant import *

std2FWHM = 2 * np.sqrt(2*np.log(2))

def func_heat_capacity(x, a0, a1, a2, T_damp):

    return a0 * x + (a1*x**2+ a2*x**3) * np.exp(-T_damp*x) 

def func_el_ph_coupling(x, gamma0, delta_gamma, T0, a, b):
   
    return gamma0 + delta_gamma * pow(x/T0, a) / ( 1+ pow(x/T0, b) )
    
class TTMMDSys():

    def __init__(self, mass_mole, a, latt_type):

        self.mass = mass_mole / kg2gmol
        self.a = a

        if latt_type == 'bcc':
            self.cell_atom = 2
        elif latt_type == 'fcc':
            self.cell_atom = 4

        self.num_rho = self.cell_atom /(self.a/m2A)**3   # [1/m^3]
        self.rho = self.num_rho * self.mass      # [kg/m^3]  
        print('density %.2f g/cm^3'%(self.rho * kg2g/(m2cm**3)))


        # initilization of parameters
        self.heat_capacity_param = np.array([1,0,0,1])
        self.gei_param = np.array([1,0,1,1,1])

        self.v_fermi = 0.0
        self.A_temp =  1.0
        self.B_temp =  1.0

        self.gamma_s = 0.0
        self.v0 = 0.0

        self.surface_movement = 1
        self.T_min = 300

    def set_system_size(self, L, mode='free', delta =3.55, min_act = 80):
        
        self.L = L
        #r_verlet = 6.10 + 1.0     # unit: angstrom
        #delta = r_verlet / 2
        
        self.Nx = int(self.L * 10 /self.a ) + 1
        self.delta = delta
        self.N_grid = int(L * 10 / self.delta) +1
        
        self.boundary_mode = mode
        if self.boundary_mode == 'free':
            self.N_expand = 3
        elif self.boundary_mode == 'substrate':
            self.N_expand = 2

        print('=====================================')
        print('system size setting starts')
        print('to generate %.2f nm nanofilm'%(self.L), 'with lattice constant %.2f A'%(self.a))
        
        print('lattice grid for supercell (MD): ',self.Nx, 
              '\n tempreature grid for TTM (FD) :', 
              self.N_grid, self.N_grid * self.N_expand, 
              ' (delta = %.2f ang)'% self.delta)
        
        self.l_surface = int(self.N_grid)
        self.r_surface = int(self.N_grid*2-1) 

        self.coords = np.arange(self.N_grid*self.N_expand) * self.delta /10 - L

        print('ngrid left Surface: %.d \n ngrid right Surface: %.d \n'%(self.l_surface, self.r_surface))        
        
        self.min_act = min_act
        
    def _generate_init_TE_grid(self, savedir, T_0):

        grid_x = np.arange(self.N_grid * self.N_expand).astype(int)

        T_init = np.vstack( (grid_x, np.zeros(grid_x.shape), np.zeros(grid_x.shape), np.zeros(grid_x.shape) ))

        thick = self.N_grid

        TE_field = T_init.T.astype(int)

        TE_field[self.N_grid:self.N_grid*2,-1] += T_0

        print('\n from grid %.d to %.d, electronic temperature is initialized at %.1f kelvin \n'%(self.N_grid,self.N_grid*2-1,T_0))
        
        np.savetxt(savedir+'TTM_init.dat',TE_field, delimiter='   ', fmt='%.d')

    def _obtain_heat_capacity(self, data):
    
        popt, pcov = curve_fit(func_heat_capacity, data[:,0], data[:,1])
        
        fig, ax = plt.subplots(figsize=(2,2),dpi=200)
        ax.plot(data[:,0],data[:,1],'o', ms=2, color='k',mfc='none',
                mew=0.5, linewidth=1.0, label='Raw Data')

        X = np.linspace(data[0,0], data[-1,0],100)
        Ce_fit = func_heat_capacity(X, *popt)

        ax.plot(X, Ce_fit,'-',color='blue',mfc='none',mew=0.5, 
                zorder=0, alpha=0.6, linewidth=2.5, label='Fit')
        
        ax.legend(fontsize=6)
        ax.set_xlim(0, data[-1,0])
        ax.set_ylim(0,np.max(data[:,1]))
        ax.set_ylabel('$C_e$ ($J/(m^3\cdot K)$)',fontsize=8)
        ax.set_xlabel('Temperature ($10^3K$)',fontsize=8)

        self.heat_capacity_param = popt
        self.heat_capacity_param[:-1] *= J2eV /(m2A**3)

        print('\n # eletronic heat capacity')
        print('el_heat_capacity.a0 = %.32f'%self.heat_capacity_param[0])
        print('el_heat_capacity.a1 = %.32f'%self.heat_capacity_param[1])
        print('el_heat_capacity.a2 = %.32f'%self.heat_capacity_param[2])
        print('el_heat_capacity.A  = %.32f\n'%self.heat_capacity_param[3])

        return ax

    def _obtain_ele_diff(self, v_fermi, A_temp, B_temp):
        
        self.v_fermi = v_fermi
        self.A_temp = A_temp
        self.B_temp = B_temp
        
        print('\n # eletronic heat conductivity')
        print('el_heat_cond.v_fermi = %.32f'%self.v_fermi)
        print('el_heat_cond.A_temp = %.32f'%self.A_temp)
        print('el_heat_cond.B_temp = %.32f\n'%self.B_temp)

    def _obtain_gamma(self, data, G0, is_const=False):
        
        print('=====================================')
        print('e-ph coupling setting starts')

        fig, ax = plt.subplots(figsize=(2, 2),dpi=200)

        ax.plot(data[:,0],data[:,1], 'o', ms=2, color='k',
                mfc='none',mew=0.5, linewidth=1.0, label='Raw Data')
        
        #unit : kg/s
        gamma_0_SI = self.mass  /(3 * self.num_rho * kb) * G0      
        self.gamma_0_metal = gamma_0_SI * kg2gmol / s2ps

        print('e-p coupling (G0) is transformed into friction parameter: ',self.gamma_0_metal, 'g/mol / ps')

        if is_const:

            self.gamma_param = [self.gamma_0_metal, 0, 1, 1, 1]

            G_0 = self.gamma_param[0] / (kg2gmol / s2ps) * (3 * self.num_rho * kb)  / self.mass  
            ax.axhline(y=G_0, color='red', linestyle='--', 
                       linewidth=1.5, label='Constant G0')

        else:
            # fitting the gamma(Te)

            gamma_temp_metal = self.mass  /(3 * self.num_rho * kb) * data[:,1] * kg2gmol / s2ps    

            popt, pcov = curve_fit(func_el_ph_coupling, data[:,0], gamma_temp_metal, maxfev=100000)
            self.gamma_param = popt

            X = np.linspace(data[0,0], data[-1,0],100)

            gamma_fit = func_el_ph_coupling(X, *popt)
            
            G_temp_fit = gamma_fit / (kg2gmol / s2ps) * (3 * self.num_rho * kb)  / self.mass  

            ax.plot(X, G_temp_fit,'-',color='blue',mfc='none',mew=0.5, zorder=0, alpha=0.6, linewidth=2.5, label='Fit')

        ax.legend(fontsize=6)
        ax.set_xlim(0, data[-1,0])
        ax.set_ylim(0, np.max(data[:,1]))
        ax.set_ylabel('$G (T_e)$ ($W/(m^3\cdot K)$)',fontsize=8)
        ax.set_xlabel('Temperature ($10^3K$)',fontsize=8)
  
        print('e-ph coupling setting ends')
        print('=====================================')
        print('\n # eletronic-phonon coupling')
        print('el_ph_coupling.gamma_0 = %.22f'%self.gamma_param[0])
        print('el_ph_coupling.delta_gamma = %.22f'%self.gamma_param[1])
        print('el_ph_coupling.T0 = %.22f'%self.gamma_param[2])
        print('el_ph_coupling.a = %.22f'%self.gamma_param[3])
        print('el_ph_coupling.b = %.22f\n'%self.gamma_param[4])

        return ax



    def _obtain_laser(self, tau_L, skin_layer=10, is_decay=False, epsilon=0, F_abs=0,  is_print=True):
    
        # ======================
        # tau_L [ps], skin_layer [nm], epsilon [J/kg], F [J/m^2]
        # L     [nm], rho     [kg/m^3],
        # ======================

        self.sigma = tau_L / std2FWHM                   # [ps]
        self.ngrid_skin_layer = int(skin_layer*nm2A/self.delta)     # [int]
        
        self.is_decay = is_decay
        decay_factor = 1 / (1 - np.exp( - self.L / skin_layer))
        
        fluence2intensity = 1/np.sqrt(2*np.pi)/self.sigma     # [ps^-1]
        epsilon2fluence   = self.rho * (self.L* nm2A/m2A)           # [kg/m^2]

        if epsilon != 0 :
            F_abs = epsilon * epsilon2fluence                # [J/m^2] <-- [J/kg] * [kg/m^2]
     
        if self.is_decay:
            self.I_abs = F_abs * fluence2intensity * J2eV/m2A**2 * decay_factor
        else:
            self.I_abs = F_abs * fluence2intensity * J2eV/m2A**2 # [eV/A^2/ps] 
        
        if is_print:
            if self.is_decay:
                print('exponential decay with optical peneration depth : %.d grids'%(self.ngrid_skin_layer))
            else:
                print('uniform heating condition')
                
            print('Intensity (eV/A^2 ps): ', self.I_abs)
            print('laser pulse (ps): ', self.sigma)

        print('\n # laser')
        print('laser.intensity = %.32f'%self.I_abs)
        print('laser.pulse_width = %.32f'%self.sigma)

        print('laser.ngrid_lsurface = %.d'%self.l_surface)
        print('laser.ngrid_rsurface = %.d'%self.r_surface)

        print('laser.ngrid_skin_layer = %.d'%self.ngrid_skin_layer)
        if self.is_decay:
            temp = 1
        else:
            temp = 0
        print('laser.is_decay = %.d \n'%temp)       
       
    
    def _generate_param(self, savedir, savefile):

        file = open( os.path.join(savedir,savefile),'w')

        file.writelines('# eletronic heat capacity \n')
        file.writelines('el_heat_capacity.a0 = %.32f \n'%self.heat_capacity_param[0])
        file.writelines('el_heat_capacity.a1 = %.32f \n'%self.heat_capacity_param[1])
        file.writelines('el_heat_capacity.a2 = %.32f \n'%self.heat_capacity_param[2])
        file.writelines('el_heat_capacity.A  = %.32f \n'%self.heat_capacity_param[3])

        file.writelines('\n # eletronic heat conductivity \n')
        file.writelines('el_heat_cond.v_fermi = %.32f \n'%self.v_fermi)
        file.writelines('el_heat_cond.A_temp = %.32f \n'%self.A_temp)
        file.writelines('el_heat_cond.B_temp = %.32f \n'%self.B_temp)

        file.writelines('\n # eletronic-phonon coupling \n')
        file.writelines('el_ph_coupling.gamma_0 = %.32f \n'%self.gamma_param[0])
        file.writelines('el_ph_coupling.delta_gamma = %.32f \n'%self.gamma_param[1])
        file.writelines('el_ph_coupling.T0 = %.32f \n'%self.gamma_param[2])
        file.writelines('el_ph_coupling.a = %.32f \n'%self.gamma_param[3])
        file.writelines('el_ph_coupling.b = %.22f \n'%self.gamma_param[4])


        file.writelines('\n # eletronic stopping power \n')
        file.writelines('el_stopping.gamma_s = %.32f \n'%self.gamma_s)
        file.writelines('el_stopping.v0 = %.32f \n'%self.v0)

        file.writelines('\n # laser \n')
        file.writelines('laser.intensity = %.32f \n'%self.I_abs)
        file.writelines('laser.pulse_width = %.32f \n'%self.sigma)

        file.writelines('laser.ngrid_lsurface = %.d \n'%self.l_surface)
        file.writelines('laser.ngrid_rsurface = %.d \n'%self.r_surface)

        file.writelines('laser.ngrid_skin_layer = %.d \n'%self.ngrid_skin_layer)
        if self.is_decay:
            temp = 1
        else:
            temp = 0
        file.writelines('laser.is_decay = %.d \n'%temp)
        
        file.writelines('\n # geometry \n')
        file.writelines('geometry.surface_movement = %.d \n'%self.surface_movement)
        file.writelines('geometry.T_min = %.3f \n'%self.T_min)
        file.writelines('geometry.number_particle_min = %.d \n'%self.min_act)

        file.close()
