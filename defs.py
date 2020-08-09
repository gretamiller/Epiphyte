from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *
class Simulation(object):
	def __init__(self, species_cls, atm_cls,  precip_cls, photo_cls, hydro_cls):
		self.species = species_cls #plant species
		self.atm = atm_cls #atmospheric inputs
		self.precip = precip_cls #precipitation scenario
		self.photo = photo_cls #photosynthesis type
		self.hydro = hydro_cls #ephiphte hydraulic model
	def update(self, dt, phi, ta, qa):
		self.atm.update(phi, ta, qa) 
		self.photo.update(self.atm, self.hydro.psi_l, self.hydro.tl, dt)
		self.hydro.update(self.species, self.atm, self.precip, self.photo, dt)
		self.precip.update(dt)
	def output(self):
		out = {}
		out.update(self.photo.output())
		out.update(self.hydro.output())
		out.update(self.precip.output())
		return out

class Atmosphere(object):
	ca = 400.
	def __init__(self, phi, ta, qa):
		self.phi = phi
		self.ta = ta
		self.qa = qa
		self.cs = self.ca
	def update(self, phi, ta, qa):
		self.phi = phi
		self.ta = ta
		self.qa = qa

class Photo(object):

	TO = 293.2 # Reference Temperature for photosynthetic parameters (K)
	KAPPA_2 = .3 # Quantum yield of photosynthesis (mol CO2/mol photon)
	GAMMA_1 = .0451 # Parameter for temp dependence of CO2 compensation point (1/K)
	GAMMA_2 = .000347 # Parameter for temp dependence of CO2 compensation point (1/K^2)
	KC0 = 302. # Michaelis constant for C02 at TO (umol/mol)
	KO0 = 256. # Michaelis constant for 02 at TO (mmol/mol)
	OI = .209  # Oxygen Concentration (mol/mol)
	SVC = 649. # Entropy term for carboxylation (J/mol)
	SVQ = 646. # Entropy term for e-transport (J/mol)
	HKC =  59430. # Activation Energy for Kc (J/mol)
	HKO =  36000. # Activation Energy for Ko (J/mol)
	HKR =  53000. # Activation Energy for Rd (J/mol)
	HDJ = 200000. # Deactivation Energy for Jmax (J/mol)
	HAJ = 50000. # Activation Energy for Jmax (J/mol)
	RD0 = .32 # Standard Dark respiration at 25 C (umol/(m^2s))
	HAV =  72000. # Activation Energy for Vc,max (J/mol)
	HDV =  200000. # Deactivation Energy for Vc,max (J/mol)

	def __init__(self, ptype, species):
		self.ptype = ptype
		if hasattr(species, 'RD0'): self.RD0 = species.RD0
		if hasattr(species, 'HAV'): self.HAV = species.HAV
		if hasattr(species, 'HDV'): self.HDV = species.HDV
		self.VCMAX0 = species.VCMAX0
		self.JMAX0 = species.JMAX0
		self.PSILA0 = species.PSILA0
		self.PSILA1 = species.PSILA1
		self.ared = 1.

	def a_c(self, ci, tl, ared):
	    """Rubisco-limited photosynthetic rate (umol/(m^2s^1))"""
	    return self.v_cmax(tl, ared)*(ci - self.gamma(tl))/(ci + self.k_c(tl)*(1. + (self.OI*1000.)/self.k_o(tl)))
	def v_cmax(self, tl, ared):
	    """Maximum carboxylation rate (umol/(m^2s))"""
	    return ared*self.VCMAX0*exp(self.HAV/(R*self.TO)*(1. - self.TO/tl))/(1. + exp((self.SVC*tl - self.HDV)/(R*tl)))
	def k_o(self, tl):
	    """Michaelis-menten coefficient for O2"""
	    return self.KO0*exp(self.HKO/(R*self.TO)*(1. - self.TO/tl))
	def k_c(self, tl):
		"""Michaelis-menten coefficient for CO2"""
		return self.KC0*exp(self.HKC/(R*self.TO)*(1. - self.TO/tl))
	def a_q(self, phi, ci, tl):
	    """Light-limited photosynthetic rate (umol/(m^2s^1))"""
	    return (self.j(phi, tl)*(ci - self.gamma(tl)))/(4.*(ci + 2.*self.gamma(tl)))
	def gamma(self, tl):
	    """CO2 compensation point (umol/mol)"""
	    return self.GAMMA_0*(1. + self.GAMMA_1*(tl - self.TO) + self.GAMMA_2*(tl - self.TO)**2.);
	def jmax(self, tl):
	    """Max. e- transport rate (umol/(m^2s))"""
	    return self.JMAX0*exp(self.HAJ/(R*self.TO)*(1. - self.TO/tl))/(1. + exp((self.SVQ*tl - self.HDJ)/(R*tl))) 
	def j(self, phi, tl):
	    """Electron transport rate (umol/(m^2s))"""
	    return min((phi*10.**6)/(EP*NA)*self.KAPPA_2*.5, self.jmax(tl)) 
	def jpar(self, phi, tl):
	    """Electron transport rate (umol/(m^2s), based off of PAR, not total solar radiatoion)"""
	    return min(phi*self.KAPPA_2, self.jmax(tl)) 
	def a_phiciTl(self, phi, ci, tl, ared):
	    """Net photosynthetic demand for CO2 (umol/(m^2s^1))"""
	    return max(min(self.a_c(ci, tl, ared), self.a_q(phi, ci, tl)),0)
	def a_psilc02(self, psi_l):  
	    """Vulnerability curve for water potential (-)"""
	    if psi_l < self.PSILA0:
	        return 0.
	    elif self.PSILA0 <= psi_l <= self.PSILA1 :
	        return (psi_l - self.PSILA0)/(self.PSILA1  - self.PSILA0)
	    else: 
	        return 1.
	def r_d(self, tl):
	    """Dark respiration flux (umol/(m^2s))"""
	    return self.RD0*exp(self.HKR/(R*self.TO)*(1. - self.TO/tl))
	def csNew(self, an):
		"""CO2 concentration at leaf surface (ppm)"""
		return ca - an/self.GA
	def ciNew(self, cs, ta, qa):
	    """CO2 concentration in mesophyll cytosol (ppm)""" 
	    return cs*(1.-1./(self.A1*self.fD(VPD(ta, qa))))  
	def cmNew(self, cs, ta, qa):
	    """CO2 concentration in mesophyll (ppm)"""
	    return self.ciNew(cs, ta, qa) 
	def fD(self, vpd):
	    """Stomatal response to vapor pressure deficit (-)"""
	    if vpd < 0.1:
	        return 1.
	    else:
	        return 3/13./sqrt(vpd/1000.)
	def gsc(self, phi, ta, psi_l, qa, tl, cx, ared, **kwargs):
	    """Stomatal conductance to CO2, per unit leaf area (mol/m2/s)"""
	    if self.an(phi, psi_l, tl, cx, ared, **kwargs) < 0.:
	        return 0.
	    else:
	        return self.A1*self.an(phi, psi_l, tl, cx, ared, **kwargs)/self.ca*self.fD(VPD(ta, qa))


class C3(Photo):
	A1 = 15.
	GAMMA_0 = 34.6
	RC = 0.7
	GMGSRATIO = 1.65
	def __init__(self, species, atm):
		Photo.__init__(self, "C3", species)
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.cx = self.cm
		self.a_a = []
	def update(self, atm, psi_l, tl, dt):
	 	self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
	 	self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
	 	self.a = self.an(atm.phi, psi_l, tl, self.cm, self.ared)
	 	self.a_a.append(self.a)
	def output(self):
		return {'a': self.a_a}
	def an(self, phi, psi_l, tl, ci, ared): 
		"""Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
		return self.a_psilc02(psi_l)*self.a_phiciTl(phi, ci, tl, ared)

class C4(Photo):
	A1 = 0.5*15.
	GAMMA_0 = 10.
	RC = 0.4
	GMGSRATIO = 2.65
	GBS = .013 # Conductance between bundle sheath and mesophyll (mol m^-2s^-1)
	VPMAX0 = 120. # Maximum PEP carboxylase under well-watered conditions (umol/(m^2s))
	VPR = 80. # PEP regeneration rate (umol/(m^2s))
	KP = 80. # Michaelis-Menten coefficient for C4 species (umol/(mol))
	def __init__(self, species, atm):
		Photo.__init__(self, "C4", species)
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.cbs = self.cm
		self.cx = self.cbs
		self.a_a = []

	def update(self, atm, psi_l, tl, dt):
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		self.cbs = self.cbsNew(self.an(atm.phi, psi_l, tl, self.cbs, self.ared), self.cm, psi_l)
		self.cx = self.cbs
		self.a = self.an(atm.phi, psi_l, tl, self.cbs, self.ared)
		self.a_a.append(self.a)

	def output(self):
		return {'a': self.a_a}

	def an(self, phi, psi_l, tl, ci, ared): 
		"""Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
		return self.a_psilc02(psi_l)*self.a_phiciTl(phi, self.cbs, tl, ared)
	def v_p(self, psi_l, ci):
		"""CO2 concentrating flux (umol/m2/s)"""
		return min((ci*self.VPMAX0)/(ci + self.KP), self.VPR)
	def cbsNew(self, an, cm, psi_l):
		"""CO2 concentration in bundle sheath cell (ppm)"""
		return (self.v_p(psi_l, cm) - an)/self.GBS + cm

class CAM(Photo):
	A1 = 0.6*15.
	GAMMA_0 = 34.6
	RC = 0.5
	GMGSRATIO = 1.
	TR = 90.; # Relaxation time for circadian oscillator (min)
	C0 = 3000. # parameter for decarboxylation of malic acid (umol/mol)
	ALPHA_1 = 1/100.
	ALPHA_2 = 1/7. 
	K = .003 
	TOPT = 288.65 # (K)
	VCM = 0.0027 # Value controlling relative storage of malate (m)
	MU = .5 # Circadian oscillator constant
	BETA = 2.764 # Circadian oscillator constant
	CIRC_1 = .365 # Circadian oscillator constant
	CIRC_2 = .55 # Circadian oscillator constant
	CIRC_3 = 10. # Circadian oscillator constant
	Z0 = .55  # Initial value of z (-)
	M0 = 0. # Initial Malic Acid Carbon Concentration (umol/m^3)
	TH = 302.65 # High temperature for CAM model (K)
	TW = 283.15 # Low temperature for CAM model (K)
	def __init__(self, species, atm):
		Photo.__init__(self, "CAM", species)
		self.MMAX = species.MMAX
		self.AMMAX = species.AMMAX
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.z = self.Z0
		self.m = self.M0
		self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
		self.cx = self.cc
		self.a_a = []
		
	def update(self, atm, psi_l, tl, dt):
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
		self.cx = self.cc
		self.z = self.zNew(atm.phi, self.m, self.z, tl, dt) 
		self.m = self.mNew(atm.phi, psi_l, self.cc, tl, self.z, self.m, dt)
		self.a = self.an(atm.phi, psi_l, tl, self.cc, self.ared)
		self.a_a.append(self.a)

	def output(self):
		return {'a': self.a_a}

	def a_sc(self, phi, psi_l, tl, ci, z, m, ared):
	    """Flux from stomata to Calvin cycle (umol/(m^2s))"""
	    return max(0, self.a_psilc02(psi_l)*(self.a_phiciTl(phi, ci, tl, ared) - self.r_dc(phi, tl))*(1. - self.f_c(z, m)))
	def r_dv(self, phi, tl):
	    """Flux of dark respiration to vacuole (umol/(m^2s))"""
	    return self.r_d(tl)*exp(-phi)
	def r_dc(self, phi, tl):
	    """Flux of dark respiration to calvin cycle (umol/(m^2s))"""
	    return self.r_d(tl)*(1. - exp(-phi))
	def f_o(self, z):
	    """Circadian order function (-)"""
	    return exp(-(z/self.MU)**self.CIRC_3)
	def f_m(self, z, m, tl):
	    """Malic acid storage function"""
	    return self.f_o(z)*(self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m)/(self.ALPHA_2*self.MMAX*((self.TH - tl)/\
	        (self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) + (self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m))
	def f_c(self, z, m):
	    """Carbon circadian control function"""
	    return (1. - self.f_o(z))*m/(self.ALPHA_1*self.MMAX + m)
	def a_sv(self, phi, tl, psi_l, z, m):
		"""Flux from stomata to vacuole (umol/(m^2s))"""
		if self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) > m and (1. - self.K*(tl - self.TOPT)**2.) >0:
			return (self.AMMAX*(1. - self.K*(tl - self.TOPT)**2.) - self.r_dv(phi, tl))*self.f_m(z, m, tl)*self.a_psilc02(psi_l)
		else:
			return 0.
	def a_vc(self, phi, cc, tl, z, m):
	    """Flux from vacuole to calvin cycle (umol/(m^2s))"""
	    return (self.a_phiciTl(phi, cc, tl, 1.) - self.r_dc(phi, tl))*self.f_c(z, m)
	def m_e(self, z, m, tl, phi): 
	    """Malic acid equilibrium value"""
	    if phi>0.:
	        return self.MMAX*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
	            self.CIRC_2*(self.TH - tl)/(self.TH - self.TW) -(1- self.f_o(z))*(1-m/(m+self.ALPHA_1*self.MMAX))) 
	    else:
	        return self.MMAX*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
	            self.CIRC_2*(self.TH - tl)/(self.TH - self.TW)+ (1-self.f_o(z))) 
	def zNew(self, phi, m, z, tl, dt):
	    return max(0, dt*(m - self.m_e(z, m, tl, phi))/(self.MMAX*60.*self.TR) + z)
	def an(self, phi, psi_l, tl, ci, ared): 
	    """Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
	    return self.a_sc(phi, psi_l, tl, ci, self.z, self.m, ared) + self.a_sv(phi, tl, psi_l, self.z, self.m) 
	def ccNew(self, cs, ta, qa, z, m):
		"""CO2 concentration in mesophyll cytosol resulting from malic acid decarboxylation (ppm)"""
		return self.cmNew(cs, ta, qa) + self.f_c(z, m)*self.C0
	def mNew(self, phi, psi_l, cc, tl, z, m, dt): 
		"""Malic acid concentration"""
		return max(((dt/ self.VCM)*(self.a_sv(phi, tl, psi_l, z, m) - self.a_vc(phi, cc, tl, z, m) + self.r_dv(phi, tl))) + m, 0.)

class Precip(object):
	def __init__(self, dynamics):
		self.r = 0.
		self.r_a = []
		self.dynamics = dynamics
	def update(self, dt):
		self.r = self.dynamics.rain(dt)
		self.r_a.append(self.r)
	def output(self):
		return {'r': self.r_a}

class DrydownPrecip(object):
	def __init__(self):
		pass
	def rain(self, dt):
		return 0.

class StochasticPrecip(object):
	"""takes alpha in cm, lda in 1/d"""
	def __init__(self, alpha, lda):
		self.alpha = alpha
		self.lambda_r = lda
	def rain(self, dt):
		if np.random.random() > self.lambda_r*dt/(3600.*24.):
			return 0.
		else:
			return np.random.exponential(self.alpha)

class Epiphyte(object): #for epiphytes
	F_CAP = 0.5 # "f" relative height of water storage xylem connection within plant [-]
	def __init__(self, species, atm, precip, photo, vwi, txi, spinup):
		self.GPMAX = species.GPMAX #max plant conductance [um/MPa s]
		self.GA = species.GA #atmospheric conductance [mm/s]
		self.GCUT = species.GCUT # cuticular conductance, used to find stomatal conductance [mm/s]
		self.GH = species.GH #absorbtion of atmospheric humidity conductance [mm/s]
		self.GWMAX = species.GWMAX #storage conductance [um/MPa s]
		self.ZW = species.VWT #max storage depth [m]
		self.ZT = species.ZT # max tank depth [m]
		self.H = species.H
		self.J = species.J
		self.d1 = species.d1
		self.d2 = species.d2
		self.ns = species.ns
		self.PSI_T = 0. #tank water potential [MPa]
		self.LAI = species.LAI # leaf area index [-]
		self.TAI = species.TAI #tank area index [-]
		self.tx = txi #fraction of water in tank [-]
		self.vw = vwi*self.ZW #initial depth of water in storage [m]
		self.gp = species.GPMAX #first initial gp since we dont know intial psi_l [um/MPa s]
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (photo, atm.phi, atm.ta, atm.qa, photo.cx, self.LAI, self.gp, photo.ared, self.vw)) #[MPa, K]
		self.gp = self.gpf(self.psi_l, self.H, self.J) #plant conductance [um/MPa s]
		self.gt = self.gtf() #tank conductance [um/MPa s] 
		self.t = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.LAI, photo.ared)#T transpiration [um/s]
		self.th = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.LAI, photo.ared)-self.hum(self.LAI, self.tl, self.psi_l, atm.qa)#T-qh [um/s]
		self.e = self.evap(self.LAI, self.tl, atm.ta, atm.qa) #E free surface evaporation from tank [um/s]
		
		self.tx_a = []#fraction of water in tank [-]
		self.vw_a = []#depth of water in storage [m]
		self.gp_a = []#plant conductance [um/MPa s]
		self.gsv_a = [] #stomatal conductance [mol/m2 s]
		self.e_a = []#E free surface evaporation from tank [um/s]
		self.t_a = []#T transpiration from plant [um/s]
		self.qh_a = []#absorption from atmospheric water vapor [um/s]
		self.th_a = []#T-qh [um/s]
		self.qt_a = []#uptake from tank (soil) [um/s]
		self.qw_a = []#qw flux from storage [um/s]
		self.tl_a = [] #leaf temperature[K]
		self.psi_l_a = [] #leaf water potential [MPa]
		self.ql_a = []#leaf specific humidity [kg/kg]
		self.psi_w_a = []#plant storage water potential [MPa]
		self.diff_a = []#difference between atmospheric humidity and leaf humidity
		self.ta_a = []#atmospheric temperature [K]
		self.qa_a = []#atmospheric humidity [kg/kg]
		self.phi_a = []#solar radiation [W/m2]

		self.spinup = spinup*48 #spinup(days) * number of timesteps per day
		self.i = 0 #i is counter for spinup time

	def update(self, species, atm, precip, photo, dt):
		self.i = self.i + 1
		if self.i <= self.spinup:
			self.vw = self.vw #depth of water in storage [m]
			self.tx = self.tx #fraction of water in tank [-]
		else:
			self.vw = self.vwf(self.vw, self.th, self.gp, self.psi_l, self.LAI, dt)#depth of water in storage [m]
			self.tx = self.tnk(dt, precip.r, self.qt, self.tl, atm.ta, atm.qa)#fraction of water in tank [-]
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (photo, atm.phi, atm.ta, atm.qa, photo.cx, self.LAI, self.gp, photo.ared, self.vw))#[MPa, K]
		self.gp = self.gpf(self.psi_l, self.H, self.J)#plant conductance [um/MPa s]
		self.gt = self.gtf()#tank conductance [um/MPa s] 
		self.gs = self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared)#stomatal conductance [mol/m2 s]
		self.psi_w = self.psi_wf(self.vw, self.d1, self.d1, self.ns, self.tl)#plant storage water potential [MPa]
		self.t = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.LAI, photo.ared)#T transpiration [um/s]
		self.qh = self.hum(self.LAI, self.tl, self.psi_l, atm.qa)#absorption from atmospheric water vapor [um/s]
		self.th = self.t-self.qh #T-qh [um/s]
		self.e = self.evap(self.LAI, self.tl, atm.ta, atm.qa)#E free surface evaporation from tank [um/s]
		self.qt = self.qtf(self.vw, self.th, self.gp, self.psi_l, self.LAI, dt)#uptake from tank [um/s]
		self.qw= self.qwf(self.vw, self.th, self.gp, self.psi_l, self.LAI, dt)#qw flux from storage [um/s]
		self.diff = max((atm.qa - self.qi(self.tl, self.psi_l)),-0.0001)#difference between atmospheric humidity and leaf humidity
		self.ta = atm.ta#atmospheric temperature [K]
		self.qa = atm.qa#atmospheric humidity [kg/kg]
		self.phi = atm.phi#solar radiation [W/m2]

		self.tx_a.append(self.tx)#fraction of water in tank [-]
		self.vw_a.append(self.vw)#depth of water in storage [m]
		self.gp_a.append(self.gp)#plant conductance [um/MPa s]
		self.gsv_a.append(self.gs)#stomatal conductance [mol/m2 s]
		self.e_a.append(self.e)#E free surface evaporation from tank [um/s]
		self.t_a.append(self.t)#T transpiration [um/s]
		self.qh_a.append(self.qh)#absorption from atmospheric water vapor [um/s]
		self.th_a.append(self.th)#T-qh [um/s]
		self.qt_a.append(self.qt)#uptake from tank (soil) [um/s]
		self.qw_a.append(self.qw)#qw flux from storage [um/s]
		self.tl_a.append(self.tl)#leaf temp [K] 
		self.psi_l_a.append(self.psi_l)#leaf water potential [MPa]
		self.ql_a.append(self.qi(self.tl, self.psi_l))#leaf specific humidity [kg/kg]
		self.psi_w_a.append(self.psi_w)#plant storage water potential [MPa]
		self.diff_a.append(self.diff)#difference between atmospheric humidity and leaf humidity (for plotting to compare to transporation)
		self.ta_a.append(self.ta)#atmospheric temperature [K]
		self.qa_a.append(atm.qa)#atmospheric humidity [kg/kg]
		self.phi_a.append(self.phi)#solar radiation [W/m2]

	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 't': self.t_a, 'vw': self.vw_a, \
		  'tank':self.tx_a, 'e':self.e_a, 'qw':self.qw_a, 'ql': self.ql_a,'psi_w':self.psi_w_a, 'qt':self.qt_a,\
		  'qa':self.qa_a, 'ta':self.ta_a, 'phi':self.phi_a, 'qh':self.qh_a, 'diff':self.diff_a, 'th':self.th_a}

	def evap(self, lai, tl, ta, qa):
	    """Free-surface evaporation from tank, per unit ground area [um/sec]"""#average ta and tl!!
	    fsev = max(self.TAI*self.GA*1000.*RHO_A/RHO_W*((.622*esat((tl+ta)/2)/P_ATM)-qa), 0.)
		#if the amount of water in tank is less than amount that will evaporate in timestep dt, then what's left will evaporate 
	    if self.tx*self.ZT*self.TAI*10**6 <= 0:
	        return 0.
	    elif self.tx*self.ZT*self.TAI*10**6 <= fsev*dt:
	        return (self.tx*self.ZT*self.TAI*10**6/dt)
	    else:
	        return fsev
	def qtf(self, vw, th, gp, psi_l, lai, dt):
	    """Uptake from tank, per unit ground area [um/s]"""
		#if the amount of water in tank is less than amount that will be absorbed by plant in timestep dt, then what's left will be absorbed 
	    qtt = th - self.qwf(vw, th, gp, psi_l, lai, dt)
	    if self.tx*self.ZT*10**6 <= 0:
	        return 0.
	    elif self.tx*self.ZT*10**6 <= qtt*dt:
	        return (self.tx*self.ZT*10**6/dt)
	    else:
	        return qtt
	def tLoss(self, dt, qt, tl, ta, qa):
	    """Losses from tank (uptake and free surface evap) [um/s]"""
	    return 1./(10.**6)*(qt+(self.evap(self.LAI, tl, ta, qa)))
	def tGain(self, dt, r):
	    """Input to tank (rainfall) [um/s]"""
	    return 1./(10.**6)*r#*precip.r
	def tnk(self, dt, r, qt, tl, ta, qa):
	    """Tank storage fraction[-]"""
	    return min(1., dt/(self.ZT*self.TAI)*(self.tGain(dt, r) - self.tLoss(dt, qt, tl, ta, qa)) + self.tx)
	def evf(self, photo, phi, ta, psi_l, qa, tl, ci, lai, ared, **kwargs):
	    """Transpiration from plant, per unit ground area (um/sec)"""
	    return max(lai*(1./(self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared, **kwargs)*R*ta/P_ATM*1000000.)+1./(self.GA*1000.))**(-1.)\
	    *RHO_A/RHO_W*(self.qi(tl, psi_l)-qa), 0.)
	def hum(self, lai, tl, psi_l, qa): 
	    """Absorption of atmospheric humidity, per unit ground area [um/sec]"""
	    return max((1./((1./(self.GA*1000.))+(1./(self.GH*1000.))))*lai*1000.*RHO_A/RHO_W*(qa-self.qi(tl, psi_l)), 0.)
	def qi(self, tl, psi_l):
	    """Specific humidity internal to leaf (kg/kg)"""
	    try: 
	        ans =  .622*esat(tl)/P_ATM*exp(psi_l*1000000.*VW/R/tl)
	    except OverflowError:
	        ans = 0.
	    return ans
	def gpf(self, psi_l, H, J):
	    """Plant conductance, per unit leaf area (um/(s-MPa))"""
	    if psi_l<-10:
	        return 0.
	    else:
	        return self.GPMAX*exp(-(-psi_l/J)**H)
	def shf(self, tl, ta, lai):
		"""Sensible heat flux (W/m^2), per unit ground area"""
		return CP_A*RHO_A*self.GA*(tl-ta)/1000.*lai
	def gtf(self):
	    """Tank Conductance, per unit ground area (um/(s-MPa))"""
	    #if tank is empty, conductance is 0
	    if self.tx <= 0:
	        return 0.
		#returns 0.5, as a function of TAI
	    else:
	        return 0.5
	def gsw(self, photo, phi, ta, psi_l, qa, tl, ci, ared): 
	    """Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
	    return photo.gsc(phi, ta, psi_l, qa, tl, ci, ared)*1.6 + (self.GCUT*P_ATM/(1000.*R*ta))
	def psi_wf(self, vw, d1, d2, ns, tl): 
	    """Water potential of stored water [MPa]"""
	    osmotic = (R*299./VW)*np.log((((vw/self.ZW)*self.ZW)/(VW))/((((vw/self.ZW)*self.ZW)/(VW))+ns))/10**6 #MPa
	    turgor = ((vw/self.ZW) - d1)**d2#MPa
	    return turgor+osmotic #MPa 
	def vwf(self, vw, ev, gp, psi_l, lai, dt):
	    """Stored water depth (volume, per unit leaf area) [m]"""
	    return min(vw - self.qwf(self.vw, self.th, self.gp, self.psi_l, self.LAI, dt)*dt/10.**6, self.ZW)
	def psi_xf(self, ev, gp, psi_l):
	    """Water potential at connection node x [MPa]"""
	    return ev*(1. - self.F_CAP)/(lai*gp) + psi_l
	def qwf(self, vw, ev, gp, psi_l, lai, dt):
	    """Stored water flux, per unit ground area [um/s]"""
		#if the amount of water in storage is less than amount that will be absorbed by plant in timestep dt, then what's left will be absorbed 
	    qw = (self.gwf(self.psi_wf(self.vw,self.d1, self.d1, self.ns, self.tl), self.H, self.J)*(self.psi_wf(self.vw, self.d1, self.d1, self.ns, self.tl) - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*lai)
	    if self.vw == 0:
	        return 0.
	    elif self.vw*10**6 <= qw*dt:
	        return (self.vw*10**6/dt)
	    else:
	        return qw
	def gwf(self, psi_w, H, J):
	    """Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
	    if self.vw <= 0.:
	        return 0.00001
	    else:
	        return self.GWMAX*exp(-(-psi_w/J)**H)
	def gsrfp(self, gp, lai):
	    """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
	    return (lai*self.gtf()*gp/self.F_CAP)/(self.gtf() +  lai*gp/self.F_CAP)
	def fBal(self, params, photo, phi, ta, qa, c1, lai, gp, ared, vw):
	    psi_l, tl = params
	    psi_w = self.psi_wf(self.vw, self.d1, self.d1, self.ns, 299.)
	    if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops ##not updated for epiphytes yet
	        return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*(self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)+self.evap(lai, tl, qa))/1000000.,  \
	            self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
		           -(self.gsrfp(gp, lai)*(0 - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
		           (1. + (self.gsrfp(gp, lai)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))
	    else:
	        return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*(self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)-self.hum(lai, tl, psi_l, qa))/1000000., #+self.evap(lai, tl, qa)
		        (self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared))-self.hum(lai, tl, psi_l, qa)
		           -(self.gsrfp(gp, lai)*(0 - psi_l) + lai*self.gwf(psi_w, self.H, self.J)*(psi_w - psi_l))/
		           (1. + (self.gsrfp(gp, lai)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w, self.H, self.J)*(1. - self.F_CAP))/gp))



class Gmono(object):
	NAME = 'G. mono'
	TAI = 0.17 #tank area index [-] 
	ZT = 0.015 #max tank depth [m] this value is normalized to the tank
	LAI = 3. #leaf area index [-] 
	GCUT = 0.01 #cuticular conductance, to find stomatal conductance [mm/s]
	GH = 0.01 #absorption of atmospheric humidity conductance [mm/s]
	GA = 13 #atmospheric conductance [mm/s]
	GPMAX = .076 # max plant conductance [um/MPa s] 
	GWMAX = 0.0045#storage conductance [um/MPa s]
	VWT = .0027#max plant storage depth [m]
	H = 2.0 #shape parameter for calcuating gp
	J = 2.0 #shape parameter for calcuating gp
	d1 = 0.0; #unitless constant for calculating psi w
	d2 = 2.0; # unitless constant for calculating psi w
	ns = 1.1 #mol of solute for calculating psi w
	VCMAX0 = 13./1.25 #Maximum carboxylation rate (umol/(m^2 s))
	JMAX0 = 26./1.25 #electron transprt rate (umol/(m^2 s))
	PSILA0 = -1.0 # max plant water stress [MPa]
	PSILA1 = -0.3 # onset of plant water stress [MPa]
	MMAX = 190000000./1.25 # max concentration of malic acid (umol/m^3)
	AMMAX = 13.5  # rate of malic acid storage flux (umol/(m^2 s)



class FacCAM(Photo):
	# A1 = 0.6*15.
	GAMMA_0 = 34.6
	RC = 0.5
	GMGSRATIO = 1.
	TR = 90.; # Relaxation time for circadian oscillator (min)
	C0 = 3000. # parameter for decarboxylation of malic acid (umol/mol)
	ALPHA_1 = 1/100.
	ALPHA_2 = 1/7. 
	K = .003 
	TOPT = 288.65 # (K)
	VCM = 0.0027 # Value controlling relative storage of malate (m)
	MU = .5 # Circadian oscillator constant
	BETA = 2.764 # Circadian oscillator constant
	CIRC_1 = .365 # Circadian oscillator constant
	CIRC_2 = .55 # Circadian oscillator constant
	CIRC_3 = 10. # Circadian oscillator constant
	Z0 = .55  # Initial value of z (-)
	M0 = 0. # Initial Malic Acid Carbon Concentration (umol/m^3)
	TH = 302.65 # High temperature for CAM model (K)
	TW = 283.15 # Low temperature for CAM model (K)
	def __init__(self, species, atm):
		Photo.__init__(self, "FacCAM", species)
		self.pmode = "C3"
		self.A1 = 15.
		self.MUPPER = species.MMAX
		self.AMMAX = species.AMMAX
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.z = self.Z0
		self.m = self.M0
		# self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
		self.cx = self.cm
		self.a_a = []
		
	def update(self, atm, psi_l, tl, dt):
		if psi_l < -1.:
			self.pmode = "CAM"
			self.A1 = 0.6*15.
			self.mmax = self.MUPPER
			if psi_l < -2.5:
				self.mmax = MUPPER
		else:
			pass
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		if self.pmode == "CAM":
			self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
			self.z = self.zNew(atm.phi, self.m, self.z, tl, dt) 
			self.m = self.mNew(atm.phi, psi_l, self.cc, tl, self.z, self.m, dt)
		else:
			pass
		self.cx = self.cxNew()
		self.a = self.an(atm.phi, psi_l, tl, self.cx, self.ared)
		self.a_a.append(self.a)

	def output(self):
		return {'a': self.a_a}

	def cxNew(self):
		if self.pmode == "CAM":
			return self.cc
		else:
			return self.cm
	def an(self, phi, psi_l, tl, ci, ared): 
	    """Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
	    if self.pmode == "CAM":
	    	return self.a_sc(phi, psi_l, tl, ci, self.z, self.m, ared) + self.a_sv(phi, tl, psi_l, self.z, self.m) 
	    else:
	    	return self.a_psilc02(psi_l)*self.a_phiciTl(phi, ci, tl, ared)
		
	def a_sc(self, phi, psi_l, tl, ci, z, m, ared):
	    """Flux from stomata to Calvin cycle (umol/(m^2s))"""
	    return max(0, self.a_psilc02(psi_l)*(self.a_phiciTl(phi, ci, tl, ared) - self.r_dc(phi, tl))*(1. - self.f_c(z, m)))
	def r_dv(self, phi, tl):
	    """Flux of dark respiration to vacuole (umol/(m^2s))"""
	    # return self.r_d(tl)*exp(-phi)
	    return 0.
	def r_dc(self, phi, tl):
	    """Flux of dark respiration to calvin cycle (umol/(m^2s))"""
	    # return self.r_d(tl)*(1. - exp(-phi))
	    return 0.
	def f_o(self, z):
	    """Circadian order function (-)"""
	    return exp(-(z/self.MU)**self.CIRC_3)
	def f_m(self, z, m, tl):
	    """Malic acid storage function"""
	    return self.f_o(z)*(self.mmax*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m)/(self.ALPHA_2*self.mmax*((self.TH - tl)/\
	        (self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) + (self.mmax*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m))
	def f_c(self, z, m):
	    """Carbon circadian control function"""
	    return (1. - self.f_o(z))*m/(self.ALPHA_1*self.mmax + m)
	def a_sv(self, phi, tl, psi_l, z, m):
		"""Flux from stomata to vacuole (umol/(m^2s))"""
		if self.mmax*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) > m and (1. - self.K*(tl - self.TOPT)**2.) >0:
			return (self.AMMAX*(1. - self.K*(tl - self.TOPT)**2.) - self.r_dv(phi, tl))*self.f_m(z, m, tl)*self.a_psilc02(psi_l)
		else:
			return 0.
	def a_vc(self, phi, cc, tl, z, m):
	    """Flux from vacuole to calvin cycle (umol/(m^2s))"""
	    return (self.a_phiciTl(phi, cc, tl, 1.) - self.r_dc(phi, tl))*self.f_c(z, m)
	def m_e(self, z, m, tl, phi): 
	    """Malic acid equilibrium value"""
	    if phi>0.:
	        return self.mmax*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
	            self.CIRC_2*(self.TH - tl)/(self.TH - self.TW) -(1- self.f_o(z))*(1-m/(m+self.ALPHA_1*self.mmax))) 
	    else:
	        return self.mmax*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
	            self.CIRC_2*(self.TH - tl)/(self.TH - self.TW)+ (1-self.f_o(z))) 
	def zNew(self, phi, m, z, tl, dt):
	    return max(0, dt*(m - self.m_e(z, m, tl, phi))/(self.mmax*60.*self.TR) + z)
	def ccNew(self, cs, ta, qa, z, m):
		"""CO2 concentration in mesophyll cytosol resulting from malic acid decarboxylation (ppm)"""
		return self.cmNew(cs, ta, qa) + self.f_c(z, m)*self.C0
	def mNew(self, phi, psi_l, cc, tl, z, m, dt): 
		"""Malic acid concentration"""
		return max(((dt/ self.VCM)*(self.a_sv(phi, tl, psi_l, z, m) - self.a_vc(phi, cc, tl, z, m) + self.r_dv(phi, tl))) + m, 0.)