from scipy.optimize import fsolve
from sympy import *
import numpy as np
import pandas as pd
from math import exp, pi, sqrt, log
import matplotlib.pyplot as plt
from dics import *
from functions import *
from defs import *

#model time
spinup  = 2
duration = 7+spinup #days
timestepM = 30 # model change in time at each step (min)
timestepD = 30 # timestep of input data 
dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations

#environmental inputs to model
weatherFile = "sample_data/Panama30_constant_humidity.xlsx"
resultsFile = "sample_output/test.csv" # default value for the location where results are saved
df = pd.read_excel(weatherFile)
tempC = df['Temperature']
taInp = tempC + 273. # convert to K
rh = df['Relative Humidity'] # extracts rh column (%)
psat = A_SAT*np.exp((B_SAT*(tempC))/(C_SAT + tempC)) # saturated vapor pressure in Pa
qaInp = 0.622*rh/100.*psat/P_ATM # needs to be in kg/kg
qaInp = list(qaInp.values)
taInp = list(taInp.values)
phiInp = list(df['GHI'].values)  # extracts global solar radiation column from Excel Worksheet in W/m^2

#enter values manually
sinit = 0.5 #initial soil moisture
vwi  = 1.0 # initial water content in plant capacitance, as a fraction of 1
txi = 1.0 # initial water content in tank, as a fraction of 1

atmosphere = Atmosphere(phiInp[0], taInp[0], qaInp[0])
species = Gmono()
photo = CAM(species, atmosphere)
precip = Precip(DrydownPrecip())
hydro = Epiphyte(species, atmosphere, precip, photo, vwi, txi, spinup)
plant = Simulation(species, atmosphere, precip, photo, hydro)


for i in range(int(steps(duration, timestepM))):

	plant.update(dt, phiInp[i], taInp[i], qaInp[i])

results = plant.output()



data = pd.DataFrame.from_dict(results)
# Save data as pandas file
data.to_pickle("sample_output/resultsFile")
#OR
# Save data as csv file
data.to_csv(resultsFile)


#Plot results, with a spinup time before plotting data

startDay = spinup #this starts the graph on the day after the spinup period
endDay = duration
duration = duration-spinup #this fixes the domain of the plots below with respect to the spinup
dispDuration = endDay-startDay
daySteps = int(60/timestepM*24)
timevec = np.linspace(0,duration,duration*daySteps)
timevecHr = np.linspace(0,duration*24,duration*daySteps)


anp = plt.figure()
plt.title("Carbon assimilation")
plt.xlabel("Time (d)")
plt.ylabel("An (umol/m2/s)")
plt.plot(timevec[0:daySteps*dispDuration], results['a'][daySteps*startDay:daySteps*endDay])
for i in range(0, duration):
 	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.axvspan(duration-.25, duration, facecolor = 'k', alpha = 0.2)
plt.xlim(0, duration)

anp = plt.figure()
plt.title("Transpiration")
plt.xlabel("Time (d)")
plt.ylabel("T (mm/d)")
plt.plot(timevec[0:daySteps*dispDuration], [x*3.6*24 for x in results['t'][daySteps*startDay:daySteps*endDay]])
for i in range(0, duration):
 	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.axvspan(duration-.25, duration, facecolor = 'k', alpha = 0.2)
plt.xlim(0, duration)

anp = plt.figure()
plt.title("Relative Tank Storage")
plt.xlabel("Time (d)")
plt.ylabel("(fraction)")
plt.plot(timevec[0:daySteps*dispDuration], results['tank'][daySteps*startDay:daySteps*endDay])
for i in range(0, duration):
 	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.axvspan(duration-.25, duration, facecolor = 'k', alpha = 0.2)
plt.xlim(0, duration)

anp = plt.figure()
plt.title("Free-surface Evaporation from Tank, per unit ground area")
plt.xlabel("Time (d)")
plt.ylabel("E (mm/d)")
plt.plot(timevec[0:daySteps*dispDuration], [x*3.6*24 for x in results['e'][daySteps*startDay:daySteps*endDay]])
for i in range(0, duration):
 	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.axvspan(duration-.25, duration, facecolor = 'k', alpha = 0.2)
plt.xlim(0, duration)

anp = plt.figure()
plt.title("Leaf Water Potential")
plt.xlabel("Time (d)")
plt.ylabel("psi_l (MPa)")
plt.plot(timevec[0:daySteps*dispDuration], results['psi_l'][daySteps*startDay:daySteps*endDay])
for i in range(0, duration):
 	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.axvspan(duration-.25, duration, facecolor = 'k', alpha = 0.2)
plt.xlim(0, duration)
