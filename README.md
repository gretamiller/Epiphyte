#Epiphyte model

The epiphyte model is a water balance model for vascular epiphytes. To capture the dynamics of the water fluxes through epiphytes, the epiphyte model follows an approach similar to the resistor-capacitor model of the soil-plant-atmosphere continuum, without the soil component, and with the addition of schemes for external tank water storage and absorption of atmospheric water vapor through trichomes and velamen in humid conditions. Given inputs of air temperature, solar radiation, and specific humidity, the model estimates transpiration, carbon assimilation, and other hydraulic and photosynthetic variables. The model is currently parameterized to Guzmania monostachia, an epiphytic tank bromeliad that uses both C3 and CAM photosynthesis.


# Photo3

The epiphyte model is built off of Photo3 for plants that root in the soil and uses much of the same code. To download the Photo3 model, go to https://github.com/samhartz/Photo3.
An academic article describing the Photo3 model details is available in Ecological Modelling:
  *Hartzell, S., Bartlett, M.S. and A. Porporato (2018) Unified representation of the C3, C4, and CAM photosynthetic pathways with the Photo3 model. Ecological Modelling, doi: 10.1016/j.ecolmodel.2018.06.012.*

The Photo3 model describes C3, C4, and CAM photosynthesis in a consistent manner using a model which is built on the Farquhar et al. model for carbon assimilation. The model incorporates soil and atmospheric conditions through a representation of the soil-plant-atmosphere continuum. Given soil moisture, air temperature, humidity and solar radiation, the model calculates net carbon assimilation and transpiration, among other variables of interest. Photo3 is currently parameterized for three representative species, one from each photosynthetic type: *Triticum aestivum* (C3), *Sorghum bicolor* (C4), and *Opuntia ficus-indica* (CAM).

# Model Execution

To run the model, simply download the essential files (main.py, defs.py, functions.py, dics.py) to the same folder and run the main file (main.py). The user chooses a plant species, photosynthesis type, initial conditions, duration of the simulation, and a data file containing weather inputs (solar radiation, temperature, and humidity). Results will be exported to the file location chosen in the user interface, and/or may be viewed directly from the command prompt or IDE. The sample_data folder contains sample weather inputs for a location in Barro Colorado Island, Panama, and the sample_output folder contains results from this simulation. 

# Model Structure

The model is structured in an object-oriented fashion, using mixins to combine different photosynthetic, hydraulic, and atmosphere sub-components. The main.py script is the engine which runs the model, creating a Simulation() object which is updated at each timestep according to the model inputs. The defs.py file contains the definitions of the classes and their functions, the dics.py file contains the model global variables, and the functions.py file contains the model global functions. 


# Model requirements
The epiphyte model was developed in Python 3.7 with the following packages: SciPy, NumPy, Pandas, tkinter, SymPy, Matplotlib. We suggest intalling a Python distribution such as [Anaconda][An] to meet these requirements. 

[An]: https://www.continuum.io/downloads

# Instructions for formatting weather data for the model input
Half hourly data for solar radiation, temperature, and specific humidity is needed to run the epiphyte model. Hourly data may be obtained from a variety of sources and then interpolated to the model timestep of 30 minutes using the script cleanData.py, which will produce a new excel file which can be read into the model. In order to use cleanData.py, the data must first be formatted with the following columns: Year, Month, Day, Hour, Minute, GHI, Temperature, Relative Humidity. Data cannot contain gaps larger than 23 hours.

The final weather data file supplied to the model should have the headings: Temperature, Relative Humidity, GHI, and should match the model timestep of 30 minutes (see example files in the sample_data folder).



