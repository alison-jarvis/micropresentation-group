import numpy as np
import pandas as pd
import networkx as nx
import os

import os 
import pickle

from jauregui.ex_path import load_spectrum, normalize_spectrum, get_absorption_emission_files

emission = []
abs = []


database_path = './Natural Chlorophylls/'
chlorophyll_spectra = pickle.load(open('./natural_chlorophyll_spectra.pkl', 'rb'))
abs_files, ems_files = get_absorption_emission_files(database_path)

print(ems_files)
for file in ems_files:
    nm, intensity = load_spectrum(file)
    intensity_norm = normalize_spectrum(intensity, nm)
    emission.append(intensity_norm)
print(emission)