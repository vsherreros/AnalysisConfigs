"""
import os
import sys

from coffea.util import load
import click

@click.command()
@click.option("-i", "--inputfile", type=str, help="Input file", required=False)


def read_coffea(inputfile):
	if os.path.isfile( inputfile ): accumulator = load(inputfile)
	variables = accumulator['variables'].keys()
"""

import matplotlib.pyplot as plt
import coffea.processor as processor
from coffea.util import load
import awkward as ak
import numpy as np
import yaml

## Load yaml file data

yaml_file_path = "separation_config.yaml"

with open(yaml_file_path, 'r') as stream:
    config = yaml.safe_load(stream)

signal = config["signal"]
bckg = config["background"]


## Load coffea file 

file_path = "output_all.coffea"

output = load(file_path)

variables = output["variables"].keys()
years = output["datasets_metadata"]["by_datataking_period"].keys() #List of years
hist = output["variables"] #["variables"]["TTbar_dR"]["ttZToQQ"]["ttZToQQ_2018"]


def calculate_S2(signal_pdf, bckg_pdf):
	S2 = 0.0
	for s, b in zip(signal_pdf, bckg_pdf):
		if s != 0 and b != 0:  # To avoid dividing by zero
			S2 += (s - b) ** 2 / (s + b)
	S2 *= 0.5
	return S2

s2_values = {}

for var in variables:
	s2_values[var] = {}  #dict to save the S2 values
	for year in years:
		signal_counts = ak.flatten(ak.flatten((hist[var][signal][f"{signal}_{year}"].values())))
		bckg_counts = ak.flatten(ak.flatten((hist[var][bckg][f"{bckg}_{year}"].values())))

		#signal_range = ak.flatten(ak.flatten(hist[var][signal][f"{signal}_{year}"].axes.edges[2]))
		#width = signal_range[1] - signal_range[0]
		
		signal_pdf = signal_counts / np.sum(signal_counts)
		bckg_pdf = bckg_counts / np.sum(bckg_counts)
		
		S2 = calculate_S2(signal_pdf, bckg_pdf)
		
		s2_values[var][year] = S2
		
		print(f"\n{var}:") 
		print("S2:", S2)

yaml_file_path = "s2_values2.yaml"
with open(yaml_file_path, 'w') as yaml_file:
    yaml.dump(s2_values, yaml_file)
    
    
		




"""
counts = ak.flatten(ak.flatten(hist.values()))
x_range = ak.flatten(ak.flatten(hist.axes.edges[2]))

counts = counts.tolist()
x_range = x_range.tolist()

plt.hist(x_range[:-1], bins=x_range, weights=counts, histtype='step', label='Histograma', edgecolor='black', density=True)
plt.xlabel('Valor del eje x')
plt.ylabel('Frecuencia')
plt.title('Histograma')
plt.legend()
plt.show()
"""
