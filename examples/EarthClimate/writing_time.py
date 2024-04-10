import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import pathlib
import re
import sys
import scipy
from matplotlib import ticker
from scipy.signal import welch

path = pathlib.Path(__file__).parents[0].absolute()
sys.path.insert(1, str(path.parents[0]))

name_of_files_hund = ['/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.0',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.100',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.200',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.300',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.400',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.500',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.600',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.700',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.800',
                 '/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.900']
name_of_files_thous= ['/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.1900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.2900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.3900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.4900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.5900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.6900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.7900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.8900',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9000',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9100',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9200',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9300',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9400',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9500',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9600',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9700',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9800',
'/home/karengarcia/vplanet_runs/run24/SeasonalClimateFiles/solarsys.Earth.DailyInsol.9900']
name_of_files_e04 = []
name_of_files_e05 = []
name_of_files_e06 = []
files_org =[]
for f in glob.glob(str(path / "SeasonalClimateFiles" / "solarsys.Earth.DailyInsol.*e+04")):
    name_of_files_e04.append(f)
for f in glob.glob(str(path / "SeasonalClimateFiles" / "solarsys.Earth.DailyInsol.*e+05")):
    name_of_files_e05.append(f)
for f in glob.glob(str(path / "SeasonalClimateFiles" / "solarsys.Earth.DailyInsol.*e+06")):
    name_of_files_e06.append(f)

name_of_files_thous= sorted(name_of_files_thous)
name_of_files_e04 = sorted(name_of_files_e04)
name_of_files_e05 = sorted(name_of_files_e05)
name_of_files_e06 = sorted(name_of_files_e06)
files = [name_of_files_hund,name_of_files_thous,name_of_files_e04,name_of_files_e05,name_of_files_e06]
time_text = open("time.txt", "a")

for file_list in files:
    print(len(file_list))
    for f in file_list:
        time_text.write(f + '\n')
        files_org.append(f)
time_text.close()