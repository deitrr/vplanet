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
import math
from matplotlib import ticker
from scipy.signal import welch
from scipy.optimize import curve_fit
import pdb 
from statsmodels.nonparametric.smoothers_lowess import lowess

d18O_time, d18O, std_error = np.loadtxt('d18O_data.txt',unpack=True)

#Applying linear regression to individual million year to get ride of noise
zeroth_mil =d18O[:800]
first_mil =d18O[800:1250]
second_mil =d18O[1250:1650]
third_mil = d18O[1650:1850]
fourth_milp = d18O[1850:]

smoothed0 = lowess(zeroth_mil, range(len(zeroth_mil)), frac=0.01)
smoothed1 = lowess(first_mil, range(len(first_mil)), frac=0.01)
smoothed2 = lowess(second_mil, range(len(second_mil)), frac=0.01)
smoothed3= lowess(third_mil, range(len(third_mil)), frac=0.01)
smoothed4= lowess(fourth_milp, range(len(fourth_milp)), frac=0.01)

clean_d18O = list(smoothed0[:, 1])+list(smoothed1[:, 1])+list(smoothed2[:, 1])+list(smoothed3[:, 1])+list(smoothed4[:, 1])

dO_dt = np.diff(clean_d18O)/np.diff(d18O_time)
new_time = (-1*d18O_time[1:])*1000

time = -1*(np.loadtxt("time.txt", unpack=True))

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


for file_list in files:
    # print(len(file_list))
    for f in file_list:
        files_org.append(f)

tot_sinsol = []
y_insol = []

for f in files_org:
    insol = np.loadtxt(f, unpack=True)
    N65_insol = insol[-7]
    yearly = np.sum(N65_insol)
    sum_insol = np.sum(N65_insol[153:245])
    tot_sinsol.append(sum_insol)
    y_insol.append(yearly)


tot_sinsol = np.array(tot_sinsol)/90

fig, ax = plt.subplots(5,1,figsize=(15, 10))
ax[0].plot(time, tot_sinsol)
ax[0].set_xlim([-1e6,0])
ax[0].set_title("Integrated summer insolation over 90 days as a function of time")

ax[1].plot(time, tot_sinsol)
ax[1].set_xlim([-2e6,-1e6])

ax[2].plot(time, tot_sinsol)
ax[2].set_xlim([-3e6,-2e6])

ax[3].plot(time, tot_sinsol)
ax[3].set_xlim([-4e6,-3e6])

ax[4].plot(time, tot_sinsol)
ax[4].set_xlim([-5e6,-4e6])
ax[4].set_xlabel("Time (Myr)")

fig.supylabel(f"W/m$^{2}$")
fig.savefig(path / f"90day_sum_insolation.png", dpi=300)


# ###This is to plot the way that Huybers plotted theirs
avg_insol = []
sum_en = []
num_days =[]

for f in files_org:
    insol = np.loadtxt(f, unpack=True)
    insol[insol<275] = np.nan

    cleanedList = [x for x in insol[-7] if str(x) != 'nan']

    num_days.append(len(cleanedList))
    insol2 = np.nanmean(insol[-7])

    avg_insol.append(insol2)
    energy = np.nansum(insol[-7]*86400)/(1e9) ##giga-joules/w^2

    sum_en.append(energy)

    # pdb.set_trace()

fig, ax = plt.subplots(5,1,figsize=(15, 10))

ax[0].plot(time, avg_insol)
ax[0].set_xlim([-1e6,0])

ax[0].set_title("Integrated summer insolation as a function of time")

ax[1].plot(time, avg_insol)
ax[1].set_xlim([-2e6,-1e6])

ax[2].plot(time, avg_insol)
ax[2].set_xlim([-3e6,-2e6])

ax[3].plot(time, avg_insol)
ax[3].set_xlim([-4e6,-3e6])

ax[4].plot(time, avg_insol)
ax[4].set_xlim([-5e6,-4e6])
ax[4].set_xlabel("Time (Myr)")

fig.supylabel(f"W/m$^{2}$")
fig.savefig(path / f"Huybers_copy.png", dpi=300)


fig, ax = plt.subplots(5,1,figsize=(15, 10))

ax[0].plot(time, sum_en)
ax[0].set_xlim([-1e6,0])

ax[0].set_title("Integrated summer insolation as a function of time")

ax[1].plot(time, sum_en)
ax[1].set_xlim([-2e6,-1e6])

ax[2].plot(time, sum_en)
ax[2].set_xlim([-3e6,-2e6])

ax[3].plot(time, sum_en)
ax[3].set_xlim([-4e6,-3e6])

ax[4].plot(time, sum_en)
ax[4].set_xlim([-5e6,-4e6])
ax[4].set_xlabel("Time (Myr)")

fig.supylabel(f"giga-Joules/m$^{2}$")
fig.savefig(path / f"Summer_energy.png", dpi=300)

#making the total insolation with d180_data on top
fig, (ax0, ax2, ax4, ax6, ax8) = plt.subplots(5, figsize = (15,10),sharey=True)  # share the primary y-axis
ax1 = ax0.twinx()
ax3 = ax2.twinx()
ax5 = ax4.twinx()
ax7 = ax6.twinx()
ax9 = ax8.twinx()
color1= 'tab:red'
color2= 'tab:blue'
#first plot
ax0.set_xlabel('time (Myr)')
ax0.set_ylabel(f"giga-Joules/m$^{2}$")
ax0.plot(time, sum_en, color=color1)
ax0.tick_params(axis='y', labelcolor=color1)
ax0.set_xlim([-1e6,-0])

ax1.set_ylabel(f'd$\delta^{18}O$/dt')
ax1.plot(new_time, dO_dt, color=color2)
ax1.tick_params(axis='y', labelcolor=color2)
ax1.set_xlim([-1e6,-0])
ax1.set_ylim([-0.32,0.32])

#second plot
ax2.set_xlabel('time (Myr)')
ax2.set_ylabel(f"giga-Joules/m$^{2}$")
ax2.plot(time, sum_en, color=color1)
ax2.tick_params(axis='y', labelcolor=color1)
ax2.set_xlim([-2e6,-1e6])


ax3.set_ylabel(f'd$\delta^{18}O$/dt')
ax3.plot(new_time, dO_dt, color=color2)
ax3.tick_params(axis='y', labelcolor=color2)
ax3.set_xlim([-2e6,-1e6])
ax3.set_ylim([-0.2,0.2])
#third plot
ax4.set_xlabel('time (Myr)')
ax4.set_ylabel(f"giga-Joules/m$^{2}$")
ax4.plot(time, sum_en, color=color1)
ax4.tick_params(axis='y', labelcolor=color1)
ax4.set_xlim([-3e6,-2e6])

ax5.set_ylabel(f'd$\delta^{18}O$/dt')
ax5.plot(new_time, dO_dt, color=color2)
ax5.tick_params(axis='y', labelcolor=color2)
ax5.set_xlim([-3e6,-2e6])
ax5.set_ylim([-0.2,0.2])
#fourth plot
ax6.set_xlabel('time (Myr)')
ax6.set_ylabel(f"giga-Joules/m$^{2}$")
ax6.plot(time, sum_en, color=color1)
ax6.tick_params(axis='y', labelcolor=color1)
ax6.set_xlim([-4e6,-3e6])

ax7.set_ylabel(f'd$\delta^{18}O$/dt')
ax7.plot(new_time, dO_dt, color=color2)
ax7.tick_params(axis='y', labelcolor=color2)
ax7.set_xlim([-4e6,-3e6])
ax7.set_ylim([-0.1,0.1])
#fifth plot
ax8.set_xlabel('time (Myr)')
ax8.set_ylabel(f"giga-Joules/m$^{2}$")
ax8.plot(time, sum_en, color=color1)
ax8.tick_params(axis='y', labelcolor=color1)
ax8.set_xlim([-5e6,-4e6])

ax9.set_ylabel(f'd$\delta^{18}O$/dt')
ax9.plot(new_time, dO_dt, color=color2)
ax9.tick_params(axis='y', labelcolor=color2)
ax9.set_xlim([-5e6,-4e6])
ax9.set_ylim([-0.06,0.06])

fig.tight_layout()
fig.savefig(path / f"Sum_en_d18O.png", dpi=300)


# #making the total insolation with d180_data on top
fig, (ax0, ax2, ax4, ax6, ax8) = plt.subplots(5, figsize = (15,10),sharey=True)  # share the primary y-axis
ax1 = ax0.twinx()
ax3 = ax2.twinx()
ax5 = ax4.twinx()
ax7 = ax6.twinx()
ax9 = ax8.twinx()
color1= 'tab:red'
color2= 'tab:blue'
#first plot
ax0.set_xlabel('time (Myr)')
ax0.set_ylabel(f"W/m$^{2}$")
ax0.plot(time, tot_sinsol, color=color1)
ax0.tick_params(axis='y', labelcolor=color1)
ax0.set_xlim([-1e6,-0])

ax1.set_ylabel(f'Days')
ax1.plot(time, num_days, color=color2)
ax1.tick_params(axis='y', labelcolor=color2)
ax1.set_xlim([-1e6,-0])
# ax1.set_ylim([-0.32,0.32])

#second plot
ax2.set_xlabel('time (Myr)')
ax2.set_ylabel(f"W/m$^{2}$")
ax2.plot(time, tot_sinsol, color=color1)
ax2.tick_params(axis='y', labelcolor=color1)
ax2.set_xlim([-2e6,-1e6])


ax3.set_ylabel(f'Days')
ax3.plot(time, num_days, color=color2)
ax3.tick_params(axis='y', labelcolor=color2)
ax3.set_xlim([-2e6,-1e6])
# ax3.set_ylim([-0.2,0.2])
#third plot
ax4.set_xlabel('time (Myr)')
ax4.set_ylabel(f"W/m$^{2}$")
ax4.plot(time, tot_sinsol, color=color1)
ax4.tick_params(axis='y', labelcolor=color1)
ax4.set_xlim([-3e6,-2e6])

ax5.set_ylabel(f'Days')
ax5.plot(time, num_days, color=color2)
ax5.tick_params(axis='y', labelcolor=color2)
ax5.set_xlim([-3e6,-2e6])
# ax5.set_ylim([-0.2,0.2])
#fourth plot
ax6.set_xlabel('time (Myr)')
ax6.set_ylabel(f"W/m$^{2}$")
ax6.plot(time, tot_sinsol, color=color1)
ax6.tick_params(axis='y', labelcolor=color1)
ax6.set_xlim([-4e6,-3e6])

ax7.set_ylabel(f'Days')
ax7.plot(time, num_days, color=color2)
ax7.tick_params(axis='y', labelcolor=color2)
ax7.set_xlim([-4e6,-3e6])
# ax7.set_ylim([-0.1,0.1])
#fifth plot
ax8.set_xlabel('time (Myr)')
ax8.set_ylabel(f"W/m$^{2}$")
ax8.plot(time, tot_sinsol, color=color1)
ax8.tick_params(axis='y', labelcolor=color1)
ax8.set_xlim([-5e6,-4e6])

ax9.set_ylabel(f'Days')
ax9.plot(time, num_days, color=color2)
ax9.tick_params(axis='y', labelcolor=color2)
ax9.set_xlim([-5e6,-4e6])
# ax9.set_ylim([-0.06,0.06])

fig.tight_layout()
fig.savefig(path / f"Wm2_days.png", dpi=300)

##SPECTRAL DENSITY PLOT
NFFT = 2**9
window = scipy.signal.windows.hann(NFFT)
fs = 1/100
##Yearly
ft, pxx = scipy.signal.welch(y_insol, fs=fs, window=window, nfft=NFFT, return_onesided=True)
gamma = ft*pxx
fig , ax = plt.subplots(figsize =(10,10))
ax.semilogy(ft*10, pxx)
ax.set_xlim([0.005,0.05])
ax.set_title("Power Spectrum of Yearly Insolation at 65 N ")
ax.set_xlabel("Frequncy (1/Kyr)")
ax.set_ylabel("Power (dB)")
fig.savefig(path /"Spectral_yearly.png", dpi=300)

##90 day summer
ft, pxx = scipy.signal.welch(tot_sinsol, fs=fs, window=window, nfft=NFFT, return_onesided=True)
gamma = ft*pxx
fig , ax = plt.subplots(figsize =(10,10))
ax.semilogy(ft*10, pxx)
ax.set_xlim([0.005,0.05])
ax.set_title("Power Spectrum of Yearly Insolation at 65 N ")
ax.set_xlabel("Frequncy (1/Kyr)")
ax.set_ylabel("Power (dB)")
fig.savefig(path /"Spectral_90day.png", dpi=300)

