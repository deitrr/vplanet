import numpy as np
import matplotlib.pyplot  as plt
import seaborn as sns
import matplotlib as mpl

sns.set_style("whitegrid")
plt.close('all')

# Set style for plot #
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.labelsize'] = 13
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 13

cmap=plt.get_cmap('nipy_spectral')
A = 0
# 1: plot individual graphs, 0: plot all graphs together
individual = 0
# 1: plot time axis log scale, 0: plot time axis linear
log_plot = 1
Initial_water = 1
# read data
data = np.loadtxt("1800K_1rad.forward")
data_sch = np.loadtxt("1800K_schaefer.forward")
data_4 = np.loadtxt("1800K_2rad.forward")

time        = data[:,0]  # time (yr)
Tpot        = data[:,1]  # Potential temp magma ocean (K)
Tsurf       = data[:,2]  # Surface temp (K)
r_sol       = data[:,3]  # solidification radius (R_earth)
M_water_mo  = data[:,4] # water mass in magma ocean + atmosphere (TO)
M_water_sol = data[:,5] # water mass in solid mantle (kg)
M_O_mo      = data[:,6] # mass of oxygen in magma ocean + atmosphere (kg)
M_O_sol     = data[:,7] # mass of oxygen in solid mantle (kg)
Press_H2O   = data[:,8] # partial pressure water in atmopshere (bar)
Press_O     = data[:,9] # partial pressure oxygen in atmosphere (bar)
M_H_Space   = data[:,10] # partial pressure oxygen in atmosphere (bar)
M_O_Space   = data[:,11] # partial pressure oxygen in atmosphere (bar)
Frac_Fe2O3  = data[:,12] # partial pressure oxygen in atmosphere (bar)
NetFluxAtmo = data[:,13] # atmospheric net flux (W/m^2)
Frac_H2O    = data[:,14] # Water fraction in magma ocean
RadioHeat   = data[:,15] # Radiogenic Heating (W/kg)

time_sch        = data_sch[:,0]  # time (yr)
Tpot_sch        = data_sch[:,1]  # Potential temp magma ocean (K)
Tsurf_sch       = data_sch[:,2]  # Surface temp (K)
r_sol_sch       = data_sch[:,3]  # solidification radius (R_earth)
M_water_mo_sch  = data_sch[:,4] # water mass in magma ocean + atmosphere (TO)
M_water_sol_sch = data_sch[:,5] # water mass in solid mantle (kg)
M_O_mo_sch      = data_sch[:,6] # mass of oxygen in magma ocean + atmosphere (kg)
M_O_sol_sch     = data_sch[:,7] # mass of oxygen in solid mantle (kg)
Press_H2O_sch   = data_sch[:,8] # partial pressure water in atmopshere (bar)
Press_O_sch     = data_sch[:,9] # partial pressure oxygen in atmosphere (bar)
M_H_Space_sch   = data_sch[:,10] # partial pressure oxygen in atmosphere (bar)
M_O_Space_sch   = data_sch[:,11] # partial pressure oxygen in atmosphere (bar)
Frac_Fe2O3_sch  = data_sch[:,12] # partial pressure oxygen in atmosphere (bar)
NetFluxAtmo_sch = data_sch[:,13] # atmospheric net flux (W/m^2)
Frac_H2O_sch    = data_sch[:,14] # Water fraction in magma ocean
RadioHeat_sch   = data_sch[:,15] # Radiogenic Heating (W/kg)

time_4        = data_4[:,0]  # time (yr)
Tpot_4        = data_4[:,1]  # Potential temp magma ocean (K)
Tsurf_4       = data_4[:,2]  # Surface temp (K)
r_sol_4       = data_4[:,3]  # solidification radius (R_earth)
M_water_mo_4  = data_4[:,4] # water mass in magma ocean + atmosphere (TO)
M_water_sol_4 = data_4[:,5] # water mass in solid mantle (kg)
M_O_mo_4      = data_4[:,6] # mass of oxygen in magma ocean + atmosphere (kg)
M_O_sol_4     = data_4[:,7] # mass of oxygen in solid mantle (kg)
Press_H2O_4   = data_4[:,8] # partial pressure water in atmopshere (bar)
Press_O_4     = data_4[:,9] # partial pressure oxygen in atmosphere (bar)
M_H_Space_4   = data_4[:,10] # partial pressure oxygen in atmosphere (bar)
M_O_Space_4   = data_4[:,11] # partial pressure oxygen in atmosphere (bar)
Frac_Fe2O3_4  = data_4[:,12] # partial pressure oxygen in atmosphere (bar)
NetFluxAtmo_4 = data_4[:,13] # atmospheric net flux (W/m^2)
Frac_H2O_4    = data_4[:,14] # Water fraction in magma ocean
RadioHeat_4   = data_4[:,15] # Radiogenic Heating (W/kg)

n_time = len(time)
n_time_sch = len(time_sch)
n_time_4 = len(time_4)
i_end  = n_time-1

M_water_atm = np.zeros(n_time)
M_O_atm     = np.zeros(n_time)

M_water_atm_sch = np.zeros(n_time_sch)
M_O_atm_sch     = np.zeros(n_time_sch)

M_water_atm_4 = np.zeros(n_time_4)
M_O_atm_4     = np.zeros(n_time_4)

# N_H_sol = np.zeros(n_time) # number of H atoms in solid mantle
# N_H_space = np.zeros(n_time) # number of H atoms in solid mantle
# N_H_mo  = np.zeros(n_time) # number of H atoms in liquid mantle
# N_H_atm = np.zeros(n_time) # number of H atoms in atmosphere
# N_O_sol = np.zeros(n_time) # number of O atoms in solid mantle
# N_O_mo  = np.zeros(n_time) # number of O atoms in liquid mantle
# N_O_atm = np.zeros(n_time) # number of O atoms in atmosphere
# N_O_space = np.zeros(n_time) # number of H atoms in solid mantle

# N_H_tot = np.zeros(n_time) # number of O atoms in atmosphere
# N_O_tot = np.zeros(n_time) # number of O atoms in atmosphere

round = 1e45

TO        = 1.39e21         # mass of 1 Terr. Ocean [kg]
AVOGADROCONST = 6.022e23

REARTH = 6.3781e6        # m
MEARTH = 5.972186e24     # kg
BIGG   = 6.67428e-11     # m**3/kg/s**2
r_p    = 1.15*REARTH
m_p    = 1.62*MEARTH
g      = (BIGG * m_p) / (r_p ** 2)

for i in range(n_time):
    M_water_atm[i] = Press_H2O[i]*1e5 * 4 * np.pi * r_p**2 / g
    M_O_atm[i]     = Press_O[i]*1e5 * 4 * np.pi * r_p**2 / g

for i in range(n_time_sch):
    M_water_atm_sch[i] = Press_H2O_sch[i]*1e5 * 4 * np.pi * r_p**2 / g
    M_O_atm_sch[i]     = Press_O_sch[i]*1e5 * 4 * np.pi * r_p**2 / g

for i in range(n_time_4):
    M_water_atm_4[i] = Press_H2O_4[i]*1e5 * 4 * np.pi * r_p**2 / g
    M_O_atm_4[i]     = Press_O_4[i]*1e5 * 4 * np.pi * r_p**2 / g

    # N_H_space[i] = M_H_Space[i] * AVOGADROCONST / (0.001 * round)
    # N_H_sol[i] = 2 * M_water_sol[i]*TO * AVOGADROCONST / (0.018 * round)
    # N_H_mo[i]  = 2 * (M_water_mo[i]*TO - M_water_atm[i]) * AVOGADROCONST / (0.018 * round)
    # N_H_atm[i] = 2 * M_water_atm[i] * AVOGADROCONST / (0.018 * round)
    # N_H_tot[i] = N_H_sol[i] + N_H_mo[i] + N_H_atm[i] + N_H_space[i]
    #
    # N_O_space[i] = M_O_Space[i] * AVOGADROCONST / (0.016 * round)
    # N_O_sol[i] = M_water_sol[i]*TO * AVOGADROCONST / (0.018 * round) \
    #              + M_O_sol[i] * AVOGADROCONST / (0.016 * round)
    # N_O_mo[i]  = (M_water_mo[i]*TO - M_water_atm[i]) * AVOGADROCONST / (0.018 * round) \
    #              + (M_O_mo[i] - M_O_atm[i]) * AVOGADROCONST / (0.016 * round)
    # N_O_atm[i] = M_water_atm[i] * AVOGADROCONST / (0.018 * round) \
    #              + M_O_atm[i] * AVOGADROCONST / (0.016 * round)
    # N_O_tot[i] = N_O_sol[i] + N_O_mo[i] + N_O_atm[i] + N_O_space[i]

# print('Solidification Time           = ',time[n_time-1]*1e-6,  ' Myr')
# print('Water mass locked in mantle   = ',M_water_sol[n_time-1], ' TO')
# print('Oxygen mass locked in mantle  = ',M_O_sol[n_time-1],     ' kg')
# print('Water pressure in atmosphere  = ',Press_H2O[n_time-1],  ' bar')
# print('Oxygen pressure in atmosphere = ',Press_O[n_time-1],    ' bar')
# print('Fe2O3 mass frac in mantle     = ',Frac_Fe2O3[n_time-1])

# integrate volume of magma ocean over time:
# MagmOc_Vol_Time = 0
# for i in range(1,n_time):
#     MagmOc_Vol_Time = MagmOc_Vol_Time + 4./3 * np.pi * (r_p**3 - r_sol[i]**3) * (time[i]-time[i-1])
#
# MagmOc_Vol_Time_sch = 0
# for i in range(1,n_time_sch):
#     MagmOc_Vol_Time_sch = MagmOc_Vol_Time_sch + 4./3 * np.pi * (r_p**3 - r_sol_sch[i]**3) * (time_sch[i]-time_sch[i-1])
#
# MagmOc_Vol_Time_sch = 0
# for i in range(1,n_time_sch):
#     MagmOc_Vol_Time_sch = MagmOc_Vol_Time_sch + 4./3 * np.pi * (r_p**3 - r_sol_sch[i]**3) * (time_sch[i]-time_sch[i-1])

print('Magma ocean solidification time (0.5 Schaefer) = ' , time[n_time-1]/1e6, ' Myr')
print('Magma ocean solidification time (1.0 Schaefer) = ' , time_sch[n_time_sch-1]/1e6, ' Myr')
print('Magma ocean solidification time (2.0 Schaefer) = ' , time_4[n_time_4-1]/1e6, ' Myr')

# print('Time integrated magma ocean volume (Schaefer) = ', MagmOc_Vol_Time_sch, ' m^3*yr')


### Plot ###

fig = plt.figure()
fig.suptitle('GJ1132b: Initial water content '+str(M_water_mo[0])+' terrestrial oceans', fontsize=16, fontweight='bold')
fig.text(x=0.5, y=0.02, s='Different Radiogenic Heating rates: Solid = 0.8 Schaefer, Dashed = 1 Schaefer, Dotted = 1.6 Schaefer', fontsize=14, ha="center")
# plt.text(x=0.5, y=0.88, s= "My title 2 in different size", fontsize=12, ha="center", transform=fig.transFigure)

ax1 = fig.add_subplot(231)
ax1.plot(time*10**-6, Tpot, label='$T_p$', color=cmap(0))
ax1.plot(time*10**-6, Tsurf, label='$T_{surf}$', color=cmap(220))

ax1.plot(time_sch*10**-6, Tpot_sch, linestyle='--', color=cmap(0))
ax1.plot(time_sch*10**-6, Tsurf_sch, linestyle='--', color=cmap(220))

ax1.plot(time_4*10**-6, Tpot_4, linestyle=':', color=cmap(0))
ax1.plot(time_4*10**-6, Tsurf_4, linestyle=':', color=cmap(220))

ax1.legend(loc='best', frameon=True)
ax1.set_ylabel('Temperature (K)')
# ax1.set_xscale('log')
ax1.set_xlim([0,37])
ax1.set_ylim([500,2000])

ax2 = fig.add_subplot(232, sharex=ax1)
ax2.plot(time*10**-6, r_sol/1.15, color=cmap(0))

ax2.plot(time_sch*10**-6, r_sol_sch/1.15, color=cmap(0), linestyle='--')

ax2.plot(time_4*10**-6, r_sol_4/1.15, color=cmap(0), linestyle=':')

ax2.set_ylim([0.989,1])
# ax2.legend(loc='best', frameon=True)
ax2.set_ylabel('Solidification radius ($r_p$)')

ax3 = fig.add_subplot(233, sharex=ax1)
ax3.plot(time*10**-6, M_water_atm/TO, label='atmosphere', color=cmap(0))
ax3.plot(time*10**-6, M_water_mo-M_water_atm/TO, label='magma ocean', color=cmap(120))
ax3.plot(time*10**-6, M_water_sol, label='solid', color=cmap(220))

ax3.plot(time_sch*10**-6, M_water_atm_sch/TO, linestyle='--', color=cmap(0))
ax3.plot(time_sch*10**-6, M_water_mo_sch-M_water_atm_sch/TO, linestyle='--', color=cmap(120))
ax3.plot(time_sch*10**-6, M_water_sol_sch, linestyle='--', color=cmap(220))

ax3.plot(time_4*10**-6, M_water_atm_4/TO, linestyle=':', color=cmap(0))
ax3.plot(time_4*10**-6, M_water_mo_4-M_water_atm_4/TO, linestyle=':', color=cmap(120))
ax3.plot(time_4*10**-6, M_water_sol_4, linestyle=':', color=cmap(220))

ax3.set_ylim([0.01,100])
ax3.legend(loc='best', frameon=True)
ax3.set_ylabel('Water Mass (TO)')
ax3.set_yscale('log')

ax4 = fig.add_subplot(234, sharex=ax1)
ax4.plot(time*10**-6, Press_H2O, label='$H_2O$', color=cmap(0))
ax4.plot(time*10**-6, Press_O, label='$O$', color=cmap(220))

ax4.plot(time_sch*10**-6, Press_H2O_sch, linestyle='--', color=cmap(0))
ax4.plot(time_sch*10**-6, Press_O_sch, linestyle='--', color=cmap(220))

ax4.plot(time_4*10**-6, Press_H2O_4, linestyle=':', color=cmap(0))
ax4.plot(time_4*10**-6, Press_O_4, linestyle=':', color=cmap(220))

ax4.legend(loc='best', frameon=True)
ax4.set_ylabel('Atmospheric pressure (bar)')
ax4.set_yscale('log')
ax4.set_xlabel('Time (Myrs)')
ax4.set_ylim([0.1,10000])

ax5 = fig.add_subplot(235, sharex=ax1)
ax5.plot(time*10**-6, Frac_H2O, color=cmap(0), label='$H_2O$')
ax5.plot(time*10**-6, Frac_Fe2O3, color=cmap(220), label='$Fe_2O_3$')

ax5.plot(time_sch*10**-6, Frac_H2O_sch, linestyle='--', color=cmap(0))
ax5.plot(time_sch*10**-6, Frac_Fe2O3_sch, linestyle='--', color=cmap(220))

ax5.plot(time_4*10**-6, Frac_H2O_4, linestyle=':', color=cmap(0))
ax5.plot(time_4*10**-6, Frac_Fe2O3_4, linestyle=':', color=cmap(220))

ax5.legend(loc='best', frameon=True)
ax5.set_ylabel('Mass frac in magma ocean')
# ax5.set_yscale('log')
# ax5.set_ylim([1e16,5e21])
ax5.set_xlabel('Time (Myrs)')
ax5.set_ylim([0.02,0.08])
# ax5.set_xlim([10,28])

# ax6 = fig.add_subplot(236, sharex=ax1)
# ax6.plot(time*10**-6, M_O_atm, label='atmosphere', color=cmap(0))
# ax6.plot(time*10**-6, M_O_mo-M_O_atm, label='magma ocean', color=cmap(120))
# ax6.plot(time*10**-6, M_O_sol, label='solid', color=cmap(220))
#
# ax6.plot(time_sch*10**-6, M_O_atm_sch, linestyle='--', color=cmap(0))
# ax6.plot(time_sch*10**-6, M_O_mo_sch-M_O_atm_sch, linestyle='--', color=cmap(120))
# ax6.plot(time_sch*10**-6, M_O_sol_sch, linestyle='--', color=cmap(220))
#
# ax6.legend(loc='best', frameon=True)
# ax6.set_ylabel('Oxygen Mass (kg)')
# ax6.set_yscale('log')
# ax6.set_ylim([1e16,1e21])

# ax6 = fig.add_subplot(236, sharex=ax1)
# ax6.plot(time*10**-6, NetFluxAtmo, color=cmap(0))
#
# ax6.plot(time_sch*10**-6, NetFluxAtmo_sch, linestyle='--', color=cmap(0))
#
# ax6.set_ylabel('Atmospheric net flux ($W/m^2$)')
# ax6.set_yscale('log')

ax6 = fig.add_subplot(236, sharex=ax1)
ax6.plot(time*10**-6, RadioHeat, color=cmap(0), label='Radiogenic')
ax6.plot(time_sch*10**-6, RadioHeat_sch, linestyle='--', color=cmap(0))
ax6.plot(time_4*10**-6, RadioHeat_4, linestyle=':', color=cmap(0))

ax6.legend(loc='best', frameon=True)
ax6.set_ylabel('Mantle Heating (TW)')
# ax6.set_yscale('log')
ax6.set_xlabel('Time (Myrs)')
ax6.set_ylim([50,300])

# ax7 = fig.add_subplot(337, sharex=ax1)
# ax7.plot(time*10**-6, N_H_space, label='space')
# ax7.plot(time*10**-6, N_H_atm, label='atmosphere')
# ax7.plot(time*10**-6, N_H_mo, label='magma ocean')
# ax7.plot(time*10**-6, N_H_sol, label='solid')
# ax7.plot(time*10**-6, N_H_tot, label='total', linestyle=':')
# ax7.legend(loc='best', frameon=True)
# ax7.set_xlabel('Time (Myrs)')
# ax7.set_ylabel('Number of H atoms ($10^{45}$)')
# ax7.set_yscale('log')
# ax7.set_ylim([1e-1,1e5])
#
# ax8 = fig.add_subplot(338, sharex=ax1)
# ax8.plot(time*10**-6, N_O_space, label='space')
# ax8.plot(time*10**-6, N_O_atm, label='atmosphere')
# ax8.plot(time*10**-6, N_O_mo, label='magma ocean')
# ax8.plot(time*10**-6, N_O_sol, label='solid')
# ax8.plot(time*10**-6, N_O_tot, label='total', linestyle=':')
# ax8.legend(loc='best', frameon=True)
# ax8.set_xlabel('Time (Myrs)')
# ax8.set_ylabel('Number of O atoms ($10^{45}$)')
# ax8.set_yscale('log')
# ax8.set_ylim([1e-1,1e5])
#
# ax9 = fig.add_subplot(339, sharex=ax1)
# ax9.plot(time*10**-6, N_H_tot/N_H_tot[0]-1, label='H')
# ax9.plot(time*10**-6, N_O_tot/N_O_tot[0]-1, label='O', linestyle='--')
# ax9.plot(time*10**-6, (N_O_tot/N_O_tot[0]-1)-(N_H_tot/N_H_tot[0]-1), label='O-H', linestyle=':')
# ax9.legend(loc='best', frameon=True)
# ax9.set_xlabel('Time (Myrs)')
# ax9.set_ylabel('Ratio of atoms lost')
# ax8.set_yscale('log')
#
# ax9 = fig.add_subplot(339, sharex=ax1)
# ax9.plot(time*10**-6, N_H_tot-N_H_tot[0], label='H')
# ax9.plot(time*10**-6, N_O_tot-N_O_tot[0], label='O')
# ax9.plot(time*10**-6, (N_O_tot-N_O_tot[0])-(N_H_tot-N_H_tot[0])/2., label='O-H/2', linestyle='--')
# ax9.legend(loc='best', frameon=True)
# ax9.set_xlabel('Time (Myrs)')
# ax9.set_ylabel('Number of atoms lost ($10^{45}$)')

plt.show()

# plt.figure()
# plt.title('GJ1132b: $M_{water}^{ini} = $'+str(M_water_mo[0])+' TO', fontsize=18, fontweight='bold')
# # plt.plot(time*10**-6, Tpot, label='T_p', linewidth=3.0)
# # plt.plot(time*10**-6, Tsurf, label='T_surf', linestyle='--', linewidth=3.0)
# plt.plot(time*10**-6, Tpot, label='T_p', linewidth=3.0, color=cmap(0))
# plt.plot(time*10**-6, Tsurf, label='T_surf', linestyle='--', linewidth=3.0, color=cmap(0))
# plt.xlim([1e-6,time[i_end+1]*1e-6])
# # ax1.set_ylim([0,4100])
# plt.legend(loc='best', frameon=True, fontsize=16)
# plt.xlabel('Time (Myrs)', fontsize=18, fontweight='bold')
# plt.ylabel('Temperature (K)', fontsize=18, fontweight='bold')
# if log_plot == 1:
#     plt.xscale('log')
# plt.xlim([1e-6,time[i_end+1]*1e-6])
# plt.tick_params(labelsize=14)
# plt.show()

# r_core = np.zeros(n_time)+r_c
#
# if individual == 1:
#     ## plot individual figures
#
#     plt.figure()
#     plt.title('GJ1132b: $M_{water}^{ini} = $'+str(Initial_water)+' TO', fontsize=18, fontweight='bold')
#     # plt.plot(time*10**-6, Tpot, label='T_p', linewidth=3.0)
#     # plt.plot(time*10**-6, Tsurf, label='T_surf', linestyle='--', linewidth=3.0)
#     plt.plot(time*10**-6, Tpot, label='T_p', linewidth=3.0, color=cmap(0))
#     plt.plot(time*10**-6, Tsurf, label='T_surf', linestyle='--', linewidth=3.0, color=cmap(0))
#     plt.xlim([1e-6,time[i_end+1]*1e-6])
#     # ax1.set_ylim([0,4100])
#     plt.legend(loc='best', frameon=True, fontsize=16)
#     plt.xlabel('Time (Myrs)', fontsize=18, fontweight='bold')
#     plt.ylabel('Temperature (K)', fontsize=18, fontweight='bold')
#     if log_plot == 1:
#         plt.xscale('log')
#     plt.xlim([1e-6,time[i_end+1]*1e-6])
#     plt.tick_params(labelsize=14)
#     plt.show()
#
#     plt.figure()
#     plt.title('GJ1132b: $M_{water}^{ini} = $'+str(Initial_water)+' TO', fontsize=18, fontweight='bold')
#     # plt.plot(time*10**-6, r_sol/r_p, label='$r_s$', color=cmap(100), linewidth=3.0)
#     # plt.plot(time*10**-6, r_core/r_p, label='$r_c$', color=cmap(220), linestyle='--', linewidth=3.0)
#     plt.plot(time*10**-6, r_sol, label='$r_s$', color=cmap(0), linewidth=3.0)
#     # plt.plot(time*10**-6, r_core, label='$r_c$', color=cmap(0), linestyle='--', linewidth=3.0)
#     plt.ylim([0.5,1])
#     plt.legend(loc='best', frameon=True, fontsize=16)
#     plt.xlabel('Time (Myrs)', fontsize=18, fontweight='bold')
#     plt.ylabel('Solidification radius ($r_p$)', fontsize=18, fontweight='bold')
#     if log_plot == 1:
#         plt.xscale('log')
#     plt.xlim([1e-6,time[i_end+1]*1e-6])
#     plt.tick_params(labelsize=14)
#     plt.show()
#
#     plt.figure()
#     plt.title('GJ1132b: $M_{water}^{ini} = $'+str(Initial_water)+' TO', fontsize=18, fontweight='bold')
#     # plt.plot(time*10**-6, q_m, label='Mantle heat flux', linewidth=3.0)
#     # plt.plot(time*10**-6, Flux_OLR, label='Outgoing longwave radiation', linewidth=3.0)
#     # plt.plot(time*10**-6, Flux_BOL, label='Absorbed stellar radiation', linewidth=3.0)
#     # plt.plot(time*10**-6, Flux_XUV, label='Absorbed stellar XUV flux', linewidth=3.0)
#     plt.plot(time*10**-6, q_m, label='Mantle heat flux', linewidth=3.0, color=cmap(0))
#     plt.plot(time*10**-6, OLR, label='Net Flux Atmosphere', linewidth=3.0, linestyle='--', color=cmap(0))
#     # plt.plot(time*10**-6, ASR, label='Absorbed stellar radiation', linewidth=3.0, linestyle='-.', color=cmap(0))
#     plt.plot(time*10**-6, XUV, label='Absorbed stellar XUV flux', linewidth=3.0, linestyle=':', color=cmap(0))
#     plt.legend(loc='best', frameon=True, fontsize=16)
#     plt.xlabel('Time (Myrs)', fontsize=18, fontweight='bold')
#     plt.ylabel('Flux ($W/m^2$)', fontsize=18, fontweight='bold')
#     plt.yscale('log')
#     if log_plot == 1:
#         plt.xscale('log')
#     plt.xlim([1e-6,time[i_end+1]*1e-6])
#     plt.tick_params(labelsize=14)
#     plt.show()
#
#     plt.figure()
#     plt.title('GJ1132b: $M_{water}^{ini} = $'+str(Initial_water)+' TO', fontsize=18, fontweight='bold')
#     # plt.plot(time*10**-6, M_water_mo/M_water_mo[0], label='magma ocean + atm', linewidth=3.0)
#     # plt.plot(time*10**-6, M_water_atm/M_water_mo[0], label='atmosphere', linewidth=3.0)
#     # plt.plot(time*10**-6, M_water_sol/M_water_mo[0], label='solid', linewidth=3.0)
#     plt.plot(time*10**-6, M_water_mo/M_water_mo[0], label='magma ocean + atm', linewidth=3.0, color=cmap(0))
#     plt.plot(time*10**-6, M_water_atm/M_water_mo[0], label='atmosphere', linewidth=3.0, linestyle='--', color=cmap(0))
#     plt.plot(time*10**-6, M_water_sol/M_water_mo[0], label='solid', linewidth=3.0, linestyle=':', color=cmap(0))
#     plt.ylim([1e-3,1.05])
#     plt.legend(loc='best', frameon=True, fontsize=16)
#     plt.xlabel('Time (Myrs)', fontsize=18, fontweight='bold')
#     plt.ylabel('Water Mass Fraction', fontsize=18, fontweight='bold')
#     plt.yscale('log')
#     if log_plot == 1:
#         plt.xscale('log')
#     plt.xlim([1e-6,time[i_end+1]*1e-6])
#     plt.tick_params(labelsize=14)
#     plt.show()
#
#     plt.figure()
#     plt.title('GJ1132b: $M_{water}^{ini} = $'+str(Initial_water)+' TO', fontsize=18, fontweight='bold')
#     plt.plot(time*10**-6, Press_atm, label='total pressure', linewidth=3.0, color=cmap(0))
#     plt.plot(time*10**-6, Press_H2O, label='partial pressure H2O', linestyle='--', linewidth=3.0, color=cmap(0))
#     plt.plot(time*10**-6, Press_O, label='partial pressure O', linestyle=':', linewidth=3.0, color=cmap(0))
#     plt.legend(loc='best', frameon=True, fontsize=16)
#     plt.xlabel('Time (Myrs)', fontsize=18, fontweight='bold')
#     plt.ylabel('Atmospheric pressure (bar)', fontsize=18, fontweight='bold')
#     if log_plot == 1:
#         plt.xscale('log')
#     plt.xlim([1e-6,time[i_end+1]*1e-6])
#     plt.yscale('log')
#     plt.tick_params(labelsize=14)
#     plt.show()
#
# else:
#     ## plot multiple figures
