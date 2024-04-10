import numpy as np
import matplotlib.pyplot as plt 

Time ,Latitude ,TempLat ,AlbedoLat ,AnnInsol ,FluxIn ,FluxOut ,IceMass ,IceHeight, DIceMassDt, IceFlow, BedrockH, IceAccum, IceAblate, TempLandL ,TempWaterL, AlbedoLand, AlbedoWater, PlanckAAvg = np.loadtxt('solarsys.Earth.Climate',unpack=True)
timef,two, three, four, temp, five, six, seven,eight  = np.loadtxt('solarsys.Earth.forward',unpack=True)

indi = np.arange(0,150150,150)
new_Planck = []
for i in indi:
    new_Planck.append(PlanckAAvg[i])
    


fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(new_Planck,temp)
ax.set_ylabel(f"Global Mean Temperature ($^{'o'}$C)")
ax.set_xlabel("PlanckAAvg")
fig.savefig("PlanckA_vs_global_mean_temp.png", dpi=300)

print(len(temp), len(new_Planck))