from benchmark import Benchmark, benchmark 
import astropy.units as u 
import pytest 
 
@benchmark( 
   { 
       "log.initial.system.Age": {"value": 6.311520e+13, "unit": u.sec}, 
       "log.initial.system.Time": {"value": 0.000000, "unit": u.sec}, 
       "log.initial.system.TotAngMom": {"value": 5.357909e+43, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.system.TotEnergy": {"value": -1.192378e+41, "unit": u.erg}, 
       "log.initial.system.PotEnergy": {"value": -2.556201e+39, "unit": u.Joule}, 
       "log.initial.system.KinEnergy": {"value": 4.054947e+37, "unit": u.Joule}, 
       "log.initial.system.DeltaTime": {"value": 0.000000, "unit": u.sec}, 
       "log.initial.a.Mass": {"value": 1.988416e+29, "unit": u.kg}, 
       "log.initial.a.Radius": {"value": 97.114438, "unit": u.Rearth}, 
       "log.initial.a.RadGyra": {"value": 0.448345}, 
       "log.initial.a.RotAngMom": {"value": 1.115191e+42, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.a.RotVel": {"value": 4.504445e+04, "unit": u.m / u.sec}, 
       "log.initial.a.BodyType": {"value": 0.000000}, 
       "log.initial.a.RotRate": {"value": 7.272205e-05, "unit": 1 / u.sec}, 
       "log.initial.a.RotPer": {"value": 1.000000, "unit": u.day}, 
       "log.initial.a.Density": {"value": 199.752981, "unit": u.kg / u.m ** 3}, 
       "log.initial.a.HZLimitDryRunaway": {"value": 3.200490e+10, "unit": u.m}, 
       "log.initial.a.HZLimRecVenus": {"value": 1.360577e+11, "unit": u.m}, 
       "log.initial.a.HZLimRunaway": {"value": 1.790817e+11, "unit": u.m}, 
       "log.initial.a.HZLimMoistGreenhouse": {"value": 1.800268e+11, "unit": u.m}, 
       "log.initial.a.HZLimMaxGreenhouse": {"value": 3.452485e+11, "unit": u.m}, 
       "log.initial.a.HZLimEarlyMars": {"value": 3.765288e+11, "unit": u.m}, 
       "log.initial.a.Instellation": {"value": -1.000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.a.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m}, 
       "log.initial.a.LXUVTot": {"value": 2.136736e+22, "unit": u.kg / u.sec ** 3}, 
       "log.initial.a.LostEnergy": {"value": 5.562685e-309, "unit": u.Joule}, 
       "log.initial.a.LostAngMom": {"value": 5.562685e-309, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.a.Luminosity": {"value": 0.055557, "unit": u.LSUN}, 
       "log.initial.a.LXUVStellar": {"value": 2.136736e+22, "unit": u.W}, 
       "log.initial.a.Temperature": {"value": 2971.232396, "unit": u.K}, 
       "log.initial.a.LXUVFrac": {"value": 0.001000}, 
       "log.initial.a.RossbyNumber": {"value": 0.014575}, 
       "log.initial.a.DRotPerDtStellar": {"value": 4.420158e-10}, 
       "log.initial.b.Mass": {"value": 1.988416e+30, "unit": u.kg}, 
       "log.initial.b.Radius": {"value": 209.259428, "unit": u.Rearth}, 
       "log.initial.b.RadGyra": {"value": 0.451302}, 
       "log.initial.b.RotAngMom": {"value": 5.246390e+43, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.b.RotVel": {"value": 9.706049e+04, "unit": u.m / u.sec}, 
       "log.initial.b.BodyType": {"value": 0.000000}, 
       "log.initial.b.RotRate": {"value": 7.272205e-05, "unit": 1 / u.sec}, 
       "log.initial.b.RotPer": {"value": 1.000000, "unit": u.day}, 
       "log.initial.b.Density": {"value": 199.659310, "unit": u.kg / u.m ** 3}, 
       "log.initial.b.HZLimitDryRunaway": {"value": 3.200490e+10, "unit": u.m}, 
       "log.initial.b.HZLimRecVenus": {"value": 1.360577e+11, "unit": u.m}, 
       "log.initial.b.HZLimRunaway": {"value": 1.790817e+11, "unit": u.m}, 
       "log.initial.b.HZLimMoistGreenhouse": {"value": 1.800268e+11, "unit": u.m}, 
       "log.initial.b.HZLimMaxGreenhouse": {"value": 3.452485e+11, "unit": u.m}, 
       "log.initial.b.HZLimEarlyMars": {"value": 3.765288e+11, "unit": u.m}, 
       "log.initial.b.Instellation": {"value": 75.978415, "unit": u.kg / u.sec ** 3}, 
       "log.initial.b.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m}, 
       "log.initial.b.LXUVTot": {"value": 4.556110e+23, "unit": u.kg / u.sec ** 3}, 
       "log.initial.b.LostEnergy": {"value": 5.562685e-309, "unit": u.Joule}, 
       "log.initial.b.LostAngMom": {"value": 5.562685e-309, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.b.Luminosity": {"value": 1.184636, "unit": u.LSUN}, 
       "log.initial.b.LXUVStellar": {"value": 4.556110e+23, "unit": u.W}, 
       "log.initial.b.Temperature": {"value": 4349.796199, "unit": u.K}, 
       "log.initial.b.LXUVFrac": {"value": 0.001000}, 
       "log.initial.b.RossbyNumber": {"value": 0.029572}, 
       "log.initial.b.DRotPerDtStellar": {"value": -4.686689e-10}, 
       "log.final.system.Age": {"value": 3.218875e+15, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.system.Time": {"value": 3.155760e+15, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.system.TotAngMom": {"value": 5.362479e+43, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.system.TotEnergy": {"value": -1.191823e+41, "unit": u.erg, "rtol": 1e-4}, 
       "log.final.system.PotEnergy": {"value": -1.234473e+40, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.system.KinEnergy": {"value": 2.087504e+37, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.a.Mass": {"value": 1.988416e+29, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.a.Radius": {"value": 20.109309, "unit": u.Rearth, "rtol": 1e-4}, 
       "log.final.a.RadGyra": {"value": 0.464900, "rtol": 1e-4}, 
       "log.final.a.RotAngMom": {"value": 1.718027e+41, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.a.RotVel": {"value": 3.116850e+04, "unit": u.m / u.sec, "rtol": 1e-4}, 
       "log.final.a.BodyType": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.a.RotRate": {"value": 0.000243, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.a.RotPer": {"value": 0.299253, "unit": u.day, "rtol": 1e-4}, 
       "log.final.a.Density": {"value": 2.249848e+04, "unit": u.kg / u.m ** 3, "rtol": 1e-4}, 
       "log.final.a.HZLimitDryRunaway": {"value": 6.713550e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.a.HZLimRecVenus": {"value": 1.015610e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.a.HZLimRunaway": {"value": 1.336980e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.a.HZLimMoistGreenhouse": {"value": 1.343820e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.a.HZLimMaxGreenhouse": {"value": 2.575080e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.a.HZLimEarlyMars": {"value": 2.808387e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.a.Instellation": {"value": -1.000000, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.a.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.a.LXUVTot": {"value": 9.175805e+20, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.a.LostEnergy": {"value": 9.808345e+39, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.a.LostAngMom": {"value": 9.438647e+41, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.a.Luminosity": {"value": 0.002445, "unit": u.LSUN, "rtol": 1e-4}, 
       "log.final.a.LXUVStellar": {"value": 9.175805e+20, "unit": u.W, "rtol": 1e-4}, 
       "log.final.a.Temperature": {"value": 2992.330594, "unit": u.K, "rtol": 1e-4}, 
       "log.final.a.LXUVFrac": {"value": 0.000976, "rtol": 1e-4}, 
       "log.final.a.RossbyNumber": {"value": 0.004409, "rtol": 1e-4}, 
       "log.final.a.DRotPerDtStellar": {"value": -3.016808e-12, "rtol": 1e-4}, 
       "log.final.b.Mass": {"value": 1.988416e+30, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.b.Radius": {"value": 98.456438, "unit": u.Rearth, "rtol": 1e-4}, 
       "log.final.b.RadGyra": {"value": 0.298250, "rtol": 1e-4}, 
       "log.final.b.RotAngMom": {"value": 1.649161e+42, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.b.RotVel": {"value": 1.484774e+04, "unit": u.m / u.sec, "rtol": 1e-4}, 
       "log.final.b.BodyType": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.b.RotRate": {"value": 2.364421e-05, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.b.RotPer": {"value": 3.075681, "unit": u.day, "rtol": 1e-4}, 
       "log.final.b.Density": {"value": 1916.956727, "unit": u.kg / u.m ** 3, "rtol": 1e-4}, 
       "log.final.b.HZLimitDryRunaway": {"value": 6.713550e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimRecVenus": {"value": 1.015610e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimRunaway": {"value": 1.336980e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimMoistGreenhouse": {"value": 1.343820e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimMaxGreenhouse": {"value": 2.575080e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimEarlyMars": {"value": 2.808387e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.Instellation": {"value": 3.343195, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.b.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.LXUVTot": {"value": 2.586460e+23, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.b.LostEnergy": {"value": 1.354496e+41, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.b.LostAngMom": {"value": 5.085996e+43, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.b.Luminosity": {"value": 0.689088, "unit": u.LSUN, "rtol": 1e-4}, 
       "log.final.b.LXUVStellar": {"value": 2.586460e+23, "unit": u.W, "rtol": 1e-4}, 
       "log.final.b.Temperature": {"value": 5539.189837, "unit": u.K, "rtol": 1e-4}, 
       "log.final.b.LXUVFrac": {"value": 0.000976, "rtol": 1e-4}, 
       "log.final.b.RossbyNumber": {"value": 0.187127, "rtol": 1e-4}, 
       "log.final.b.DRotPerDtStellar": {"value": 2.479932e-10, "rtol": 1e-4}, 
   } 
)
class Test_StellarEvol(Benchmark): 
   pass 
