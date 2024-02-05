from benchmark import Benchmark, benchmark 
import astropy.units as u 
import pytest 
 
@benchmark( 
   { 
       "log.initial.system.Age": {"value": 0.000000, "unit": u.sec}, 
       "log.initial.system.Time": {"value": 0.000000, "unit": u.Gyr}, 
       "log.initial.system.TotAngMom": {"value": 4.416946e+33, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.system.TotEnergy": {"value": -2.237790e+32, "unit": u.Joule}, 
       "log.initial.system.PotEnergy": {"value": -2.239397e+32, "unit": u.Joule}, 
       "log.initial.system.KinEnergy": {"value": 1.606047e+29, "unit": u.Joule}, 
       "log.initial.system.DeltaTime": {"value": 0.000000, "unit": u.sec}, 
       "log.initial.earth.Mass": {"value": 5.972186e+24, "unit": u.kg}, 
       "log.initial.earth.Radius": {"value": 6.378100e+06, "unit": u.m}, 
       "log.initial.earth.RadGyra": {"value": 0.500000}, 
       "log.initial.earth.BodyType": {"value": 0.000000}, 
       "log.initial.earth.Density": {"value": 5495.038549, "unit": u.kg / u.m ** 3}, 
       "log.initial.earth.SurfEnFluxTotal": {"value": 0.278724, "unit": u.kg / u.sec ** 3}, 
       "log.initial.earth.HZLimitDryRunaway": {"value": -1.000000, "unit": u.m}, 
       "log.initial.earth.HZLimRecVenus": {"value": -1.000000}, 
       "log.initial.earth.HZLimRunaway": {"value": -1.000000}, 
       "log.initial.earth.HZLimMoistGreenhouse": {"value": -1.000000}, 
       "log.initial.earth.HZLimMaxGreenhouse": {"value": -1.000000}, 
       "log.initial.earth.HZLimEarlyMars": {"value": -1.000000}, 
       "log.initial.earth.Instellation": {"value": -1.000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.earth.D26AlPowerDt": {"value": -1.000000}, 
       "log.initial.earth.D26AlNumManDt": {"value": 0.000000, "unit": 1 / u.sec}, 
       "log.initial.earth.D40KPowerDt": {"value": -1.000000}, 
       "log.initial.earth.D40KNumManDt": {"value": -1.694597e+26, "unit": 1 / u.sec}, 
       "log.initial.earth.D232ThNumManDt": {"value": -9.534855e+23, "unit": 1 / u.sec}, 
       "log.initial.earth.D238UNumManDt": {"value": -1.408805e+24, "unit": 1 / u.sec}, 
       "log.initial.earth.D235UNumManDt": {"value": -3.089012e+24, "unit": 1 / u.sec}, 
       "log.initial.earth.RadPowerMan": {"value": 74.591573, "unit": u.TW}, 
       "log.initial.earth.RadPowerCore": {"value": 34.624713, "unit": u.TW}, 
       "log.initial.earth.RadPowerCrust": {"value": 33.267825, "unit": u.TW}, 
       "log.initial.earth.RadPowerTotal": {"value": 142.484111, "unit": u.TW}, 
       "log.initial.earth.SurfEnFluxRadTotal": {"value": 0.278724, "unit": u.kg / u.sec ** 3}, 
       "log.final.system.Age": {"value": 1.420092e+17, "unit": u.sec}, 
       "log.final.system.Time": {"value": 4.500000, "unit": u.Gyr}, 
       "log.final.system.TotAngMom": {"value": 4.416946e+33, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.final.system.TotEnergy": {"value": -2.237790e+32, "unit": u.Joule}, 
       "log.final.system.PotEnergy": {"value": -2.239397e+32, "unit": u.Joule}, 
       "log.final.system.KinEnergy": {"value": 1.606047e+29, "unit": u.Joule}, 
       "log.final.system.DeltaTime": {"value": 2.874680e+14, "unit": u.sec}, 
       "log.final.earth.Mass": {"value": 5.972186e+24, "unit": u.kg}, 
       "log.final.earth.Radius": {"value": 6.378100e+06, "unit": u.m}, 
       "log.final.earth.RadGyra": {"value": 0.500000}, 
       "log.final.earth.BodyType": {"value": 0.000000}, 
       "log.final.earth.Density": {"value": 5495.038549, "unit": u.kg / u.m ** 3}, 
       "log.final.earth.SurfEnFluxTotal": {"value": 0.047491, "unit": u.kg / u.sec ** 3}, 
       "log.final.earth.HZLimitDryRunaway": {"value": -1.000000, "unit": u.m}, 
       "log.final.earth.HZLimRecVenus": {"value": -1.000000}, 
       "log.final.earth.HZLimRunaway": {"value": -1.000000}, 
       "log.final.earth.HZLimMoistGreenhouse": {"value": -1.000000}, 
       "log.final.earth.HZLimMaxGreenhouse": {"value": -1.000000}, 
       "log.final.earth.HZLimEarlyMars": {"value": -1.000000}, 
       "log.final.earth.Instellation": {"value": -1.000000, "unit": u.kg / u.sec ** 3}, 
       "log.final.earth.D26AlPowerDt": {"value": -1.000000}, 
       "log.final.earth.D26AlNumManDt": {"value": 0.000000, "unit": 1 / u.sec}, 
       "log.final.earth.D40KPowerDt": {"value": -1.000000}, 
       "log.final.earth.D40KNumManDt": {"value": -1.425474e+25, "unit": 1 / u.sec}, 
       "log.final.earth.D232ThNumManDt": {"value": -7.630887e+23, "unit": 1 / u.sec}, 
       "log.final.earth.D238UNumManDt": {"value": -7.013429e+23, "unit": 1 / u.sec}, 
       "log.final.earth.D235UNumManDt": {"value": -3.671211e+22, "unit": 1 / u.sec}, 
       "log.final.earth.RadPowerMan": {"value": 14.306031, "unit": u.TW}, 
       "log.final.earth.RadPowerCore": {"value": 3.030301, "unit": u.TW}, 
       "log.final.earth.RadPowerCrust": {"value": 6.926660, "unit": u.TW}, 
       "log.final.earth.RadPowerTotal": {"value": 24.277332, "unit": u.TW}, 
       "log.final.earth.SurfEnFluxRadTotal": {"value": 0.047491, "unit": u.kg / u.sec ** 3}, 
   } 
)
class Test_InertEarth(Benchmark): 
   pass 
