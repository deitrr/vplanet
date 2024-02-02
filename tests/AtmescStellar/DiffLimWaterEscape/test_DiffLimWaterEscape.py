from benchmark import Benchmark, benchmark 
import astropy.units as u 
import pytest 
 
@benchmark( 
   { 
       "log.initial.system.Age": {"value": 3.155760e+13, "unit": u.sec}, 
       "log.initial.system.Time": {"value": 0.000000, "unit": u.sec}, 
       "log.initial.system.TotAngMom": {"value": 1.151432e+42, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.system.TotEnergy": {"value": -1.877843e+39, "unit": u.Joule}, 
       "log.initial.system.PotEnergy": {"value": -1.919599e+39, "unit": u.Joule}, 
       "log.initial.system.KinEnergy": {"value": 4.178988e+37, "unit": u.Joule}, 
       "log.initial.system.DeltaTime": {"value": 0.000000, "unit": u.sec}, 
       "log.initial.star.Mass": {"value": 1.789574e+29, "unit": u.kg}, 
       "log.initial.star.Radius": {"value": 6.681051e+08, "unit": u.m}, 
       "log.initial.star.RadGyra": {"value": 0.444800}, 
       "log.initial.star.RotAngMom": {"value": 1.149304e+42, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.star.RotVel": {"value": 4.858597e+04, "unit": u.m / u.sec}, 
       "log.initial.star.BodyType": {"value": 0.000000}, 
       "log.initial.star.RotRate": {"value": 7.272205e-05, "unit": 1 / u.sec}, 
       "log.initial.star.RotPer": {"value": 8.640000e+04, "unit": u.sec}, 
       "log.initial.star.Density": {"value": 143.260634, "unit": u.kg / u.m ** 3}, 
       "log.initial.star.HZLimitDryRunaway": {"value": 3.305062e+10, "unit": u.m}, 
       "log.initial.star.HZLimRecVenus": {"value": 2.977014e+10, "unit": u.m}, 
       "log.initial.star.HZLimRunaway": {"value": 3.916594e+10, "unit": u.m}, 
       "log.initial.star.HZLimMoistGreenhouse": {"value": 3.939084e+10, "unit": u.m}, 
       "log.initial.star.HZLimMaxGreenhouse": {"value": 7.572440e+10, "unit": u.m}, 
       "log.initial.star.HZLimEarlyMars": {"value": 8.258537e+10, "unit": u.m}, 
       "log.initial.star.Instellation": {"value": -1.000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.star.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m}, 
       "log.initial.star.LXUVTot": {"value": 5.924722e-05, "unit": u.LSUN}, 
       "log.initial.star.LostEnergy": {"value": 5.562685e-309, "unit": u.Joule}, 
       "log.initial.star.LostAngMom": {"value": 5.562685e-309, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.star.Luminosity": {"value": 0.059247, "unit": u.LSUN}, 
       "log.initial.star.LXUVStellar": {"value": 2.278648e+22, "unit": u.W}, 
       "log.initial.star.Temperature": {"value": 2907.334487, "unit": u.K}, 
       "log.initial.star.LXUVFrac": {"value": 0.001000}, 
       "log.initial.star.RossbyNumber": {"value": 0.014106}, 
       "log.initial.star.DRotPerDtStellar": {"value": 2.355449e-09}, 
       "log.initial.b.Mass": {"value": 8.223700e+24, "unit": u.kg}, 
       "log.initial.b.Radius": {"value": 7.124338e+06, "unit": u.m}, 
       "log.initial.b.RadGyra": {"value": 0.500000}, 
       "log.initial.b.BodyType": {"value": 0.000000}, 
       "log.initial.b.Density": {"value": 5429.316562, "unit": u.kg / u.m ** 3}, 
       "log.initial.b.HZLimitDryRunaway": {"value": 3.305089e+10, "unit": u.m}, 
       "log.initial.b.HZLimRecVenus": {"value": 2.977014e+10, "unit": u.m}, 
       "log.initial.b.HZLimRunaway": {"value": 3.916594e+10, "unit": u.m}, 
       "log.initial.b.HZLimMoistGreenhouse": {"value": 3.939084e+10, "unit": u.m}, 
       "log.initial.b.HZLimMaxGreenhouse": {"value": 7.572440e+10, "unit": u.m}, 
       "log.initial.b.HZLimEarlyMars": {"value": 8.258537e+10, "unit": u.m}, 
       "log.initial.b.Instellation": {"value": 6.076083e+05, "unit": u.kg / u.sec ** 3}, 
       "log.initial.b.MeanMotion": {"value": 4.813397e-05, "unit": 1 / u.sec}, 
       "log.initial.b.OrbPeriod": {"value": 1.305354e+05, "unit": u.sec}, 
       "log.initial.b.SemiMajorAxis": {"value": 1.727522e+09, "unit": u.m}, 
       "log.initial.b.LXUVTot": {"value": -1.000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.b.SurfWaterMass": {"value": 10.000000, "unit": u.TO}, 
       "log.initial.b.EnvelopeMass": {"value": 0.000000, "unit": u.kg}, 
       "log.initial.b.OxygenMass": {"value": 0.000000, "unit": u.bar}, 
       "log.initial.b.RGLimit": {"value": 3.770574e+10, "unit": u.m}, 
       "log.initial.b.XO": {"value": 0.333333}, 
       "log.initial.b.EtaO": {"value": 0.902643}, 
       "log.initial.b.PlanetRadius": {"value": 7.124338e+06, "unit": u.m}, 
       "log.initial.b.OxygenMantleMass": {"value": 0.000000, "unit": u.kg}, 
       "log.initial.b.RadXUV": {"value": -1.000000, "unit": u.m}, 
       "log.initial.b.RadSolid": {"value": -1.000000, "unit": u.m}, 
       "log.initial.b.PresXUV": {"value": 5.000000}, 
       "log.initial.b.ScaleHeight": {"value": -1.000000, "unit": u.m}, 
       "log.initial.b.ThermTemp": {"value": 400.000000, "unit": u.K}, 
       "log.initial.b.AtmGasConst": {"value": 4124.000000}, 
       "log.initial.b.PresSurf": {"value": -1.000000, "unit": u.Pa}, 
       "log.initial.b.DEnvMassDt": {"value": 0.000000, "unit": u.kg / u.sec}, 
       "log.initial.b.FXUV": {"value": 607.608327, "unit": u.W / u.m ** 2}, 
       "log.initial.b.AtmXAbsEffH2O": {"value": 0.010000}, 
       "log.initial.b.RocheRadius": {"value": 4.290313e+07, "unit": u.m}, 
       "log.initial.b.BondiRadius": {"value": 3.107839e+08, "unit": u.m}, 
       "log.initial.b.HEscapeRegime": {"value": 8.000000}, 
       "log.initial.b.RRCriticalFlux": {"value": 41.454587, "unit": u.W / u.m ** 2}, 
       "log.initial.b.CrossoverMass": {"value": 2.575027e-25, "unit": u.kg}, 
       "log.initial.b.WaterEscapeRegime": {"value": 3.000000}, 
       "log.initial.b.FXUVCRITDRAG": {"value": 7.247839, "unit": u.W / u.m ** 2}, 
       "log.initial.b.HREFFLUX": {"value": 1.178796e+19, "unit": 1 / u.m ** 2 / u.sec}, 
       "log.initial.b.XO2": {"value": 0.000000}, 
       "log.initial.b.XH2O": {"value": 1.000000}, 
       "log.initial.b.HDiffFlux": {"value": 1.395966e+17, "unit": 1 / u.m ** 2 / u.sec}, 
       "log.initial.b.HRefODragMod": {"value": 0.121638}, 
       "log.initial.b.KTide": {"value": 0.753205}, 
       "log.initial.b.RGDuration": {"value": 0.00000e+00, "unit": u.yr}, 
       "log.initial.e.Mass": {"value": 4.138725e+24, "unit": u.kg}, 
       "log.initial.e.Radius": {"value": 5.880608e+06, "unit": u.m}, 
       "log.initial.e.RadGyra": {"value": 0.500000}, 
       "log.initial.e.BodyType": {"value": 0.000000}, 
       "log.initial.e.Density": {"value": 4858.600773, "unit": u.kg / u.m ** 3}, 
       "log.initial.e.HZLimitDryRunaway": {"value": 3.305129e+10, "unit": u.m}, 
       "log.initial.e.HZLimRecVenus": {"value": 2.977014e+10, "unit": u.m}, 
       "log.initial.e.HZLimRunaway": {"value": 3.916594e+10, "unit": u.m}, 
       "log.initial.e.HZLimMoistGreenhouse": {"value": 3.939084e+10, "unit": u.m}, 
       "log.initial.e.HZLimMaxGreenhouse": {"value": 7.572440e+10, "unit": u.m}, 
       "log.initial.e.HZLimEarlyMars": {"value": 8.258537e+10, "unit": u.m}, 
       "log.initial.e.Instellation": {"value": 9.448996e+04, "unit": u.kg / u.sec ** 3}, 
       "log.initial.e.MeanMotion": {"value": 1.191967e-05, "unit": 1 / u.sec}, 
       "log.initial.e.OrbPeriod": {"value": 5.271275e+05, "unit": u.sec}, 
       "log.initial.e.SemiMajorAxis": {"value": 4.380718e+09, "unit": u.m}, 
       "log.initial.e.LXUVTot": {"value": -1.000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.e.SurfWaterMass": {"value": 10.000000, "unit": u.TO}, 
       "log.initial.e.EnvelopeMass": {"value": 0.000000, "unit": u.kg}, 
       "log.initial.e.OxygenMass": {"value": 0.000000, "unit": u.bar}, 
       "log.initial.e.RGLimit": {"value": 3.832506e+10, "unit": u.m}, 
       "log.initial.e.XO": {"value": 0.333333}, 
       "log.initial.e.EtaO": {"value": 0.776474}, 
       "log.initial.e.PlanetRadius": {"value": 5.880608e+06, "unit": u.m}, 
       "log.initial.e.OxygenMantleMass": {"value": 0.000000, "unit": u.kg}, 
       "log.initial.e.RadXUV": {"value": -1.000000, "unit": u.m}, 
       "log.initial.e.RadSolid": {"value": -1.000000, "unit": u.m}, 
       "log.initial.e.PresXUV": {"value": 5.000000}, 
       "log.initial.e.ScaleHeight": {"value": -1.000000, "unit": u.m}, 
       "log.initial.e.ThermTemp": {"value": 400.000000, "unit": u.K}, 
       "log.initial.e.AtmGasConst": {"value": 4124.000000}, 
       "log.initial.e.PresSurf": {"value": -1.000000, "unit": u.Pa}, 
       "log.initial.e.DEnvMassDt": {"value": 0.000000, "unit": u.kg / u.sec}, 
       "log.initial.e.FXUV": {"value": 94.489956, "unit": u.W / u.m ** 2}, 
       "log.initial.e.AtmXAbsEffH2O": {"value": 0.011065}, 
       "log.initial.e.RocheRadius": {"value": 8.653872e+07, "unit": u.m}, 
       "log.initial.e.BondiRadius": {"value": 9.821934e+07, "unit": u.m}, 
       "log.initial.e.HEscapeRegime": {"value": 8.000000}, 
       "log.initial.e.RRCriticalFlux": {"value": 26.551159, "unit": u.W / u.m ** 2}, 
       "log.initial.e.CrossoverMass": {"value": 1.130931e-25, "unit": u.kg}, 
       "log.initial.e.WaterEscapeRegime": {"value": 3.000000}, 
       "log.initial.e.FXUVCRITDRAG": {"value": 2.949984, "unit": u.W / u.m ** 2}, 
       "log.initial.e.HREFFLUX": {"value": 3.326851e+18, "unit": 1 / u.m ** 2 / u.sec}, 
       "log.initial.e.XO2": {"value": 0.000000}, 
       "log.initial.e.XH2O": {"value": 1.000000}, 
       "log.initial.e.HDiffFlux": {"value": 1.031142e+17, "unit": 1 / u.m ** 2 / u.sec}, 
       "log.initial.e.HRefODragMod": {"value": 0.138662}, 
       "log.initial.e.KTide": {"value": 0.898227}, 
       "log.initial.e.RGDuration": {"value": 0.00000e+00, "unit": u.yr}, 
       "log.final.system.Age": {"value": 3.001128e+16, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.system.Time": {"value": 2.997972e+16, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.system.TotAngMom": {"value": 1.140605e+42, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.system.TotEnergy": {"value": -1.890023e+39, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.system.PotEnergy": {"value": -1.631379e+40, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.system.KinEnergy": {"value": 9.851807e+35, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.star.Mass": {"value": 1.789574e+29, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.star.Radius": {"value": 7.861410e+07, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.RadGyra": {"value": 0.465100, "rtol": 1e-4}, 
       "log.final.star.RotAngMom": {"value": 2.171174e+40, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.star.RotVel": {"value": 7134.305972, "unit": u.m / u.sec, "rtol": 1e-4}, 
       "log.final.star.BodyType": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.star.RotRate": {"value": 9.075097e-05, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.star.RotPer": {"value": 6.923546e+04, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.star.Density": {"value": 8.793460e+04, "unit": u.kg / u.m ** 3, "rtol": 1e-4}, 
       "log.final.star.HZLimitDryRunaway": {"value": 3.236261e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.HZLimRecVenus": {"value": 2.926995e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.HZLimRunaway": {"value": 3.845735e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.HZLimMoistGreenhouse": {"value": 3.872906e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.HZLimMaxGreenhouse": {"value": 7.519102e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.HZLimEarlyMars": {"value": 8.200580e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.Instellation": {"value": -1.000000, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.star.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.star.LXUVTot": {"value": 5.680619e-07, "unit": u.LSUN, "rtol": 1e-4}, 
       "log.final.star.LostEnergy": {"value": 1.442281e+40, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.star.LostAngMom": {"value": 1.116766e+42, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.star.Luminosity": {"value": 0.000568, "unit": u.LSUN, "rtol": 1e-4}, 
       "log.final.star.LXUVStellar": {"value": 2.184766e+20, "unit": u.W, "rtol": 1e-4}, 
       "log.final.star.Temperature": {"value": 2649.348933, "unit": u.K, "rtol": 1e-4}, 
       "log.final.star.LXUVFrac": {"value": 0.001000, "rtol": 1e-4}, 
       "log.final.star.RossbyNumber": {"value": 0.009905, "rtol": 1e-4}, 
       "log.final.star.DRotPerDtStellar": {"value": 1.447229e-12, "rtol": 1e-4}, 
       "log.final.b.Mass": {"value": 8.223700e+24, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.b.Radius": {"value": 7.124338e+06, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.RadGyra": {"value": 0.500000, "rtol": 1e-4}, 
       "log.final.b.BodyType": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.b.Density": {"value": 5429.316562, "unit": u.kg / u.m ** 3, "rtol": 1e-4}, 
       "log.final.b.HZLimitDryRunaway": {"value": 3.236287e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimRecVenus": {"value": 2.926995e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimRunaway": {"value": 3.845735e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimMoistGreenhouse": {"value": 3.872906e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimMaxGreenhouse": {"value": 7.519102e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HZLimEarlyMars": {"value": 8.200580e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.Instellation": {"value": 5825.743992, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.b.MeanMotion": {"value": 4.813397e-05, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.b.OrbPeriod": {"value": 1.305354e+05, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.b.SemiMajorAxis": {"value": 1.727522e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.LXUVTot": {"value": -1.000000, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.b.SurfWaterMass": {"value": 0.000000, "unit": u.TO, "rtol": 1e-4}, 
       "log.final.b.EnvelopeMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.b.OxygenMass": {"value": 1146.053603, "unit": u.bar, "rtol": 1e-4}, 
       "log.final.b.RGLimit": {"value": 3.704010e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.XO": {"value": 1.000000, "rtol": 1e-4}, 
       "log.final.b.EtaO": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.b.PlanetRadius": {"value": 7.124338e+06, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.OxygenMantleMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.b.RadXUV": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.RadSolid": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.PresXUV": {"value": 5.000000, "rtol": 1e-4}, 
       "log.final.b.ScaleHeight": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.ThermTemp": {"value": 400.000000, "unit": u.K, "rtol": 1e-4}, 
       "log.final.b.AtmGasConst": {"value": 4124.000000, "rtol": 1e-4}, 
       "log.final.b.PresSurf": {"value": -1.000000, "unit": u.Pa, "rtol": 1e-4}, 
       "log.final.b.DEnvMassDt": {"value": 0.000000, "unit": u.kg / u.sec, "rtol": 1e-4}, 
       "log.final.b.FXUV": {"value": 5.825744, "unit": u.W / u.m ** 2, "rtol": 1e-4}, 
       "log.final.b.AtmXAbsEffH2O": {"value": 0.041374, "rtol": 1e-4}, 
       "log.final.b.RocheRadius": {"value": 4.290313e+07, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.BondiRadius": {"value": 9.942293e+08, "unit": u.m, "rtol": 1e-4}, 
       "log.final.b.HEscapeRegime": {"value": 8.000000, "rtol": 1e-4}, 
       "log.final.b.RRCriticalFlux": {"value": 41.454587, "unit": u.W / u.m ** 2, "rtol": 1e-4}, 
       "log.final.b.CrossoverMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.b.WaterEscapeRegime": {"value": 8.000000, "rtol": 1e-4}, 
       "log.final.b.FXUVCRITDRAG": {"value": 0.000000, "unit": u.W / u.m ** 2, "rtol": 1e-4}, 
       "log.final.b.HREFFLUX": {"value": 4.676160e+17, "unit": 1 / u.m ** 2 / u.sec, "rtol": 1e-4}, 
       "log.final.b.XO2": {"value": 1.000000, "rtol": 1e-4}, 
       "log.final.b.XH2O": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.b.HDiffFlux": {"value": 0.000000, "unit": 1 / u.m ** 2 / u.sec, "rtol": 1e-4}, 
       "log.final.b.HRefODragMod": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.b.KTide": {"value": 0.753205, "rtol": 1e-4}, 
       "log.final.b.RGDuration": {"value": 0.00000e+00, "unit": u.yr, "rtol": 1e-4}, 
       "log.final.e.Mass": {"value": 4.138725e+24, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.e.Radius": {"value": 5.880608e+06, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.RadGyra": {"value": 0.500000, "rtol": 1e-4}, 
       "log.final.e.BodyType": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.e.Density": {"value": 4858.600773, "unit": u.kg / u.m ** 3, "rtol": 1e-4}, 
       "log.final.e.HZLimitDryRunaway": {"value": 3.236326e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.HZLimRecVenus": {"value": 2.926995e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.HZLimRunaway": {"value": 3.845735e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.HZLimMoistGreenhouse": {"value": 3.872906e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.HZLimMaxGreenhouse": {"value": 7.519102e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.HZLimEarlyMars": {"value": 8.200580e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.Instellation": {"value": 905.968971, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.e.MeanMotion": {"value": 1.191967e-05, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.e.OrbPeriod": {"value": 5.271275e+05, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.e.SemiMajorAxis": {"value": 4.380718e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.LXUVTot": {"value": -1.000000, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.e.SurfWaterMass": {"value": 3.053422, "unit": u.TO, "rtol": 1e-4}, 
       "log.final.e.EnvelopeMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.e.OxygenMass": {"value": 1031.243496, "unit": u.bar, "rtol": 1e-4}, 
       "log.final.e.RGLimit": {"value": 3.765018e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.XO": {"value": 0.554278, "rtol": 1e-4}, 
       "log.final.e.EtaO": {"value": 0.000000, "rtol": 1e-4}, 
       "log.final.e.PlanetRadius": {"value": 5.880608e+06, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.OxygenMantleMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.e.RadXUV": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.RadSolid": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.PresXUV": {"value": 5.000000, "rtol": 1e-4}, 
       "log.final.e.ScaleHeight": {"value": -1.000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.ThermTemp": {"value": 400.000000, "unit": u.K, "rtol": 1e-4}, 
       "log.final.e.AtmGasConst": {"value": 4124.000000, "rtol": 1e-4}, 
       "log.final.e.PresSurf": {"value": -1.000000, "unit": u.Pa, "rtol": 1e-4}, 
       "log.final.e.DEnvMassDt": {"value": 0.000000, "unit": u.kg / u.sec, "rtol": 1e-4}, 
       "log.final.e.FXUV": {"value": 0.905969, "unit": u.W / u.m ** 2, "rtol": 1e-4}, 
       "log.final.e.AtmXAbsEffH2O": {"value": 0.073282, "rtol": 1e-4}, 
       "log.final.e.RocheRadius": {"value": 8.653872e+07, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.BondiRadius": {"value": 3.142136e+08, "unit": u.m, "rtol": 1e-4}, 
       "log.final.e.HEscapeRegime": {"value": 8.000000, "rtol": 1e-4}, 
       "log.final.e.RRCriticalFlux": {"value": 26.551159, "unit": u.W / u.m ** 2, "rtol": 1e-4}, 
       "log.final.e.CrossoverMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.e.WaterEscapeRegime": {"value": 8.000000, "rtol": 1e-4}, 
       "log.final.e.FXUVCRITDRAG": {"value": 0.297807, "unit": u.W / u.m ** 2, "rtol": 1e-4}, 
       "log.final.e.HREFFLUX": {"value": 2.112523e+17, "unit": 1 / u.m ** 2 / u.sec, "rtol": 1e-4}, 
       "log.final.e.XO2": {"value": 0.426457, "rtol": 1e-4}, 
       "log.final.e.XH2O": {"value": 0.573543, "rtol": 1e-4}, 
       "log.final.e.HDiffFlux": {"value": 6.894047e+16, "unit": 1 / u.m ** 2 / u.sec, "rtol": 1e-4}, 
       "log.final.e.HRefODragMod": {"value": 0.003131, "rtol": 1e-4}, 
       "log.final.e.KTide": {"value": 0.898227, "rtol": 1e-4}, 
       "log.final.e.RGDuration": {"value": 3.80345e+08, "unit": u.yr, "rtol": 1e-4}, 
   } 
)
class Test_DiffLimWaterEscape(Benchmark): 
   pass 
