import astropy.units as u
from benchmark import Benchmark, benchmark


@benchmark(
    {
        "log.initial.system.Age": {"value": 1.577880e13, "unit": u.sec},
        "log.initial.system.Time": {"value": 0.000000, "unit": u.sec},
        "log.initial.system.TotAngMom": {
            "value": 8.296053e43,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.system.TotEnergy": {"value": -9.104482e40, "unit": u.Joule},
        "log.initial.system.PotEnergy": {"value": -9.406038e40, "unit": u.Joule},
        "log.initial.system.KinEnergy": {"value": 3.015563e39, "unit": u.Joule},
        "log.initial.system.DeltaTime": {"value": 0.000000, "unit": u.sec},
        "log.initial.Sun.Mass": {"value": 1.988416e30, "unit": u.kg},
        "log.initial.Sun.Radius": {"value": 263.919878, "unit": u.Rearth},
        "log.initial.Sun.RadGyra": {"value": 0.449900},
        "log.initial.Sun.RotAngMom": {
            "value": 8.293392e43,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.Sun.RotVel": {"value": 1.224136e05, "unit": u.m / u.sec},
        "log.initial.Sun.BodyType": {"value": 0.000000},
        "log.initial.Sun.RotRate": {"value": 7.272205e-05, "unit": 1 / u.sec},
        "log.initial.Sun.RotPer": {"value": 8.640000e04, "unit": u.sec},
        "log.initial.Sun.Density": {"value": 99.524124, "unit": u.kg / u.m**3},
        "log.initial.Sun.HZLimitDryRunaway": {"value": 1.887931e11, "unit": u.m},
        "log.initial.Sun.HZLimRecVenus": {"value": 1.640992e11, "unit": u.m},
        "log.initial.Sun.HZLimRunaway": {"value": 2.178515e11, "unit": u.m},
        "log.initial.Sun.HZLimMoistGreenhouse": {"value": 2.171294e11, "unit": u.m},
        "log.initial.Sun.HZLimMaxGreenhouse": {"value": 3.929566e11, "unit": u.m},
        "log.initial.Sun.HZLimEarlyMars": {"value": 4.286401e11, "unit": u.m},
        "log.initial.Sun.Instellation": {"value": -1.000000, "unit": u.kg / u.sec**3},
        "log.initial.Sun.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m},
        "log.initial.Sun.LXUVTot": {"value": 7.435159e23, "unit": u.kg / u.sec**3},
        "log.initial.Sun.LostEnergy": {"value": 5.562685e-309, "unit": u.Joule},
        "log.initial.Sun.LostAngMom": {
            "value": 5.562685e-309,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.Sun.Luminosity": {"value": 1.933219, "unit": u.LSUN},
        "log.initial.Sun.LXUVStellar": {"value": 0.001933, "unit": u.LSUN},
        "log.initial.Sun.Temperature": {"value": 4377.256537, "unit": u.K},
        "log.initial.Sun.LXUVFrac": {"value": 0.001000},
        "log.initial.Sun.RossbyNumber": {"value": 0.029996},
        "log.initial.Sun.DRotPerDtStellar": {"value": 9.510363e-10},
        "log.initial.Earth.Mass": {"value": 5.972186e24, "unit": u.kg},
        "log.initial.Earth.Radius": {"value": 6.378100e06, "unit": u.m},
        "log.initial.Earth.RadGyra": {"value": 0.500000},
        "log.initial.Earth.BodyType": {"value": 0.000000},
        "log.initial.Earth.Density": {"value": 5495.038549, "unit": u.kg / u.m**3},
        "log.initial.Earth.HZLimitDryRunaway": {"value": 1.888204e11, "unit": u.m},
        "log.initial.Earth.HZLimRecVenus": {"value": 1.640992e11, "unit": u.m},
        "log.initial.Earth.HZLimRunaway": {"value": 2.178515e11, "unit": u.m},
        "log.initial.Earth.HZLimMoistGreenhouse": {"value": 2.171294e11, "unit": u.m},
        "log.initial.Earth.HZLimMaxGreenhouse": {"value": 3.929566e11, "unit": u.m},
        "log.initial.Earth.HZLimEarlyMars": {"value": 4.286401e11, "unit": u.m},
        "log.initial.Earth.Instellation": {
            "value": 2644.188161,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.Earth.MeanMotion": {"value": 1.990987e-07, "unit": 1 / u.sec},
        "log.initial.Earth.OrbPeriod": {"value": 3.155815e07, "unit": u.sec},
        "log.initial.Earth.SemiMajorAxis": {"value": 1.495979e11, "unit": u.m},
        "log.initial.Earth.LXUVTot": {"value": -1.000000, "unit": u.kg / u.sec**3},
        "log.initial.Earth.SurfWaterMass": {"value": 6.950000e21, "unit": u.kg},
        "log.initial.Earth.EnvelopeMass": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.OxygenMass": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.RGLimit": {"value": 2.116533e11, "unit": u.m},
        "log.initial.Earth.XO": {"value": 0.333333},
        "log.initial.Earth.EtaO": {"value": 0.612158},
        "log.initial.Earth.PlanetRadius": {"value": 6.378100e06, "unit": u.m},
        "log.initial.Earth.OxygenMantleMass": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.RadXUV": {"value": -1.000000, "unit": u.m},
        "log.initial.Earth.RadSolid": {"value": -1.000000, "unit": u.m},
        "log.initial.Earth.PresXUV": {"value": 5.000000},
        "log.initial.Earth.ScaleHeight": {"value": -1.000000, "unit": u.m},
        "log.initial.Earth.ThermTemp": {"value": 195.374154, "unit": u.K},
        "log.initial.Earth.AtmGasConst": {"value": 4124.000000},
        "log.initial.Earth.PresSurf": {"value": -1.000000, "unit": u.Pa},
        "log.initial.Earth.DEnvMassDt": {"value": 0.000000, "unit": u.kg / u.sec},
        "log.initial.Earth.FXUV": {"value": 2.644188, "unit": u.W / u.m**2},
        "log.initial.Earth.AtmXAbsEffH2O": {"value": 0.300000},
        "log.initial.Earth.RocheRadius": {"value": 1.496558e09, "unit": u.m},
        "log.initial.Earth.BondiRadius": {"value": 1.014863e07, "unit": u.m},
        "log.initial.Earth.HEscapeRegime": {"value": 8.000000},
        "log.initial.Earth.RRCriticalFlux": {
            "value": 53.023382,
            "unit": u.W / u.m**2,
        },
        "log.initial.Earth.CrossoverMass": {"value": 6.588282e-26, "unit": u.kg},
        "log.initial.Earth.WaterEscapeRegime": {"value": 3.000000},
        "log.initial.Earth.FXUVCRITDRAG": {"value": 0.177574, "unit": u.W / u.m**2},
        "log.initial.Earth.HREFFLUX": {
            "value": 1.923284e18,
            "unit": 1 / u.m**2 / u.sec,
        },
        "log.initial.Earth.XO2": {"value": 0.000000},
        "log.initial.Earth.XH2O": {"value": 1.000000},
        "log.initial.Earth.HDiffFlux": {
            "value": 1.264874e17,
            "unit": 1 / u.m**2 / u.sec,
        },
        "log.initial.Earth.HRefODragMod": {"value": 0.169570},
        "log.initial.Earth.KTide": {"value": 0.993607},
        "log.initial.Earth.RGDuration": {"value": 0.00000e00, "unit": u.yr},
        "log.initial.Earth.WaterMassMOAtm": {"value": 5.000000, "unit": u.TO},
        "log.initial.Earth.PotTemp": {"value": 4000.000000, "unit": u.K},
        "log.initial.Earth.SurfTemp": {"value": 4000.000000, "unit": u.K},
        "log.initial.Earth.WaterMassSol": {"value": 0.000000, "unit": u.TO},
        "log.initial.Earth.SolidRadius": {"value": 0.533074, "unit": u.Rearth},
        "log.initial.Earth.OxygenMassMOAtm": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.OxygenMassSol": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.PressWaterAtm": {"value": 2.468645e06},
        "log.initial.Earth.PressOxygenAtm": {"value": 0.000000},
        "log.initial.Earth.HydrogenMassSpace": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.OxygenMassSpace": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.FracFe2O3Man": {"value": 0.000000},
        "log.initial.Earth.NetFluxAtmo": {"value": 1.340582e05},
        "log.initial.Earth.WaterFracMelt": {"value": 0.001849},
        "log.initial.Earth.RadioPower": {"value": 9.823813e13, "unit": u.W},
        "log.initial.Earth.TidalPower": {"value": 0.000000, "unit": u.W},
        "log.initial.Earth.HZInnerEdge": {"value": 2.116533e11, "unit": u.m},
        "log.initial.Earth.MeltFraction": {"value": 1.000000},
        "log.initial.Earth.CO2MassMOAtm": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.CO2MassSol": {"value": 0.000000, "unit": u.kg},
        "log.initial.Earth.PressCO2Atm": {"value": 0.000000},
        "log.initial.Earth.CO2FracMelt": {"value": 0.000000},
        "log.final.system.Age": {"value": 1.893456e13, "unit": u.sec, "rtol": 1e-4},
        "log.final.system.Time": {"value": 3.155760e12, "unit": u.sec, "rtol": 1e-4},
        "log.final.system.TotAngMom": {
            "value": 8.296053e43,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.system.TotEnergy": {
            "value": -9.104482e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.system.PotEnergy": {
            "value": -9.406038e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.system.KinEnergy": {
            "value": 2.813173e39,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.Sun.Mass": {"value": 1.988416e30, "unit": u.kg, "rtol": 1e-4},
        "log.final.Sun.Radius": {"value": 263.919878, "unit": u.Rearth, "rtol": 1e-4},
        "log.final.Sun.RadGyra": {"value": 0.449900, "rtol": 1e-4},
        "log.final.Sun.RotAngMom": {
            "value": 8.010254e43,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Sun.RotVel": {
            "value": 1.182343e05,
            "unit": u.m / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Sun.BodyType": {"value": 0.000000, "rtol": 1e-4},
        "log.final.Sun.RotRate": {
            "value": 7.023931e-05,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Sun.RotPer": {"value": 8.945398e04, "unit": u.sec, "rtol": 1e-4},
        "log.final.Sun.Density": {
            "value": 99.524124,
            "unit": u.kg / u.m**3,
            "rtol": 1e-4,
        },
        "log.final.Sun.HZLimitDryRunaway": {
            "value": 1.887931e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Sun.HZLimRecVenus": {
            "value": 1.640992e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Sun.HZLimRunaway": {"value": 2.178515e11, "unit": u.m, "rtol": 1e-4},
        "log.final.Sun.HZLimMoistGreenhouse": {
            "value": 2.171294e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Sun.HZLimMaxGreenhouse": {
            "value": 3.929566e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Sun.HZLimEarlyMars": {
            "value": 4.286401e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Sun.Instellation": {
            "value": -1.000000,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.Sun.CriticalSemiMajorAxis": {
            "value": -1.000000,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Sun.LXUVTot": {
            "value": 7.435159e23,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.Sun.LostEnergy": {
            "value": 2.023891e38,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.Sun.LostAngMom": {
            "value": 2.831382e42,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Sun.Luminosity": {"value": 1.933219, "unit": u.LSUN, "rtol": 1e-4},
        "log.final.Sun.LXUVStellar": {"value": 0.001933, "unit": u.LSUN, "rtol": 1e-4},
        "log.final.Sun.Temperature": {"value": 4377.256537, "unit": u.K, "rtol": 1e-4},
        "log.final.Sun.LXUVFrac": {"value": 0.001000, "rtol": 1e-4},
        "log.final.Sun.RossbyNumber": {"value": 0.031056, "rtol": 1e-4},
        "log.final.Sun.DRotPerDtStellar": {"value": 9.846526e-10, "rtol": 1e-4},
        "log.final.Earth.Mass": {"value": 5.972186e24, "unit": u.kg, "rtol": 1e-4},
        "log.final.Earth.Radius": {"value": 6.378100e06, "unit": u.m, "rtol": 1e-4},
        "log.final.Earth.RadGyra": {"value": 0.500000, "rtol": 1e-4},
        "log.final.Earth.BodyType": {"value": 0.000000, "rtol": 1e-4},
        "log.final.Earth.Density": {
            "value": 5495.038549,
            "unit": u.kg / u.m**3,
            "rtol": 1e-4,
        },
        "log.final.Earth.HZLimitDryRunaway": {
            "value": 1.888204e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.HZLimRecVenus": {
            "value": 1.640992e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.HZLimRunaway": {
            "value": 2.178515e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.HZLimMoistGreenhouse": {
            "value": 2.171294e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.HZLimMaxGreenhouse": {
            "value": 3.929566e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.HZLimEarlyMars": {
            "value": 4.286401e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.Instellation": {
            "value": 2644.188161,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.Earth.MeanMotion": {
            "value": 1.990987e-07,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Earth.OrbPeriod": {
            "value": 3.155815e07,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.Earth.SemiMajorAxis": {
            "value": 1.495979e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.LXUVTot": {
            "value": -1.000000,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.Earth.SurfWaterMass": {
            "value": 6.950000e21,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.EnvelopeMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4},
        "log.final.Earth.OxygenMass": {"value": 0.000000, "unit": u.kg, "rtol": 1e-4},
        "log.final.Earth.RGLimit": {"value": 2.116533e11, "unit": u.m, "rtol": 1e-4},
        "log.final.Earth.XO": {"value": 0.333333, "rtol": 1e-4},
        "log.final.Earth.EtaO": {"value": 0.612158, "rtol": 1e-4},
        "log.final.Earth.PlanetRadius": {
            "value": 6.378100e06,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.OxygenMantleMass": {
            "value": 0.000000,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.RadXUV": {"value": -1.000000, "unit": u.m, "rtol": 1e-4},
        "log.final.Earth.RadSolid": {"value": -1.000000, "unit": u.m, "rtol": 1e-4},
        "log.final.Earth.PresXUV": {"value": 5.000000, "rtol": 1e-4},
        "log.final.Earth.ScaleHeight": {"value": -1.000000, "unit": u.m, "rtol": 1e-4},
        "log.final.Earth.ThermTemp": {"value": 195.374154, "unit": u.K, "rtol": 1e-4},
        "log.final.Earth.AtmGasConst": {"value": 4124.000000, "rtol": 1e-4},
        "log.final.Earth.PresSurf": {"value": -1.000000, "unit": u.Pa, "rtol": 1e-4},
        "log.final.Earth.DEnvMassDt": {
            "value": 0.000000,
            "unit": u.kg / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Earth.FXUV": {
            "value": 2.644188,
            "unit": u.W / u.m**2,
            "rtol": 1e-4,
        },
        "log.final.Earth.AtmXAbsEffH2O": {"value": 0.300000, "rtol": 1e-4},
        "log.final.Earth.RocheRadius": {
            "value": 1.496558e09,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.BondiRadius": {
            "value": 1.014863e07,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.HEscapeRegime": {"value": 8.000000, "rtol": 1e-4},
        "log.final.Earth.RRCriticalFlux": {
            "value": 53.023382,
            "unit": u.W / u.m**2,
            "rtol": 1e-4,
        },
        "log.final.Earth.CrossoverMass": {
            "value": 6.588282e-26,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.WaterEscapeRegime": {"value": 3.000000, "rtol": 1e-4},
        "log.final.Earth.FXUVCRITDRAG": {
            "value": 0.177574,
            "unit": u.W / u.m**2,
            "rtol": 1e-4,
        },
        "log.final.Earth.HREFFLUX": {
            "value": 1.923284e18,
            "unit": 1 / u.m**2 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Earth.XO2": {"value": 0.000000, "rtol": 1e-4},
        "log.final.Earth.XH2O": {"value": 1.000000, "rtol": 1e-4},
        "log.final.Earth.HDiffFlux": {
            "value": 1.264874e17,
            "unit": 1 / u.m**2 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.Earth.HRefODragMod": {"value": 0.169570, "rtol": 1e-4},
        "log.final.Earth.KTide": {"value": 0.993607, "rtol": 1e-4},
        "log.final.Earth.RGDuration": {"value": 0.00000e00, "unit": u.yr, "rtol": 1e-4},
        "log.final.Earth.WaterMassMOAtm": {
            "value": 4.892788,
            "unit": u.TO,
            "rtol": 1e-4,
        },
        "log.final.Earth.PotTemp": {"value": 2173.620200, "unit": u.K, "rtol": 1e-4},
        "log.final.Earth.SurfTemp": {"value": 4000.000000, "unit": u.K, "rtol": 1e-4},
        "log.final.Earth.WaterMassSol": {"value": 0.101555, "unit": u.TO, "rtol": 1e-4},
        "log.final.Earth.SolidRadius": {
            "value": 0.919480,
            "unit": u.Rearth,
            "rtol": 1e-4,
        },
        "log.final.Earth.OxygenMassMOAtm": {
            "value": 1.836758e18,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.OxygenMassSol": {
            "value": 8.739505e17,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.PressWaterAtm": {"value": 3.567104e07, "rtol": 1e-4},
        "log.final.Earth.PressOxygenAtm": {"value": 0.000000, "rtol": 1e-4},
        "log.final.Earth.HydrogenMassSpace": {
            "value": 8.798414e17,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.OxygenMassSpace": {
            "value": 4.272139e18,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.FracFe2O3Man": {"value": 1.894182e-05, "rtol": 1e-4},
        "log.final.Earth.NetFluxAtmo": {"value": 815.884123, "rtol": 1e-4},
        "log.final.Earth.WaterFracMelt": {"value": 0.013344, "rtol": 1e-4},
        "log.final.Earth.RadioPower": {"value": 9.823253e13, "unit": u.W, "rtol": 1e-4},
        "log.final.Earth.TidalPower": {"value": 0.000000, "unit": u.W, "rtol": 1e-4},
        "log.final.Earth.HZInnerEdge": {
            "value": 2.116533e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.Earth.MeltFraction": {"value": 0.376273, "rtol": 1e-4},
        "log.final.Earth.CO2MassMOAtm": {
            "value": 1.753694e-296,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.CO2MassSol": {
            "value": 1.753694e-296,
            "unit": u.kg,
            "rtol": 1e-4,
        },
        "log.final.Earth.PressCO2Atm": {"value": 0.000000, "rtol": 1e-4},
        "log.final.Earth.CO2FracMelt": {"value": 0.000000, "rtol": 1e-4},
    }
)
class Test_MagmOc_Earth(Benchmark):
    pass
