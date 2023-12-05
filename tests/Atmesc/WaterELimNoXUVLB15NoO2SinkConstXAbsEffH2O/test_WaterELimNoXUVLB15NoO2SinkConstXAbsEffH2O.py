import astropy.units as u
import pytest
from benchmark import Benchmark, benchmark


@benchmark(
    {
        "log.initial.system.Age": {"value": 0.000000, "unit": u.sec},
        "log.initial.system.Time": {"value": 0.000000, "unit": u.sec},
        "log.initial.system.TotAngMom": {
            "value": 4.416946e33,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.system.TotEnergy": {"value": -2.237790e32, "unit": u.Joule},
        "log.initial.system.PotEnergy": {"value": -2.239397e32, "unit": u.Joule},
        "log.initial.system.KinEnergy": {"value": 1.606047e29, "unit": u.Joule},
        "log.initial.system.DeltaTime": {"value": 0.000000, "unit": u.sec},
        "log.initial.earth.Mass": {"value": 5.972186e24, "unit": u.kg},
        "log.initial.earth.Radius": {"value": 6.378100e06, "unit": u.m},
        "log.initial.earth.RadGyra": {"value": 0.500000},
        "log.initial.earth.BodyType": {"value": 0.000000},
        "log.initial.earth.Density": {"value": 5495.038549, "unit": u.kg / u.m**3},
        "log.initial.earth.HZLimitDryRunaway": {"value": -1.000000, "unit": u.m},
        "log.initial.earth.HZLimRecVenus": {"value": -1.000000},
        "log.initial.earth.HZLimRunaway": {"value": -1.000000},
        "log.initial.earth.HZLimMoistGreenhouse": {"value": -1.000000},
        "log.initial.earth.HZLimMaxGreenhouse": {"value": -1.000000},
        "log.initial.earth.HZLimEarlyMars": {"value": -1.000000},
        "log.initial.earth.Instellation": {
            "value": -1.000000,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.earth.MeanMotion": {"value": -1.000000, "unit": 1 / u.sec},
        "log.initial.earth.OrbPeriod": {"value": -1.000000, "unit": u.sec},
        "log.initial.earth.SemiMajorAxis": {"value": -1.000000, "unit": u.m},
        "log.initial.earth.LXUVTot": {"value": -1.000000, "unit": u.kg / u.sec**3},
        "log.initial.earth.SurfWaterMass": {"value": 3.000000, "unit": u.TO},
        "log.initial.earth.EnvelopeMass": {"value": 0.000000, "unit": u.kg},
        "log.initial.earth.OxygenMass": {"value": 0.000000, "unit": u.bar},
        "log.initial.earth.RGLimit": {"value": 0.000000, "unit": u.au},
        "log.initial.earth.XO": {"value": 0.333333},
        "log.initial.earth.EtaO": {"value": 0.000000},
        "log.initial.earth.PlanetRadius": {"value": 6.378100e06, "unit": u.m},
        "log.initial.earth.OxygenMantleMass": {"value": 0.000000, "unit": u.kg},
        "log.initial.earth.RadXUV": {"value": -1.000000, "unit": u.m},
        "log.initial.earth.RadSolid": {"value": -1.000000, "unit": u.m},
        "log.initial.earth.PresXUV": {"value": 5.000000},
        "log.initial.earth.ScaleHeight": {"value": -1.000000, "unit": u.m},
        "log.initial.earth.ThermTemp": {"value": 400.000000, "unit": u.K},
        "log.initial.earth.AtmGasConst": {"value": 4124.000000},
        "log.initial.earth.PresSurf": {"value": -1.000000, "unit": u.Pa},
        "log.initial.earth.DEnvMassDt": {"value": 0.000000, "unit": u.kg / u.sec},
        "log.initial.earth.FXUV": {"value": 0.000000, "unit": u.W / u.m**2},
        "log.initial.earth.AtmXAbsEffH2O": {"value": 0.100000},
        "log.initial.earth.RocheRadius": {"value": 1.037254e11, "unit": u.m},
        "log.initial.earth.BondiRadius": {"value": 1.249016e08, "unit": u.m},
        "log.initial.earth.HEscapeRegime": {"value": 8.000000},
        "log.initial.earth.RRCriticalFlux": {
            "value": 53.697959,
            "unit": u.W / u.m**2,
        },
        "log.initial.earth.CrossoverMass": {"value": 0.000000, "unit": u.kg},
        "log.initial.earth.WaterEscapeRegime": {"value": 8.000000},
        "log.initial.earth.FXUVCRITDRAG": {"value": 0.532721, "unit": u.W / u.m**2},
        "log.initial.earth.HREFFLUX": {"value": 0.000000, "unit": 1 / u.m**2 / u.sec},
        "log.initial.earth.XO2": {"value": 0.000000},
        "log.initial.earth.XH2O": {"value": 1.000000},
        "log.initial.earth.HDiffFlux": {
            "value": 1.264874e17,
            "unit": 1 / u.m**2 / u.sec,
        },
        "log.initial.earth.HRefODragMod": {"value": 0.007752},
        "log.initial.earth.KTide": {"value": 0.999908},
        "log.initial.earth.RGDuration": {"value": 0.00000e00, "unit": u.yr},
        "log.final.system.Age": {"value": 3.155760e15, "unit": u.sec},
        "log.final.system.Time": {"value": 3.155760e15, "unit": u.sec},
        "log.final.system.TotAngMom": {
            "value": 4.416946e33,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.final.system.TotEnergy": {"value": -2.237790e32, "unit": u.Joule},
        "log.final.system.PotEnergy": {"value": -2.239397e32, "unit": u.Joule},
        "log.final.system.KinEnergy": {"value": 1.606047e29, "unit": u.Joule},
        "log.final.system.DeltaTime": {"value": 3.155760e15, "unit": u.sec},
        "log.final.earth.Mass": {"value": 5.972186e24, "unit": u.kg},
        "log.final.earth.Radius": {"value": 6.378100e06, "unit": u.m},
        "log.final.earth.RadGyra": {"value": 0.500000},
        "log.final.earth.BodyType": {"value": 0.000000},
        "log.final.earth.Density": {"value": 5495.038549, "unit": u.kg / u.m**3},
        "log.final.earth.HZLimitDryRunaway": {"value": -1.000000, "unit": u.m},
        "log.final.earth.HZLimRecVenus": {"value": -1.000000},
        "log.final.earth.HZLimRunaway": {"value": -1.000000},
        "log.final.earth.HZLimMoistGreenhouse": {"value": -1.000000},
        "log.final.earth.HZLimMaxGreenhouse": {"value": -1.000000},
        "log.final.earth.HZLimEarlyMars": {"value": -1.000000},
        "log.final.earth.Instellation": {"value": -1.000000, "unit": u.kg / u.sec**3},
        "log.final.earth.MeanMotion": {"value": -1.000000, "unit": 1 / u.sec},
        "log.final.earth.OrbPeriod": {"value": -1.000000, "unit": u.sec},
        "log.final.earth.SemiMajorAxis": {"value": -1.000000, "unit": u.m},
        "log.final.earth.LXUVTot": {"value": -1.000000, "unit": u.kg / u.sec**3},
        "log.final.earth.SurfWaterMass": {"value": 3.000000, "unit": u.TO},
        "log.final.earth.EnvelopeMass": {"value": 0.000000, "unit": u.kg},
        "log.final.earth.OxygenMass": {"value": 0.000000, "unit": u.bar},
        "log.final.earth.RGLimit": {"value": 0.000000, "unit": u.au},
        "log.final.earth.XO": {"value": 0.333333},
        "log.final.earth.EtaO": {"value": 0.000000},
        "log.final.earth.PlanetRadius": {"value": 6.378100e06, "unit": u.m},
        "log.final.earth.OxygenMantleMass": {"value": 0.000000, "unit": u.kg},
        "log.final.earth.RadXUV": {"value": -1.000000, "unit": u.m},
        "log.final.earth.RadSolid": {"value": -1.000000, "unit": u.m},
        "log.final.earth.PresXUV": {"value": 5.000000},
        "log.final.earth.ScaleHeight": {"value": -1.000000, "unit": u.m},
        "log.final.earth.ThermTemp": {"value": 400.000000, "unit": u.K},
        "log.final.earth.AtmGasConst": {"value": 4124.000000},
        "log.final.earth.PresSurf": {"value": -1.000000, "unit": u.Pa},
        "log.final.earth.DEnvMassDt": {"value": 0.000000, "unit": u.kg / u.sec},
        "log.final.earth.FXUV": {"value": 0.000000, "unit": u.W / u.m**2},
        "log.final.earth.AtmXAbsEffH2O": {"value": 0.100000},
        "log.final.earth.RocheRadius": {"value": 1.037254e11, "unit": u.m},
        "log.final.earth.BondiRadius": {"value": 1.249016e08, "unit": u.m},
        "log.final.earth.HEscapeRegime": {"value": 8.000000},
        "log.final.earth.RRCriticalFlux": {"value": 53.697959, "unit": u.W / u.m**2},
        "log.final.earth.CrossoverMass": {"value": 0.000000, "unit": u.kg},
        "log.final.earth.WaterEscapeRegime": {"value": 8.000000},
        "log.final.earth.FXUVCRITDRAG": {"value": 0.532721, "unit": u.W / u.m**2},
        "log.final.earth.HREFFLUX": {"value": 0.000000, "unit": 1 / u.m**2 / u.sec},
        "log.final.earth.XO2": {"value": 0.000000},
        "log.final.earth.XH2O": {"value": 1.000000},
        "log.final.earth.HDiffFlux": {
            "value": 1.264874e17,
            "unit": 1 / u.m**2 / u.sec,
        },
        "log.final.earth.HRefODragMod": {"value": 0.007752},
        "log.final.earth.KTide": {"value": 0.999908},
        "log.final.earth.RGDuration": {"value": 1.00000e08, "unit": u.yr},
    }
)
class Test_WaterELimNoXUVLB15NoO2SinkConstXAbsEffH2O(Benchmark):
    pass
