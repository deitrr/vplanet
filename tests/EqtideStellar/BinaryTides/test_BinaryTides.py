import astropy.units as u
from benchmark import Benchmark, benchmark


@benchmark(
    {
        "log.initial.system.Age": {"value": 3.155760e13, "unit": u.sec},
        "log.initial.system.Time": {"value": 0.000000, "unit": u.sec},
        "log.initial.system.TotAngMom": {
            "value": 1.704021e45,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.system.TotEnergy": {"value": -1.981914e41, "unit": u.Joule},
        "log.initial.system.PotEnergy": {"value": -9.406038e40, "unit": u.Joule},
        "log.initial.system.KinEnergy": {"value": 1.081273e39, "unit": u.Joule},
        "log.initial.system.DeltaTime": {"value": 0.000000, "unit": u.sec},
        "log.initial.primary.Mass": {"value": 1.988416e30, "unit": u.kg},
        "log.initial.primary.Obliquity": {"value": 0.000000, "unit": u.rad},
        "log.initial.primary.PrecA": {"value": 0.000000, "unit": u.rad},
        "log.initial.primary.Xobl": {"value": 0.000000},
        "log.initial.primary.Yobl": {"value": 0.000000},
        "log.initial.primary.Zobl": {"value": 1.000000},
        "log.initial.primary.Radius": {"value": 263.919878, "unit": u.Rearth},
        "log.initial.primary.RadGyra": {"value": 0.449900},
        "log.initial.primary.RotAngMom": {
            "value": 4.966103e43,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.primary.RotKinEnergy": {"value": 1.081273e39, "unit": u.Joule},
        "log.initial.primary.RotVel": {"value": 7.330154e04, "unit": u.m / u.sec},
        "log.initial.primary.BodyType": {"value": 1.000000},
        "log.initial.primary.RotRate": {"value": 3.762386, "unit": 1 / u.day},
        "log.initial.primary.RotPer": {"value": 1.670000, "unit": u.day},
        "log.initial.primary.Density": {"value": 99.524124, "unit": u.kg / u.m**3},
        "log.initial.primary.SurfEnFluxTotal": {
            "value": 1.592014e07,
            "unit": u.W / u.m**2,
        },
        "log.initial.primary.TidalQ": {"value": 1.250000e05},
        "log.initial.primary.ImK2": {"value": -4.000000e-06},
        "log.initial.primary.K2": {"value": 0.500000},
        "log.initial.primary.K2Man": {"value": 0.010000},
        "log.initial.primary.Imk2Man": {"value": 0.000000},
        "log.initial.primary.TidalQMantle": {"value": 100.000000},
        "log.initial.primary.HEcc": {"value": 0.000000},
        "log.initial.primary.HZLimitDryRunaway": {"value": 1.887931e11, "unit": u.m},
        "log.initial.primary.HZLimRecVenus": {"value": 2.320714e11, "unit": u.m},
        "log.initial.primary.HZLimRunaway": {"value": 3.080886e11, "unit": u.m},
        "log.initial.primary.HZLimMoistGreenhouse": {"value": 3.070674e11, "unit": u.m},
        "log.initial.primary.HZLimMaxGreenhouse": {"value": 5.557245e11, "unit": u.m},
        "log.initial.primary.HZLimEarlyMars": {"value": 6.061886e11, "unit": u.m},
        "log.initial.primary.Instellation": {
            "value": -1.000000,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.primary.KEcc": {"value": 0.000000},
        "log.initial.primary.Eccentricity": {"value": -1.000000},
        "log.initial.primary.OrbEnergy": {"value": 0.000000, "unit": u.Joule},
        "log.initial.primary.MeanMotion": {"value": -1.000000, "unit": 1 / u.sec},
        "log.initial.primary.OrbPeriod": {"value": -1.000000, "unit": u.sec},
        "log.initial.primary.SemiMajorAxis": {"value": -1.000000, "unit": u.m},
        "log.initial.primary.CriticalSemiMajorAxis": {"value": -1.000000, "unit": u.m},
        "log.initial.primary.COPP": {"value": 0.000000},
        "log.initial.primary.OrbAngMom": {
            "value": 0.000000,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.primary.LongP": {"value": 0.000000, "unit": u.rad},
        "log.initial.primary.LXUVTot": {
            "value": 7.435159e23,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.primary.TotOrbEnergy": {"value": -1.223317e40, "unit": u.Joule},
        "log.initial.primary.OrbPotEnergy": {"value": -1.000000, "unit": u.Joule},
        "log.initial.primary.LostEnergy": {"value": 5.562685e-309, "unit": u.Joule},
        "log.initial.primary.LostAngMom": {
            "value": 5.562685e-309,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.primary.LockTime": {"value": -1.000000, "unit": u.sec},
        "log.initial.primary.BodyDsemiDtEqtide": {"value": -1.000000},
        "log.initial.primary.BodyDeccDt": {"value": -1.000000},
        "log.initial.primary.DOblDtEqtide": {"value": 0.000000, "unit": u.rad / u.sec},
        "log.initial.primary.DRotPerDtEqtide": {"value": 6.614815e-08},
        "log.initial.primary.DRotRateDtEqtide": {
            "value": -1.996352e-17,
            "unit": 1 / u.sec**2,
        },
        "log.initial.primary.EqRotRateDiscrete": {
            "value": 2.181662e-05,
            "unit": 1 / u.sec,
        },
        "log.initial.primary.EqRotPerDiscrete": {"value": 2.880000e05, "unit": u.sec},
        "log.initial.primary.EqRotRateCont": {"value": 2.697988e-05, "unit": 1 / u.sec},
        "log.initial.primary.EqRotPerCont": {"value": 2.328841e05, "unit": u.sec},
        "log.initial.primary.EqRotPer": {"value": 2.880000e05, "unit": u.sec},
        "log.initial.primary.EqTidePower": {"value": 0.000000, "unit": 1 / u.sec},
        "log.initial.primary.GammaRot": {"value": -1.000000, "unit": u.sec},
        "log.initial.primary.GammaOrb": {"value": -1.000000, "unit": u.sec},
        "log.initial.primary.OceanK2": {"value": 0.010000},
        "log.initial.primary.EnvTidalQ": {"value": -1.000000},
        "log.initial.primary.OceanTidalQ": {"value": -1.000000},
        "log.initial.primary.TideLock": {"value": 0.000000},
        "log.initial.primary.RotTimeEqtide": {"value": 2.181285e12, "unit": u.sec},
        "log.initial.primary.EnvK2": {"value": 0.010000},
        "log.initial.primary.OblTimeEqtide": {"value": -1.000000},
        "log.initial.primary.PowerEqtide": {"value": 5.668700e26, "unit": u.W},
        "log.initial.primary.SurfEnFluxEqtide": {
            "value": 1.592014e07,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.primary.Luminosity": {"value": 7.435159e26, "unit": u.W},
        "log.initial.primary.LXUVStellar": {"value": 7.435159e23, "unit": u.W},
        "log.initial.primary.Temperature": {"value": 4377.256537, "unit": u.K},
        "log.initial.primary.LXUVFrac": {"value": 0.001000},
        "log.initial.primary.RossbyNumber": {"value": 0.050093},
        "log.initial.primary.DRotPerDtStellar": {"value": -4.478836e-10},
        "log.initial.secondary.Mass": {"value": 1.988416e30, "unit": u.kg},
        "log.initial.secondary.Obliquity": {"value": 0.000000, "unit": u.rad},
        "log.initial.secondary.PrecA": {"value": 0.000000, "unit": u.rad},
        "log.initial.secondary.Xobl": {"value": 0.000000},
        "log.initial.secondary.Yobl": {"value": 0.000000},
        "log.initial.secondary.Zobl": {"value": 1.000000},
        "log.initial.secondary.Radius": {"value": 263.919878, "unit": u.Rearth},
        "log.initial.secondary.RadGyra": {"value": 0.449900},
        "log.initial.secondary.RotAngMom": {
            "value": 4.966103e43,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.secondary.RotKinEnergy": {"value": 1.081273e39, "unit": u.Joule},
        "log.initial.secondary.RotVel": {"value": 7.330154e04, "unit": u.m / u.sec},
        "log.initial.secondary.BodyType": {"value": 1.000000},
        "log.initial.secondary.RotRate": {"value": 3.762386, "unit": 1 / u.day},
        "log.initial.secondary.RotPer": {"value": 1.670000, "unit": u.day},
        "log.initial.secondary.Density": {"value": 99.524124, "unit": u.kg / u.m**3},
        "log.initial.secondary.SurfEnFluxTotal": {
            "value": 1.592014e07,
            "unit": u.W / u.m**2,
        },
        "log.initial.secondary.TidalQ": {"value": 1.250000e05},
        "log.initial.secondary.ImK2": {"value": -4.000000e-06},
        "log.initial.secondary.K2": {"value": 0.500000},
        "log.initial.secondary.K2Man": {"value": 0.010000},
        "log.initial.secondary.Imk2Man": {"value": 0.000000},
        "log.initial.secondary.TidalQMantle": {"value": 100.000000},
        "log.initial.secondary.HEcc": {"value": 0.000000},
        "log.initial.secondary.HZLimitDryRunaway": {"value": 1.979090e11, "unit": u.m},
        "log.initial.secondary.HZLimRecVenus": {"value": 2.320714e11, "unit": u.m},
        "log.initial.secondary.HZLimRunaway": {"value": 3.080886e11, "unit": u.m},
        "log.initial.secondary.HZLimMoistGreenhouse": {
            "value": 3.070674e11,
            "unit": u.m,
        },
        "log.initial.secondary.HZLimMaxGreenhouse": {"value": 5.557245e11, "unit": u.m},
        "log.initial.secondary.HZLimEarlyMars": {"value": 6.061886e11, "unit": u.m},
        "log.initial.secondary.Instellation": {
            "value": 5.331631e05,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.secondary.KEcc": {"value": 0.300000},
        "log.initial.secondary.Eccentricity": {"value": 0.300000},
        "log.initial.secondary.OrbEnergy": {"value": -1.223317e40, "unit": u.Joule},
        "log.initial.secondary.MeanMotion": {"value": 1.256637, "unit": 1 / u.day},
        "log.initial.secondary.OrbPeriod": {"value": 5.000000, "unit": u.day},
        "log.initial.secondary.SemiMajorAxis": {"value": 0.072098, "unit": u.au},
        "log.initial.secondary.CriticalSemiMajorAxis": {
            "value": 3.430862e10,
            "unit": u.m,
        },
        "log.initial.secondary.COPP": {"value": 0.000000},
        "log.initial.secondary.OrbAngMom": {
            "value": 1.604699e45,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.secondary.LongP": {"value": 0.000000, "unit": u.rad},
        "log.initial.secondary.LXUVTot": {
            "value": 7.435159e23,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.secondary.TotOrbEnergy": {"value": -1.223317e40, "unit": u.Joule},
        "log.initial.secondary.OrbPotEnergy": {"value": -2.446635e40, "unit": u.Joule},
        "log.initial.secondary.LostEnergy": {"value": 5.562685e-309, "unit": u.Joule},
        "log.initial.secondary.LostAngMom": {
            "value": 5.562685e-309,
            "unit": (u.kg * u.m**2) / u.sec,
        },
        "log.initial.secondary.TidalRadius": {"value": 1.683307e09, "unit": u.m},
        "log.initial.secondary.DsemiDtEqtide": {"value": 0.000749, "unit": u.m / u.sec},
        "log.initial.secondary.DeccDtEqtide": {
            "value": 2.302822e-14,
            "unit": 1 / u.sec,
        },
        "log.initial.secondary.DMeanMotionDtEqtide": {
            "value": -1.514244e-18,
            "unit": 1 / u.sec**2,
        },
        "log.initial.secondary.DOrbPerDtEqtide": {"value": 4.497629e-08},
        "log.initial.secondary.EccTimeEqtide": {"value": 1.302750e13, "unit": u.sec},
        "log.initial.secondary.SemiTimeEqtide": {"value": 1.440759e13, "unit": u.sec},
        "log.initial.secondary.DHEccDtEqtide": {"value": 0.000000, "unit": 1 / u.sec},
        "log.initial.secondary.DKEccDtEqtide": {
            "value": 2.302822e-14,
            "unit": 1 / u.sec,
        },
        "log.initial.secondary.DXoblDtEqtide": {"value": 0.000000, "unit": 1 / u.sec},
        "log.initial.secondary.DYoblDtEqtide": {"value": 0.000000, "unit": 1 / u.sec},
        "log.initial.secondary.DZoblDtEqtide": {"value": 0.000000, "unit": 1 / u.sec},
        "log.initial.secondary.LockTime": {"value": -1.000000, "unit": u.sec},
        "log.initial.secondary.BodyDsemiDtEqtide": {"value": -1.000000},
        "log.initial.secondary.BodyDeccDt": {"value": -1.000000},
        "log.initial.secondary.DOblDtEqtide": {
            "value": 0.000000,
            "unit": u.rad / u.sec,
        },
        "log.initial.secondary.DRotPerDtEqtide": {"value": 6.614815e-08},
        "log.initial.secondary.DRotRateDtEqtide": {
            "value": -1.996352e-17,
            "unit": 1 / u.sec**2,
        },
        "log.initial.secondary.EqRotRateDiscrete": {
            "value": 2.181662e-05,
            "unit": 1 / u.sec,
        },
        "log.initial.secondary.EqRotPerDiscrete": {"value": 2.880000e05, "unit": u.sec},
        "log.initial.secondary.EqRotRateCont": {
            "value": 2.697988e-05,
            "unit": 1 / u.sec,
        },
        "log.initial.secondary.EqRotPerCont": {"value": 2.328841e05, "unit": u.sec},
        "log.initial.secondary.EqRotPer": {"value": 2.880000e05, "unit": u.sec},
        "log.initial.secondary.EqTidePower": {"value": 1.897082e26, "unit": 1 / u.sec},
        "log.initial.secondary.GammaRot": {"value": -1.000000, "unit": u.sec},
        "log.initial.secondary.GammaOrb": {"value": -1.000000, "unit": u.sec},
        "log.initial.secondary.OceanK2": {"value": 0.010000},
        "log.initial.secondary.EnvTidalQ": {"value": -1.000000},
        "log.initial.secondary.OceanTidalQ": {"value": -1.000000},
        "log.initial.secondary.TideLock": {"value": 0.000000},
        "log.initial.secondary.RotTimeEqtide": {"value": 2.181285e12, "unit": u.sec},
        "log.initial.secondary.EnvK2": {"value": 0.010000},
        "log.initial.secondary.OblTimeEqtide": {"value": -1.000000},
        "log.initial.secondary.PowerEqtide": {"value": 5.668700e26, "unit": u.W},
        "log.initial.secondary.SurfEnFluxEqtide": {
            "value": 1.592014e07,
            "unit": u.kg / u.sec**3,
        },
        "log.initial.secondary.Luminosity": {"value": 7.435159e26, "unit": u.W},
        "log.initial.secondary.LXUVStellar": {"value": 7.435159e23, "unit": u.W},
        "log.initial.secondary.Temperature": {"value": 4377.256537, "unit": u.K},
        "log.initial.secondary.LXUVFrac": {"value": 0.001000},
        "log.initial.secondary.RossbyNumber": {"value": 0.050093},
        "log.initial.secondary.DRotPerDtStellar": {"value": -4.478836e-10},
        "log.final.system.Age": {"value": 6.311520e13, "unit": u.sec, "rtol": 1e-4},
        "log.final.system.Time": {"value": 3.155760e13, "unit": u.sec, "rtol": 1e-4},
        "log.final.system.TotAngMom": {
            "value": 1.533170e45,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.system.TotEnergy": {
            "value": -1.991596e41,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.system.PotEnergy": {
            "value": -1.186293e41,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.system.KinEnergy": {
            "value": 1.490449e38,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.primary.Mass": {"value": 1.988416e30, "unit": u.kg, "rtol": 1e-4},
        "log.final.primary.Obliquity": {"value": 0.000000, "unit": u.rad, "rtol": 1e-4},
        "log.final.primary.PrecA": {"value": 0.000000, "unit": u.rad, "rtol": 1e-4},
        "log.final.primary.Xobl": {"value": 0.000000, "rtol": 1e-4},
        "log.final.primary.Yobl": {"value": 0.000000, "rtol": 1e-4},
        "log.final.primary.Zobl": {"value": 1.000000, "rtol": 1e-4},
        "log.final.primary.Radius": {
            "value": 209.260399,
            "unit": u.Rearth,
            "rtol": 1e-4,
        },
        "log.final.primary.RadGyra": {"value": 0.451302, "rtol": 1e-4},
        "log.final.primary.RotAngMom": {
            "value": 1.466469e43,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.RotKinEnergy": {
            "value": 1.490449e38,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.primary.RotVel": {
            "value": 2.713019e04,
            "unit": u.m / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.BodyType": {"value": 1.000000, "rtol": 1e-4},
        "log.final.primary.RotRate": {
            "value": 1.756257,
            "unit": 1 / u.day,
            "rtol": 1e-4,
        },
        "log.final.primary.RotPer": {"value": 3.577599, "unit": u.day, "rtol": 1e-4},
        "log.final.primary.Density": {
            "value": 199.656530,
            "unit": u.kg / u.m**3,
            "rtol": 1e-4,
        },
        "log.final.primary.SurfEnFluxTotal": {
            "value": 2.137000e05,
            "unit": u.W / u.m**2,
            "rtol": 1e-4,
        },
        "log.final.primary.TidalQ": {"value": 1.250000e05, "rtol": 1e-4},
        "log.final.primary.ImK2": {"value": -4.000000e-06, "rtol": 1e-4},
        "log.final.primary.K2": {"value": 0.500000, "rtol": 1e-4},
        "log.final.primary.K2Man": {"value": 0.010000, "rtol": 1e-4},
        "log.final.primary.Imk2Man": {"value": 0.000000, "rtol": 1e-4},
        "log.final.primary.TidalQMantle": {"value": 100.000000, "rtol": 1e-4},
        "log.final.primary.HEcc": {"value": 0.000000, "rtol": 1e-4},
        "log.final.primary.HZLimitDryRunaway": {
            "value": 1.477883e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.HZLimRecVenus": {
            "value": 1.818297e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.HZLimRunaway": {
            "value": 2.413966e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.HZLimMoistGreenhouse": {
            "value": 2.405896e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.HZLimMaxGreenhouse": {
            "value": 4.359646e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.HZLimEarlyMars": {
            "value": 4.755510e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.Instellation": {
            "value": -1.000000,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.primary.KEcc": {"value": 0.000000, "rtol": 1e-4},
        "log.final.primary.Eccentricity": {"value": -1.000000, "rtol": 1e-4},
        "log.final.primary.OrbEnergy": {
            "value": 0.000000,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.primary.MeanMotion": {
            "value": -1.000000,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.OrbPeriod": {
            "value": -1.000000,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.SemiMajorAxis": {
            "value": -1.000000,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.CriticalSemiMajorAxis": {
            "value": -1.000000,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.primary.COPP": {"value": 0.000000, "rtol": 1e-4},
        "log.final.primary.OrbAngMom": {
            "value": 0.000000,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.LongP": {"value": 0.000000, "unit": u.rad, "rtol": 1e-4},
        "log.final.primary.LXUVTot": {
            "value": 4.556153e23,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.primary.TotOrbEnergy": {
            "value": -1.529178e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.primary.OrbPotEnergy": {
            "value": -1.000000,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.primary.LostEnergy": {
            "value": 2.654630e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.primary.LostAngMom": {
            "value": 4.120838e41,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.LockTime": {
            "value": 1.616555e12,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.BodyDsemiDtEqtide": {"value": -1.000000, "rtol": 1e-4},
        "log.final.primary.BodyDeccDt": {"value": -1.000000, "rtol": 1e-4},
        "log.final.primary.DOblDtEqtide": {
            "value": 0.000000,
            "unit": u.rad / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.DRotPerDtEqtide": {"value": -8.458931e-299, "rtol": 1e-4},
        "log.final.primary.DRotRateDtEqtide": {
            "value": 5.562685e-309,
            "unit": 1 / u.sec**2,
            "rtol": 1e-4,
        },
        "log.final.primary.EqRotRateDiscrete": {
            "value": 2.032705e-05,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.EqRotPerDiscrete": {
            "value": 3.091046e05,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.EqRotRateCont": {
            "value": 2.072656e-05,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.EqRotPerCont": {
            "value": 3.031466e05,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.EqRotPer": {
            "value": 3.091046e05,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.EqTidePower": {
            "value": 0.000000,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.GammaRot": {"value": -1.000000, "unit": u.sec, "rtol": 1e-4},
        "log.final.primary.GammaOrb": {"value": -1.000000, "unit": u.sec, "rtol": 1e-4},
        "log.final.primary.OceanK2": {"value": 0.010000, "rtol": 1e-4},
        "log.final.primary.EnvTidalQ": {"value": -1.000000, "rtol": 1e-4},
        "log.final.primary.OceanTidalQ": {"value": -1.000000, "rtol": 1e-4},
        "log.final.primary.TideLock": {"value": 1.000000, "rtol": 1e-4},
        "log.final.primary.RotTimeEqtide": {
            "value": 3.654180e303,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.primary.EnvK2": {"value": 0.010000, "rtol": 1e-4},
        "log.final.primary.OblTimeEqtide": {"value": -1.000000, "rtol": 1e-4},
        "log.final.primary.PowerEqtide": {
            "value": 4.783780e24,
            "unit": u.W,
            "rtol": 1e-4,
        },
        "log.final.primary.SurfEnFluxEqtide": {
            "value": 2.137000e05,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.primary.Luminosity": {
            "value": 4.556153e26,
            "unit": u.W,
            "rtol": 1e-4,
        },
        "log.final.primary.LXUVStellar": {
            "value": 4.556153e23,
            "unit": u.W,
            "rtol": 1e-4,
        },
        "log.final.primary.Temperature": {
            "value": 4349.796398,
            "unit": u.K,
            "rtol": 1e-4,
        },
        "log.final.primary.LXUVFrac": {"value": 0.001000, "rtol": 1e-4},
        "log.final.primary.RossbyNumber": {"value": 0.105795, "rtol": 1e-4},
        "log.final.primary.DRotPerDtStellar": {"value": -3.122954e-09, "rtol": 1e-4},
        "log.final.secondary.Mass": {"value": 1.988416e30, "unit": u.kg, "rtol": 1e-4},
        "log.final.secondary.Obliquity": {
            "value": 0.000000,
            "unit": u.rad,
            "rtol": 1e-4,
        },
        "log.final.secondary.PrecA": {"value": 0.000000, "unit": u.rad, "rtol": 1e-4},
        "log.final.secondary.Xobl": {"value": 0.000000, "rtol": 1e-4},
        "log.final.secondary.Yobl": {"value": 0.000000, "rtol": 1e-4},
        "log.final.secondary.Zobl": {"value": 1.000000, "rtol": 1e-4},
        "log.final.secondary.Radius": {
            "value": 209.260399,
            "unit": u.Rearth,
            "rtol": 1e-4,
        },
        "log.final.secondary.RadGyra": {"value": 0.451302, "rtol": 1e-4},
        "log.final.secondary.RotAngMom": {
            "value": 1.466469e43,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.RotKinEnergy": {
            "value": 1.490449e38,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.secondary.RotVel": {
            "value": 2.713019e04,
            "unit": u.m / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.BodyType": {"value": 1.000000, "rtol": 1e-4},
        "log.final.secondary.RotRate": {
            "value": 1.756257,
            "unit": 1 / u.day,
            "rtol": 1e-4,
        },
        "log.final.secondary.RotPer": {"value": 3.577599, "unit": u.day, "rtol": 1e-4},
        "log.final.secondary.Density": {
            "value": 199.656530,
            "unit": u.kg / u.m**3,
            "rtol": 1e-4,
        },
        "log.final.secondary.SurfEnFluxTotal": {
            "value": 2.137000e05,
            "unit": u.W / u.m**2,
            "rtol": 1e-4,
        },
        "log.final.secondary.TidalQ": {"value": 1.250000e05, "rtol": 1e-4},
        "log.final.secondary.ImK2": {"value": -4.000000e-06, "rtol": 1e-4},
        "log.final.secondary.K2": {"value": 0.500000, "rtol": 1e-4},
        "log.final.secondary.K2Man": {"value": 0.010000, "rtol": 1e-4},
        "log.final.secondary.Imk2Man": {"value": 0.000000, "rtol": 1e-4},
        "log.final.secondary.TidalQMantle": {"value": 100.000000, "rtol": 1e-4},
        "log.final.secondary.HEcc": {"value": 0.000000, "rtol": 1e-4},
        "log.final.secondary.HZLimitDryRunaway": {
            "value": 1.479415e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.HZLimRecVenus": {
            "value": 1.818297e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.HZLimRunaway": {
            "value": 2.413966e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.HZLimMoistGreenhouse": {
            "value": 2.405896e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.HZLimMaxGreenhouse": {
            "value": 4.359646e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.HZLimEarlyMars": {
            "value": 4.755510e11,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.Instellation": {
            "value": 4.875019e05,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.secondary.KEcc": {"value": 0.045484, "rtol": 1e-4},
        "log.final.secondary.Eccentricity": {"value": 0.045484, "rtol": 1e-4},
        "log.final.secondary.OrbEnergy": {
            "value": -1.529178e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.secondary.MeanMotion": {
            "value": 1.756257,
            "unit": 1 / u.day,
            "rtol": 1e-4,
        },
        "log.final.secondary.OrbPeriod": {
            "value": 3.577599,
            "unit": u.day,
            "rtol": 1e-4,
        },
        "log.final.secondary.SemiMajorAxis": {
            "value": 0.057677,
            "unit": u.au,
            "rtol": 1e-4,
        },
        "log.final.secondary.CriticalSemiMajorAxis": {
            "value": 2.174492e10,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.COPP": {"value": 0.000000, "rtol": 1e-4},
        "log.final.secondary.OrbAngMom": {
            "value": 1.503017e45,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.LongP": {"value": 0.000000, "unit": u.rad, "rtol": 1e-4},
        "log.final.secondary.LXUVTot": {
            "value": 4.556153e23,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.secondary.TotOrbEnergy": {
            "value": -1.529178e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.secondary.OrbPotEnergy": {
            "value": -3.058356e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.secondary.LostEnergy": {
            "value": 2.654630e40,
            "unit": u.Joule,
            "rtol": 1e-4,
        },
        "log.final.secondary.LostAngMom": {
            "value": 4.120838e41,
            "unit": (u.kg * u.m**2) / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.TidalRadius": {
            "value": 1.334684e09,
            "unit": u.m,
            "rtol": 1e-4,
        },
        "log.final.secondary.DsemiDtEqtide": {
            "value": -5.398513e-06,
            "unit": u.m / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DeccDtEqtide": {
            "value": -6.877840e-15,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DMeanMotionDtEqtide": {
            "value": 1.907695e-20,
            "unit": 1 / u.sec**2,
            "rtol": 1e-4,
        },
        "log.final.secondary.DOrbPerDtEqtide": {"value": -2.900948e-10, "rtol": 1e-4},
        "log.final.secondary.EccTimeEqtide": {
            "value": 0.000000,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.SemiTimeEqtide": {
            "value": 0.000000,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DHEccDtEqtide": {
            "value": -0.000000,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DKEccDtEqtide": {
            "value": -6.877840e-15,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DXoblDtEqtide": {
            "value": 0.000000,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DYoblDtEqtide": {
            "value": 0.000000,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DZoblDtEqtide": {
            "value": 0.000000,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.LockTime": {
            "value": 1.616555e12,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.BodyDsemiDtEqtide": {"value": -1.000000, "rtol": 1e-4},
        "log.final.secondary.BodyDeccDt": {"value": -1.000000, "rtol": 1e-4},
        "log.final.secondary.DOblDtEqtide": {
            "value": 0.000000,
            "unit": u.rad / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.DRotPerDtEqtide": {"value": -8.458931e-299, "rtol": 1e-4},
        "log.final.secondary.DRotRateDtEqtide": {
            "value": 5.562685e-309,
            "unit": 1 / u.sec**2,
            "rtol": 1e-4,
        },
        "log.final.secondary.EqRotRateDiscrete": {
            "value": 2.032705e-05,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.EqRotPerDiscrete": {
            "value": 3.091046e05,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.EqRotRateCont": {
            "value": 2.072656e-05,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.EqRotPerCont": {
            "value": 3.031466e05,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.EqRotPer": {
            "value": 3.091046e05,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.EqTidePower": {
            "value": 3.005173e24,
            "unit": 1 / u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.GammaRot": {
            "value": -1.000000,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.GammaOrb": {
            "value": -1.000000,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.OceanK2": {"value": 0.010000, "rtol": 1e-4},
        "log.final.secondary.EnvTidalQ": {"value": -1.000000, "rtol": 1e-4},
        "log.final.secondary.OceanTidalQ": {"value": -1.000000, "rtol": 1e-4},
        "log.final.secondary.TideLock": {"value": 1.000000, "rtol": 1e-4},
        "log.final.secondary.RotTimeEqtide": {
            "value": 3.654180e303,
            "unit": u.sec,
            "rtol": 1e-4,
        },
        "log.final.secondary.EnvK2": {"value": 0.010000, "rtol": 1e-4},
        "log.final.secondary.OblTimeEqtide": {"value": -1.000000, "rtol": 1e-4},
        "log.final.secondary.PowerEqtide": {
            "value": 4.783780e24,
            "unit": u.W,
            "rtol": 1e-4,
        },
        "log.final.secondary.SurfEnFluxEqtide": {
            "value": 2.137000e05,
            "unit": u.kg / u.sec**3,
            "rtol": 1e-4,
        },
        "log.final.secondary.Luminosity": {
            "value": 4.556153e26,
            "unit": u.W,
            "rtol": 1e-4,
        },
        "log.final.secondary.LXUVStellar": {
            "value": 4.556153e23,
            "unit": u.W,
            "rtol": 1e-4,
        },
        "log.final.secondary.Temperature": {
            "value": 4349.796398,
            "unit": u.K,
            "rtol": 1e-4,
        },
        "log.final.secondary.LXUVFrac": {"value": 0.001000, "rtol": 1e-4},
        "log.final.secondary.RossbyNumber": {"value": 0.105795, "rtol": 1e-4},
        "log.final.secondary.DRotPerDtStellar": {"value": -3.122954e-09, "rtol": 1e-4},
    }
)
class Test_BinaryTides(Benchmark):
    pass
