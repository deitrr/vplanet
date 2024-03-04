from benchmark import Benchmark, benchmark 
import astropy.units as u 
 
@benchmark( 
   { 
       "log.initial.system.Age": {"value": 0.0000000000000000, "unit": u.sec}, 
       "log.initial.system.Time": {"value": 0.0000000000000000, "unit": u.sec}, 
       "log.initial.system.TotAngMom": {"value": 8.9534036048945070e+41, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.system.TotEnergy": {"value": -3.1347047276701209e+40, "unit": u.Joule}, 
       "log.initial.system.PotEnergy": {"value": -3.1348025171406187e+40, "unit": u.Joule}, 
       "log.initial.system.KinEnergy": {"value": 9.9742042956183368e+35, "unit": u.Joule}, 
       "log.initial.system.DeltaTime": {"value": 0.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514.Mass": {"value": 1.0140921600000001e+30, "unit": u.kg}, 
       "log.initial.gl514.Obliquity": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514.PrecA": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514.Xobl": {"value": 0.0000000000000000}, 
       "log.initial.gl514.Yobl": {"value": 0.0000000000000000}, 
       "log.initial.gl514.Zobl": {"value": 1.0000000000000000}, 
       "log.initial.gl514.Radius": {"value": 1.3137125118741000e+09, "unit": u.m}, 
       "log.initial.gl514.RadGyra": {"value": 0.4484976455000000}, 
       "log.initial.gl514.RotAngMom": {"value": 8.3801799358956994e+41, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.gl514.RotKinEnergy": {"value": 9.9742042956183368e+35, "unit": u.Joule}, 
       "log.initial.gl514.RotVel": {"value": 3127.1970481244061375, "unit": u.m / u.sec}, 
       "log.initial.gl514.BodyType": {"value": 0.0000000000000000}, 
       "log.initial.gl514.RotRate": {"value": 2.3804272394903566e-06, "unit": 1 / u.sec}, 
       "log.initial.gl514.RotPer": {"value": 30.5499999999999972, "unit": u.day}, 
       "log.initial.gl514.Density": {"value": 106.7794814169220103, "unit": u.kg / u.m ** 3}, 
       "log.initial.gl514.SurfEnFluxTotal": {"value": 5.6927005838186184e-10, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514.TidalQ": {"value": 1.0000000000000000e+06}, 
       "log.initial.gl514.ImK2": {"value": -4.9999999999999998e-07}, 
       "log.initial.gl514.K2": {"value": 0.5000000000000000}, 
       "log.initial.gl514.K2Man": {"value": 0.0100000000000000}, 
       "log.initial.gl514.Imk2Man": {"value": 0.0000000000000000}, 
       "log.initial.gl514.TidalQMantle": {"value": 100.0000000000000000}, 
       "log.initial.gl514.HEcc": {"value": 0.0000000000000000}, 
       "log.initial.gl514.HZLimitDryRunaway": {"value": 1.1215301622470532e+11, "unit": u.m}, 
       "log.initial.gl514.HZLimRecVenus": {"value": 9.9104353063518631e+10, "unit": u.m}, 
       "log.initial.gl514.HZLimRunaway": {"value": 1.3136272247075568e+11, "unit": u.m}, 
       "log.initial.gl514.HZLimMoistGreenhouse": {"value": 1.3113105818748537e+11, "unit": u.m}, 
       "log.initial.gl514.HZLimMaxGreenhouse": {"value": 2.4324580025875378e+11, "unit": u.m}, 
       "log.initial.gl514.HZLimEarlyMars": {"value": 2.6530630489685046e+11, "unit": u.m}, 
       "log.initial.gl514.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514.KEcc": {"value": 0.0000000000000000}, 
       "log.initial.gl514.Eccentricity": {"value": -1.0000000000000000}, 
       "log.initial.gl514.OrbEnergy": {"value": 0.0000000000000000, "unit": u.Joule}, 
       "log.initial.gl514.MeanMotion": {"value": -1.0000000000000000, "unit": 1 / u.sec}, 
       "log.initial.gl514.OrbPeriod": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514.SemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.gl514.CriticalSemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.gl514.COPP": {"value": 0.0000000000000000}, 
       "log.initial.gl514.OrbAngMom": {"value": 0.0000000000000000, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.gl514.LongP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514.LXUVTot": {"value": 2.6238552211501671e+23, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514.TotOrbEnergy": {"value": -1.6646340895690071e+34, "unit": u.Joule}, 
       "log.initial.gl514.OrbPotEnergy": {"value": -1.0000000000000000, "unit": u.Joule}, 
       "log.initial.gl514.LostEnergy": {"value": 5.5626846462680035e-309, "unit": u.Joule}, 
       "log.initial.gl514.LostAngMom": {"value": 5.5626846462680035e-309, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.gl514.LockTime": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514.BodyDsemiDtEqtide": {"value": -1.0000000000000000}, 
       "log.initial.gl514.BodyDeccDt": {"value": -1.0000000000000000}, 
       "log.initial.gl514.DOblDtEqtide": {"value": 0.0000000000000000, "unit": u.rad / u.sec}, 
       "log.initial.gl514.DRotPerDtEqtide": {"value": 2.3670446824643342e-20}, 
       "log.initial.gl514.DRotRateDtEqtide": {"value": -2.1346978387088952e-32, "unit": 1 / u.sec ** 2}, 
       "log.initial.gl514.EqRotRateDiscrete": {"value": 7.7800437205935984e-07, "unit": 1 / u.sec}, 
       "log.initial.gl514.EqRotPerDiscrete": {"value": 8.0760282754557505e+06, "unit": u.sec}, 
       "log.initial.gl514.EqRotRateCont": {"value": 1.5164601885390357e-06, "unit": 1 / u.sec}, 
       "log.initial.gl514.EqRotPerCont": {"value": 4.1433236128887986e+06, "unit": u.sec}, 
       "log.initial.gl514.EqRotPer": {"value": 4.1433236128887986e+06, "unit": u.sec}, 
       "log.initial.gl514.EqTidePower": {"value": 0.0000000000000000, "unit": 1 / u.sec}, 
       "log.initial.gl514.GammaRot": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514.GammaOrb": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514.OceanK2": {"value": 0.0100000000000000}, 
       "log.initial.gl514.EnvTidalQ": {"value": -1.0000000000000000}, 
       "log.initial.gl514.OceanTidalQ": {"value": -1.0000000000000000}, 
       "log.initial.gl514.TideLock": {"value": 0.0000000000000000}, 
       "log.initial.gl514.RotTimeEqtide": {"value": 1.1151120295929485e+26, "unit": u.sec}, 
       "log.initial.gl514.EnvK2": {"value": 0.0100000000000000}, 
       "log.initial.gl514.OblTimeEqtide": {"value": -1.0000000000000000}, 
       "log.initial.gl514.PowerEqtide": {"value": 1.2346074076706444e+10, "unit": u.W}, 
       "log.initial.gl514.SurfEnFluxEqtide": {"value": 5.6927005838186184e-10, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514.Luminosity": {"value": 2.6238552211501671e+26, "unit": u.W}, 
       "log.initial.gl514.LXUVStellar": {"value": 2.6238552211501671e+23, "unit": u.W}, 
       "log.initial.gl514.Temperature": {"value": 3818.8020013755422042, "unit": u.K}, 
       "log.initial.gl514.LXUVFrac": {"value": 0.0010000000000000}, 
       "log.initial.gl514.RossbyNumber": {"value": 0.6873965284599931}, 
       "log.initial.gl514.DRotPerDtStellar": {"value": 2.3504479341501690e-10}, 
       "log.initial.gl514b.Mass": {"value": 3.1053088098643201e+25, "unit": u.kg}, 
       "log.initial.gl514b.Obliquity": {"value": 0.4101523742069445, "unit": u.rad}, 
       "log.initial.gl514b.PrecA": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.Xobl": {"value": 0.3987490689144955}, 
       "log.initial.gl514b.Yobl": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.Zobl": {"value": 0.9170600743897985}, 
       "log.initial.gl514b.Radius": {"value": 1.3394010000000000e+07, "unit": u.m}, 
       "log.initial.gl514b.RadGyra": {"value": 0.5000000000000000}, 
       "log.initial.gl514b.RotAngMom": {"value": 1.0128197643484414e+35, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.gl514b.RotKinEnergy": {"value": 3.6827165869069548e+30, "unit": u.Joule}, 
       "log.initial.gl514b.RotVel": {"value": 974.0398939376904082, "unit": u.m / u.sec}, 
       "log.initial.gl514b.BodyType": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.RotRate": {"value": 7.2722052166430395e-05, "unit": 1 / u.sec}, 
       "log.initial.gl514b.RotPer": {"value": 8.6400000000000000e+04, "unit": u.sec}, 
       "log.initial.gl514b.Density": {"value": 3085.2071528551509800, "unit": u.kg / u.m ** 3}, 
       "log.initial.gl514b.SurfEnFluxTotal": {"value": 1.4096431035620620, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.TidalQ": {"value": 12.0000000000000000}, 
       "log.initial.gl514b.ImK2": {"value": -0.0250000000000000}, 
       "log.initial.gl514b.K2": {"value": 0.3000000000000000}, 
       "log.initial.gl514b.K2Man": {"value": 0.0100000000000000}, 
       "log.initial.gl514b.Imk2Man": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.TidalQMantle": {"value": 100.0000000000000000}, 
       "log.initial.gl514b.HEcc": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.HZLimitDryRunaway": {"value": 1.0509802225988144e+11, "unit": u.m}, 
       "log.initial.gl514b.HZLimRecVenus": {"value": 9.9104353063518631e+10, "unit": u.m}, 
       "log.initial.gl514b.HZLimRunaway": {"value": 1.3136272247075568e+11, "unit": u.m}, 
       "log.initial.gl514b.HZLimMoistGreenhouse": {"value": 1.3113105818748537e+11, "unit": u.m}, 
       "log.initial.gl514b.HZLimMaxGreenhouse": {"value": 2.4324580025875378e+11, "unit": u.m}, 
       "log.initial.gl514b.HZLimEarlyMars": {"value": 2.6530630489685046e+11, "unit": u.m}, 
       "log.initial.gl514b.Instellation": {"value": 5866.6331509136580280, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.KEcc": {"value": 0.4500000000000000}, 
       "log.initial.gl514b.Eccentricity": {"value": 0.4500000000000000}, 
       "log.initial.gl514b.OrbEnergy": {"value": -1.6646340895690071e+34, "unit": u.Joule}, 
       "log.initial.gl514b.MeanMotion": {"value": 5.1866958137290656e-07, "unit": 1 / u.sec}, 
       "log.initial.gl514b.OrbPeriod": {"value": 1.2114042413183626e+07, "unit": u.sec}, 
       "log.initial.gl514b.SemiMajorAxis": {"value": 6.3130301435400002e+10, "unit": u.m}, 
       "log.initial.gl514b.CriticalSemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.gl514b.COPP": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.OrbAngMom": {"value": 5.7322265617904285e+40, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.gl514b.LongP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.TotOrbEnergy": {"value": -1.6646340895690071e+34, "unit": u.Joule}, 
       "log.initial.gl514b.OrbPotEnergy": {"value": -3.3292681791380142e+34, "unit": u.Joule}, 
       "log.initial.gl514b.LostEnergy": {"value": 5.5626846462680035e-309, "unit": u.Joule}, 
       "log.initial.gl514b.TidalRadius": {"value": 1.3394010000000000e+07, "unit": u.m}, 
       "log.initial.gl514b.DsemiDtEqtide": {"value": 1.2351242905116021e-10, "unit": u.m / u.sec}, 
       "log.initial.gl514b.DeccDtEqtide": {"value": 5.8376560316556904e-22, "unit": 1 / u.sec}, 
       "log.initial.gl514b.DMeanMotionDtEqtide": {"value": -1.5221408360263325e-27, "unit": 1 / u.sec ** 2}, 
       "log.initial.gl514b.DOrbPerDtEqtide": {"value": 3.5551108660842267e-14}, 
       "log.initial.gl514b.EccTimeEqtide": {"value": 7.7085733993198282e+20, "unit": u.sec}, 
       "log.initial.gl514b.SemiTimeEqtide": {"value": 5.1112509016603289e+20, "unit": u.sec}, 
       "log.initial.gl514b.DHEccDtEqtide": {"value": 0.0000000000000000, "unit": 1 / u.sec}, 
       "log.initial.gl514b.DKEccDtEqtide": {"value": 5.8376560316556904e-22, "unit": 1 / u.sec}, 
       "log.initial.gl514b.DXoblDtEqtide": {"value": 3.1641063776408136e-17, "unit": 1 / u.sec}, 
       "log.initial.gl514b.DYoblDtEqtide": {"value": 0.0000000000000000, "unit": 1 / u.sec}, 
       "log.initial.gl514b.DZoblDtEqtide": {"value": -1.3757926086469334e-17, "unit": 1 / u.sec}, 
       "log.initial.gl514b.LockTime": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514b.BodyDsemiDtEqtide": {"value": -1.0000000000000000}, 
       "log.initial.gl514b.BodyDeccDt": {"value": -1.0000000000000000}, 
       "log.initial.gl514b.DOblDtEqtide": {"value": 3.4502716517739303e-17, "unit": u.rad / u.sec}, 
       "log.initial.gl514b.DRotPerDtEqtide": {"value": 3.7660240995856781e-11}, 
       "log.initial.gl514b.DRotRateDtEqtide": {"value": -3.1698264008113854e-20, "unit": 1 / u.sec ** 2}, 
       "log.initial.gl514b.EqRotRateDiscrete": {"value": 7.7800437205935984e-07, "unit": 1 / u.sec}, 
       "log.initial.gl514b.EqRotPerDiscrete": {"value": 8.0760282754557505e+06, "unit": u.sec}, 
       "log.initial.gl514b.EqRotRateCont": {"value": 1.5164601885390357e-06, "unit": 1 / u.sec}, 
       "log.initial.gl514b.EqRotPerCont": {"value": 4.1433236128887986e+06, "unit": u.sec}, 
       "log.initial.gl514b.EqRotPer": {"value": 4.1433236128887986e+06, "unit": u.sec}, 
       "log.initial.gl514b.EqTidePower": {"value": 3.4384723793793242e+13, "unit": 1 / u.sec}, 
       "log.initial.gl514b.GammaRot": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514b.GammaOrb": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.gl514b.OceanK2": {"value": 0.0100000000000000}, 
       "log.initial.gl514b.EnvTidalQ": {"value": -1.0000000000000000}, 
       "log.initial.gl514b.OceanTidalQ": {"value": -1.0000000000000000}, 
       "log.initial.gl514b.TideLock": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.RotTimeEqtide": {"value": 2.2941966836990070e+15, "unit": u.sec}, 
       "log.initial.gl514b.EnvK2": {"value": 0.0100000000000000}, 
       "log.initial.gl514b.OblTimeEqtide": {"value": -1.0000000000000000}, 
       "log.initial.gl514b.PowerEqtide": {"value": 3.1779003342801765e+15, "unit": u.W}, 
       "log.initial.gl514b.SurfEnFluxEqtide": {"value": 1.4096431035620620, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.TGlobal": {"value": 676.2679353212223532, "unit": u.sec}, 
       "log.initial.gl514b.AlbedoGlobal": {"value": 0.2996778493802214}, 
       "log.initial.gl514b.FluxInGlobal": {"value": 1046.1316028565233864, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.FluxOutGlobal": {"value": 1046.1299848213550376, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.TotIceMass": {"value": 0.0000000000000000, "unit": u.kg}, 
       "log.initial.gl514b.TotIceFlow": {"value": 0.0000000000000000, "unit": u.kg}, 
       "log.initial.gl514b.TotIceBalance": {"value": 0.0000000000000000, "unit": u.kg}, 
       "log.initial.gl514b.SkipSeas": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.AreaIceCov": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.Latitude": {"value": -1.4552620265106593, "unit": u.rad}, 
       "log.initial.gl514b.TempLat": {"value": 557.2247797713517912, "unit": u.sec}, 
       "log.initial.gl514b.AlbedoLat": {"value": 0.3675649878041329}, 
       "log.initial.gl514b.AnnInsol": {"value": 755.3407450595735781, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.FluxMerid": {"value": -2.0352457824211096e+16, "unit": u.Joule}, 
       "log.initial.gl514b.FluxIn": {"value": 486.1623336259239636, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.FluxOut": {"value": 797.3297897221252697, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.DivFlux": {"value": -311.1690738881308107, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.IceMass": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceHeight": {"value": 0.0000000000000000, "unit": u.m}, 
       "log.initial.gl514b.DIceMassDt": {"value": 0.0000000000000000, "unit": u.m}, 
       "log.initial.gl514b.IceFlow": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.EnergyResL": {"value": -4.7634785005357116e-11, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.EnergyResW": {"value": 4.0836312109604478e-10, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.BedrockH": {"value": 0.0000000000000000, "unit": u.m}, 
       "log.initial.gl514b.TempLandLat": {"value": 545.8662422200939091, "unit": u.sec}, 
       "log.initial.gl514b.TempWaterLat": {"value": 563.0761476007876354, "unit": u.sec}, 
       "log.initial.gl514b.AlbedoLandLat": {"value": 0.4335649878041328}, 
       "log.initial.gl514b.AlbedoWaterLat": {"value": 0.3335649878041328}, 
       "log.initial.gl514b.TempMinLat": {"value": 520.6583069986143073, "unit": u.sec}, 
       "log.initial.gl514b.TempMaxLat": {"value": 598.1727520370602633, "unit": u.sec}, 
       "log.initial.gl514b.Snowball": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.PlanckBAvg": {"value": 2.0899999999999990}, 
       "log.initial.gl514b.IceAccum": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceAblate": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.TempMaxLand": {"value": 657.9976852023412448, "unit": u.sec}, 
       "log.initial.gl514b.TempMaxWater": {"value": 568.6086486017817379, "unit": u.sec}, 
       "log.initial.gl514b.PeakInsol": {"value": 4243.8792145061070187, "unit": u.kg / u.sec ** 3}, 
       "log.initial.gl514b.IceCapNorthLand": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceCapNorthSea": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceCapSouthLand": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceCapSouthSea": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceBeltLand": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceBeltSea": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.SnowballLand": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.SnowballSea": {"value": 0.0000000000000000}, 
       "log.initial.gl514b.IceFree": {"value": 1.0000000000000000}, 
       "log.initial.gl514b.IceCapNorthLatLand": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.IceCapNorthLatSea": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.IceCapSouthLatLand": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.IceCapSouthLatSea": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.IceBeltNorthLatLand": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.IceBeltNorthLatSea": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.IceBeltSouthLatLand": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.initial.gl514b.IceBeltSouthLatSea": {"value": 100.0000000000000000, "unit": u.rad}, 
       "log.final.system.Age": {"value": 3.1557600000000000e+07, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.system.Time": {"value": 3.1557600000000000e+07, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.system.TotAngMom": {"value": 8.9534036048945039e+41, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.system.TotEnergy": {"value": -3.1347047276701209e+40, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.system.PotEnergy": {"value": -3.1348025171406187e+40, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.system.KinEnergy": {"value": 9.9742042395602989e+35, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514.Mass": {"value": 1.0140921600000001e+30, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.gl514.Obliquity": {"value": 0.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514.PrecA": {"value": 0.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514.Xobl": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.Yobl": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.Zobl": {"value": 1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.Radius": {"value": 1.3137125118741000e+09, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.RadGyra": {"value": 0.4484976455000000, "rtol": 1e-4}, 
       "log.final.gl514.RotAngMom": {"value": 8.3801799123461292e+41, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.RotKinEnergy": {"value": 9.9742042395602989e+35, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514.RotVel": {"value": 3127.1970393365104428, "unit": u.m / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.BodyType": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.RotRate": {"value": 2.3804272328009969e-06, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.RotPer": {"value": 30.5500000858501117, "unit": u.day, "rtol": 1e-4}, 
       "log.final.gl514.Density": {"value": 106.7794814169220103, "unit": u.kg / u.m ** 3, "rtol": 1e-4}, 
       "log.final.gl514.SurfEnFluxTotal": {"value": 5.6927005606372488e-10, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514.TidalQ": {"value": 1.0000000000000000e+06, "rtol": 1e-4}, 
       "log.final.gl514.ImK2": {"value": -4.9999999999999998e-07, "rtol": 1e-4}, 
       "log.final.gl514.K2": {"value": 0.5000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.K2Man": {"value": 0.0100000000000000, "rtol": 1e-4}, 
       "log.final.gl514.Imk2Man": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.TidalQMantle": {"value": 100.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.HEcc": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.HZLimitDryRunaway": {"value": 1.1215301622470532e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.HZLimRecVenus": {"value": 9.9104353063518631e+10, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.HZLimRunaway": {"value": 1.3136272247075568e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.HZLimMoistGreenhouse": {"value": 1.3113105818748537e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.HZLimMaxGreenhouse": {"value": 2.4324580025875378e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.HZLimEarlyMars": {"value": 2.6530630489685046e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514.KEcc": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.Eccentricity": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.OrbEnergy": {"value": 0.0000000000000000, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514.MeanMotion": {"value": -1.0000000000000000, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.OrbPeriod": {"value": -1.0000000000000000, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.SemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.CriticalSemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514.COPP": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.OrbAngMom": {"value": 0.0000000000000000, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.LongP": {"value": 0.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514.LXUVTot": {"value": 2.6238552211501671e+23, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514.TotOrbEnergy": {"value": -1.6646340895689042e+34, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514.OrbPotEnergy": {"value": -1.0000000000000000, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514.LostEnergy": {"value": 5.6058038409805720e+27, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514.LostAngMom": {"value": 2.3549570242977578e+33, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.LockTime": {"value": -1.0000000000000000, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.BodyDsemiDtEqtide": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.BodyDeccDt": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.DOblDtEqtide": {"value": 0.0000000000000000, "unit": u.rad / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.DRotPerDtEqtide": {"value": 2.3670446957670796e-20, "rtol": 1e-4}, 
       "log.final.gl514.DRotRateDtEqtide": {"value": -2.1346978387082101e-32, "unit": 1 / u.sec ** 2, "rtol": 1e-4}, 
       "log.final.gl514.EqRotRateDiscrete": {"value": 7.7800437205928774e-07, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.EqRotPerDiscrete": {"value": 8.0760282754564993e+06, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.EqRotRateCont": {"value": 1.5164601885389766e-06, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.EqRotPerCont": {"value": 4.1433236128889602e+06, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.EqRotPer": {"value": 4.1433236128889602e+06, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.EqTidePower": {"value": 0.0000000000000000, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514.GammaRot": {"value": -1.0000000000000000, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.GammaOrb": {"value": -1.0000000000000000, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.OceanK2": {"value": 0.0100000000000000, "rtol": 1e-4}, 
       "log.final.gl514.EnvTidalQ": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.OceanTidalQ": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.TideLock": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.RotTimeEqtide": {"value": 1.1151120264596733e+26, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514.EnvK2": {"value": 0.0100000000000000, "rtol": 1e-4}, 
       "log.final.gl514.OblTimeEqtide": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514.PowerEqtide": {"value": 1.2346074026431728e+10, "unit": u.W, "rtol": 1e-4}, 
       "log.final.gl514.SurfEnFluxEqtide": {"value": 5.6927005606372488e-10, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514.Luminosity": {"value": 2.6238552211501671e+26, "unit": u.W, "rtol": 1e-4}, 
       "log.final.gl514.LXUVStellar": {"value": 2.6238552211501671e+23, "unit": u.W, "rtol": 1e-4}, 
       "log.final.gl514.Temperature": {"value": 3818.8020013755422042, "unit": u.K, "rtol": 1e-4}, 
       "log.final.gl514.LXUVFrac": {"value": 0.0010000000000000, "rtol": 1e-4}, 
       "log.final.gl514.RossbyNumber": {"value": 0.6873965303916811, "rtol": 1e-4}, 
       "log.final.gl514.DRotPerDtStellar": {"value": 2.3504479143348279e-10, "rtol": 1e-4}, 
       "log.final.gl514b.Mass": {"value": 3.1053088098643201e+25, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.gl514b.Obliquity": {"value": 0.4101523752957674, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.PrecA": {"value": 0.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.Xobl": {"value": 0.3987490699130116, "rtol": 1e-4}, 
       "log.final.gl514b.Yobl": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.Zobl": {"value": 0.9170600739556314, "rtol": 1e-4}, 
       "log.final.gl514b.Radius": {"value": 1.3394010000000000e+07, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.RadGyra": {"value": 0.5000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.RotAngMom": {"value": 1.0128197504166973e+35, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.RotKinEnergy": {"value": 3.6827164855924537e+30, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514b.RotVel": {"value": 974.0398805393790553, "unit": u.m / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.BodyType": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.RotRate": {"value": 7.2722051166109258e-05, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.RotPer": {"value": 8.6400001188466835e+04, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.Density": {"value": 3085.2071528551509800, "unit": u.kg / u.m ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.SurfEnFluxTotal": {"value": 1.4096430839727192, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.TidalQ": {"value": 12.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.ImK2": {"value": -0.0250000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.K2": {"value": 0.3000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.K2Man": {"value": 0.0100000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.Imk2Man": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.TidalQMantle": {"value": 100.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.HEcc": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.HZLimitDryRunaway": {"value": 1.0509802225887555e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.HZLimRecVenus": {"value": 9.9104353063518631e+10, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.HZLimRunaway": {"value": 1.3136272247075568e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.HZLimMoistGreenhouse": {"value": 1.3113105818748537e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.HZLimMaxGreenhouse": {"value": 2.4324580025875378e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.HZLimEarlyMars": {"value": 2.6530630489685046e+11, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.Instellation": {"value": 5866.6331509129950064, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.KEcc": {"value": 0.4500000000000184, "rtol": 1e-4}, 
       "log.final.gl514b.Eccentricity": {"value": 0.4500000000000184, "rtol": 1e-4}, 
       "log.final.gl514b.OrbEnergy": {"value": -1.6646340895689042e+34, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514b.MeanMotion": {"value": 5.1866958137285849e-07, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.OrbPeriod": {"value": 1.2114042413184749e+07, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.SemiMajorAxis": {"value": 6.3130301435403900e+10, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.CriticalSemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.COPP": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.OrbAngMom": {"value": 5.7322265617905455e+40, "unit": (u.kg * u.m ** 2) / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.LongP": {"value": 0.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.TotOrbEnergy": {"value": -1.6646340895689042e+34, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514b.OrbPotEnergy": {"value": -3.3292681791378085e+34, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514b.LostEnergy": {"value": 1.0028690689225319e+23, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514b.TidalRadius": {"value": 1.3394010000000000e+07, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.DsemiDtEqtide": {"value": 1.2351242905112557e-10, "unit": u.m / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DeccDtEqtide": {"value": 5.8376560316535877e-22, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DMeanMotionDtEqtide": {"value": -1.5221408360256703e-27, "unit": 1 / u.sec ** 2, "rtol": 1e-4}, 
       "log.final.gl514b.DOrbPerDtEqtide": {"value": 3.5551108660833406e-14, "rtol": 1e-4}, 
       "log.final.gl514b.EccTimeEqtide": {"value": 7.7085733993229202e+20, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.SemiTimeEqtide": {"value": 5.1112509016620781e+20, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DHEccDtEqtide": {"value": 0.0000000000000000, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DKEccDtEqtide": {"value": 5.8376560316535877e-22, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DXoblDtEqtide": {"value": 3.1641064275885809e-17, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DYoblDtEqtide": {"value": 0.0000000000000000, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DZoblDtEqtide": {"value": -1.3757926344613390e-17, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.LockTime": {"value": -1.0000000000000000, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.BodyDsemiDtEqtide": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.BodyDeccDt": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.DOblDtEqtide": {"value": 3.4502717078725040e-17, "unit": u.rad / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.DRotPerDtEqtide": {"value": 3.7660242031908298e-11, "rtol": 1e-4}, 
       "log.final.gl514b.DRotRateDtEqtide": {"value": -3.1698264008103682e-20, "unit": 1 / u.sec ** 2, "rtol": 1e-4}, 
       "log.final.gl514b.EqRotRateDiscrete": {"value": 7.7800437205928774e-07, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.EqRotPerDiscrete": {"value": 8.0760282754564993e+06, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.EqRotRateCont": {"value": 1.5164601885389766e-06, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.EqRotPerCont": {"value": 4.1433236128889602e+06, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.EqRotPer": {"value": 4.1433236128889602e+06, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.EqTidePower": {"value": 3.4384723793782301e+13, "unit": 1 / u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.GammaRot": {"value": -1.0000000000000000, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.GammaOrb": {"value": -1.0000000000000000, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.OceanK2": {"value": 0.0100000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.EnvTidalQ": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.OceanTidalQ": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.TideLock": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.RotTimeEqtide": {"value": 2.2941966521421430e+15, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.EnvK2": {"value": 0.0100000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.OblTimeEqtide": {"value": -1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.PowerEqtide": {"value": 3.1779002901179495e+15, "unit": u.W, "rtol": 1e-4}, 
       "log.final.gl514b.SurfEnFluxEqtide": {"value": 1.4096430839727192, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.TGlobal": {"value": 676.2684488390327715, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.AlbedoGlobal": {"value": 0.2996778493936415, "rtol": 1e-4}, 
       "log.final.gl514b.FluxInGlobal": {"value": 1046.1316028588294103, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.FluxOutGlobal": {"value": 1046.1310580735782878, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.TotIceMass": {"value": 0.0000000000000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.gl514b.TotIceFlow": {"value": 0.0000000000000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.gl514b.TotIceBalance": {"value": 0.0000000000000000, "unit": u.kg, "rtol": 1e-4}, 
       "log.final.gl514b.SkipSeas": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.AreaIceCov": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.Latitude": {"value": -1.4552620265106593, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.TempLat": {"value": 557.2252934226596608, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.AlbedoLat": {"value": 0.3675649877648788, "rtol": 1e-4}, 
       "log.final.gl514b.AnnInsol": {"value": 755.3407468728647700, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.FluxMerid": {"value": -2.0352457779997064e+16, "unit": u.Joule, "rtol": 1e-4}, 
       "log.final.gl514b.FluxIn": {"value": 486.1623348490650187, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.FluxOut": {"value": 797.3308632533586433, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.DivFlux": {"value": -311.1690732121418819, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.IceMass": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceHeight": {"value": 0.0000000000000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.DIceMassDt": {"value": 0.0000000000000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.IceFlow": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.EnergyResL": {"value": 7.9239725891966373e-11, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.EnergyResW": {"value": 6.8553163146134466e-11, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.BedrockH": {"value": 0.0000000000000000, "unit": u.m, "rtol": 1e-4}, 
       "log.final.gl514b.TempLandLat": {"value": 545.8665750487834885, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.TempWaterLat": {"value": 563.0767544031414218, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.AlbedoLandLat": {"value": 0.4335649877648788, "rtol": 1e-4}, 
       "log.final.gl514b.AlbedoWaterLat": {"value": 0.3335649877648789, "rtol": 1e-4}, 
       "log.final.gl514b.TempMinLat": {"value": 520.6587997956662548, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.TempMaxLat": {"value": 598.1732888468328611, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.Snowball": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.PlanckBAvg": {"value": 2.0899999999999990, "rtol": 1e-4}, 
       "log.final.gl514b.IceAccum": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceAblate": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.TempMaxLand": {"value": 657.9980332140321480, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.TempMaxWater": {"value": 568.6092758126177387, "unit": u.sec, "rtol": 1e-4}, 
       "log.final.gl514b.PeakInsol": {"value": 4243.8792251330132785, "unit": u.kg / u.sec ** 3, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapNorthLand": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapNorthSea": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapSouthLand": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapSouthSea": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceBeltLand": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceBeltSea": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.SnowballLand": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.SnowballSea": {"value": 0.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceFree": {"value": 1.0000000000000000, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapNorthLatLand": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapNorthLatSea": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapSouthLatLand": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.IceCapSouthLatSea": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.IceBeltNorthLatLand": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.IceBeltNorthLatSea": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.IceBeltSouthLatLand": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
       "log.final.gl514b.IceBeltSouthLatSea": {"value": 100.0000000000000000, "unit": u.rad, "rtol": 1e-4}, 
   } 
)
class Test_CPL_POISE(Benchmark): 
   pass 
