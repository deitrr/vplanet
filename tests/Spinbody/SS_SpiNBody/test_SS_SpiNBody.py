from benchmark import Benchmark, benchmark 
import astropy.units as u 
import pytest 
 
@benchmark( 
   { 
       "log.initial.system.Age": {"value": 0.0000000000000000, "unit": u.sec}, 
       "log.initial.system.Time": {"value": 0.0000000000000000, "unit": u.sec}, 
       "log.initial.system.TotAngMom": {"value": 4.9405684954906725e+40, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.system.TotEnergy": {"value": -2.4289670619261560e+40, "unit": u.Joule}, 
       "log.initial.system.PotEnergy": {"value": -2.4824390943032254e+40, "unit": u.Joule}, 
       "log.initial.system.KinEnergy": {"value": 5.3472690917973792e+38, "unit": u.Joule}, 
       "log.initial.system.DeltaTime": {"value": 0.0000000000000000, "unit": u.sec}, 
       "log.initial.Star.Mass": {"value": 1.9884154399559999e+30, "unit": u.kg}, 
       "log.initial.Star.Radius": {"value": 6.3781000000000000e+09, "unit": u.m}, 
       "log.initial.Star.RadGyra": {"value": 0.5000000000000000}, 
       "log.initial.Star.BodyType": {"value": 0.0000000000000000}, 
       "log.initial.Star.Density": {"value": 1.8295511047660136, "unit": u.kg / u.m ** 3}, 
       "log.initial.Star.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.Star.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.initial.Star.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.initial.Star.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Star.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Star.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.initial.Star.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.Star.Eccentricity": {"value": -1.0000000000000000}, 
       "log.initial.Star.MeanMotion": {"value": -1.0000000000000000, "unit": 1 / u.sec}, 
       "log.initial.Star.OrbPeriod": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.initial.Star.SemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.Star.COPP": {"value": 0.0000000000000000}, 
       "log.initial.Star.OrbAngMom": {"value": 1.1839688346281812e+35, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.Star.ArgP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Star.MeanAnomaly": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Star.Inc": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Star.LongA": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Star.MeanLongitude": {"value": -1.0000000000000000, "unit": u.rad}, 
       "log.initial.Star.LongP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Star.TotOrbEnergy": {"value": -6.1982871601954894e+33, "unit": u.Joule}, 
       "log.initial.Star.OrbPotEnergy": {"value": -1.0000000000000000, "unit": u.Joule}, 
       "log.initial.Star.PositionXSpiNBody": {"value": 2.7915940084973827e+05}, 
       "log.initial.Star.PositionYSpiNBody": {"value": -4.1014856265306339e+05}, 
       "log.initial.Star.PositionZSpiNBody": {"value": -3.7106311591366364e+04}, 
       "log.initial.Star.VelXSpiNBody": {"value": 0.0795471921909874}, 
       "log.initial.Star.VelYSpiNBody": {"value": 0.0958149953307798}, 
       "log.initial.Star.VelZSpiNBody": {"value": 0.0039116626797918}, 
       "log.initial.Star.SpiNBodyInc": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Star.SpiNBodyLongA": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Mercury.Mass": {"value": 3.3026188579999998e+23, "unit": u.kg}, 
       "log.initial.Mercury.Radius": {"value": 2.4428122999999998e+06, "unit": u.m}, 
       "log.initial.Mercury.RadGyra": {"value": 0.5000000000000000}, 
       "log.initial.Mercury.BodyType": {"value": 0.0000000000000000}, 
       "log.initial.Mercury.Density": {"value": 5408.7829368350048753, "unit": u.kg / u.m ** 3}, 
       "log.initial.Mercury.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.Mercury.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.initial.Mercury.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.initial.Mercury.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Mercury.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Mercury.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.initial.Mercury.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.Mercury.Eccentricity": {"value": 0.2056306900000000}, 
       "log.initial.Mercury.MeanMotion": {"value": 8.2667487086227897e-07, "unit": 1 / u.sec}, 
       "log.initial.Mercury.OrbPeriod": {"value": 7.6005519565700479e+06, "unit": u.sec}, 
       "log.initial.Mercury.SemiMajorAxis": {"value": 5.7909175678248352e+10, "unit": u.m}, 
       "log.initial.Mercury.COPP": {"value": 0.0001191111427544}, 
       "log.initial.Mercury.OrbAngMom": {"value": 8.9600164476519010e+38, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.Mercury.ArgP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Mercury.MeanAnomaly": {"value": 3.0507376194278546, "unit": u.rad}, 
       "log.initial.Mercury.Inc": {"value": 0.1222580451706808, "unit": u.rad}, 
       "log.initial.Mercury.LongA": {"value": 0.8435467744632575, "unit": u.rad}, 
       "log.initial.Mercury.MeanLongitude": {"value": 4.4026076987955776, "unit": u.rad}, 
       "log.initial.Mercury.LongP": {"value": 1.3518700793677227, "unit": u.rad}, 
       "log.initial.Mercury.TotOrbEnergy": {"value": -6.1982871601954894e+33, "unit": u.Joule}, 
       "log.initial.Mercury.OrbPotEnergy": {"value": -1.6734887233371770e+27, "unit": u.Joule}, 
       "log.initial.Mercury.PositionXSpiNBody": {"value": -1.9460729480525391e+10}, 
       "log.initial.Mercury.PositionYSpiNBody": {"value": -6.6914040443213783e+10}, 
       "log.initial.Mercury.PositionZSpiNBody": {"value": -3.6797570861737509e+09}, 
       "log.initial.Mercury.VelXSpiNBody": {"value": 3.6995108268405180e+04}, 
       "log.initial.Mercury.VelYSpiNBody": {"value": -1.1164070038032824e+04}, 
       "log.initial.Mercury.VelZSpiNBody": {"value": -4307.5569771203754499}, 
       "log.initial.Mercury.SpiNBodyInc": {"value": 0.1222580451706808, "unit": u.rad}, 
       "log.initial.Mercury.SpiNBodyLongA": {"value": 0.8435467744632575, "unit": u.rad}, 
       "log.initial.Venus.Mass": {"value": 4.8673315899999993e+24, "unit": u.kg}, 
       "log.initial.Venus.Radius": {"value": 6.0528168999999994e+06, "unit": u.m}, 
       "log.initial.Venus.RadGyra": {"value": 0.5000000000000000}, 
       "log.initial.Venus.BodyType": {"value": 0.0000000000000000}, 
       "log.initial.Venus.Density": {"value": 5239.9810321605646095, "unit": u.kg / u.m ** 3}, 
       "log.initial.Venus.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.Venus.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.initial.Venus.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.initial.Venus.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Venus.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Venus.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.initial.Venus.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.Venus.Eccentricity": {"value": 0.0067732300000000}, 
       "log.initial.Venus.MeanMotion": {"value": 3.2363961741823590e-07, "unit": 1 / u.sec}, 
       "log.initial.Venus.OrbPeriod": {"value": 1.9414141437016640e+07, "unit": u.sec}, 
       "log.initial.Venus.SemiMajorAxis": {"value": 1.0820892551319370e+11, "unit": u.m}, 
       "log.initial.Venus.COPP": {"value": 0.0000000000000000}, 
       "log.initial.Venus.OrbAngMom": {"value": 1.8444487737839402e+40, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.Venus.ArgP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Venus.MeanAnomaly": {"value": 0.8804618844052841, "unit": u.rad}, 
       "log.initial.Venus.Inc": {"value": 0.0592488666486833, "unit": u.rad}, 
       "log.initial.Venus.LongA": {"value": 1.3383305131628385, "unit": u.rad}, 
       "log.initial.Venus.MeanLongitude": {"value": 3.1761454602995198, "unit": u.rad}, 
       "log.initial.Venus.LongP": {"value": 2.2956835758942358, "unit": u.rad}, 
       "log.initial.Venus.TotOrbEnergy": {"value": -6.1982871601954894e+33, "unit": u.Joule}, 
       "log.initial.Venus.OrbPotEnergy": {"value": -1.2078693330250667e+28, "unit": u.Joule}, 
       "log.initial.Venus.PositionXSpiNBody": {"value": -1.0745862317120734e+11}, 
       "log.initial.Venus.PositionYSpiNBody": {"value": -4.8950065818498735e+09}, 
       "log.initial.Venus.PositionZSpiNBody": {"value": 6.1358940673992805e+09}, 
       "log.initial.Venus.VelXSpiNBody": {"value": 1383.6855202680806087}, 
       "log.initial.Venus.VelYSpiNBody": {"value": -3.5139320872477445e+04}, 
       "log.initial.Venus.VelZSpiNBody": {"value": -560.0632256368625121}, 
       "log.initial.Venus.SpiNBodyInc": {"value": 0.0592488666486833, "unit": u.rad}, 
       "log.initial.Venus.SpiNBodyLongA": {"value": 1.3383305131628385, "unit": u.rad}, 
       "log.initial.Earth.Mass": {"value": 5.9721859999999998e+24, "unit": u.kg}, 
       "log.initial.Earth.Radius": {"value": 6.3781000000000000e+06, "unit": u.m}, 
       "log.initial.Earth.RadGyra": {"value": 0.5000000000000000}, 
       "log.initial.Earth.BodyType": {"value": 0.0000000000000000}, 
       "log.initial.Earth.Density": {"value": 5495.0385490920862139, "unit": u.kg / u.m ** 3}, 
       "log.initial.Earth.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.Earth.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.initial.Earth.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.initial.Earth.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Earth.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Earth.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.initial.Earth.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.Earth.Eccentricity": {"value": 0.0167102200000000}, 
       "log.initial.Earth.MeanMotion": {"value": 1.9909861410536333e-07, "unit": 1 / u.sec}, 
       "log.initial.Earth.OrbPeriod": {"value": 3.1558156923455600e+07, "unit": u.sec}, 
       "log.initial.Earth.SemiMajorAxis": {"value": 1.4959788715576578e+11, "unit": u.m}, 
       "log.initial.Earth.COPP": {"value": 0.0064781433321986}, 
       "log.initial.Earth.OrbAngMom": {"value": 2.6606583190350611e+40, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.Earth.ArgP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Earth.MeanAnomaly": {"value": 6.2398515742010225, "unit": u.rad}, 
       "log.initial.Earth.Inc": {"value": 0.0593411945661111, "unit": u.rad}, 
       "log.initial.Earth.LongA": {"value": 6.0866500631238427, "unit": u.rad}, 
       "log.initial.Earth.MeanLongitude": {"value": 8.0366189953258491, "unit": u.rad}, 
       "log.initial.Earth.LongP": {"value": 1.7967674211248261, "unit": u.rad}, 
       "log.initial.Earth.TotOrbEnergy": {"value": -6.1982871601954894e+33, "unit": u.Joule}, 
       "log.initial.Earth.OrbPotEnergy": {"value": -9.2006090724902290e+26, "unit": u.Joule}, 
       "log.initial.Earth.PositionXSpiNBody": {"value": -2.6549739061807003e+10}, 
       "log.initial.Earth.PositionYSpiNBody": {"value": 1.4445683679635303e+11}, 
       "log.initial.Earth.PositionZSpiNBody": {"value": 8.1090710842026081e+09}, 
       "log.initial.Earth.VelXSpiNBody": {"value": -2.9782610455354614e+04}, 
       "log.initial.Earth.VelYSpiNBody": {"value": -5459.1654954847208501}, 
       "log.initial.Earth.VelZSpiNBody": {"value": -663.6106761971726655}, 
       "log.initial.Earth.SpiNBodyInc": {"value": 0.0593411945661111, "unit": u.rad}, 
       "log.initial.Earth.SpiNBodyLongA": {"value": 6.0866500631238427, "unit": u.rad}, 
       "log.initial.Mars.Mass": {"value": 6.3902390200000002e+23, "unit": u.kg}, 
       "log.initial.Mars.Radius": {"value": 3.3931492000000002e+06, "unit": u.m}, 
       "log.initial.Mars.RadGyra": {"value": 0.5000000000000000}, 
       "log.initial.Mars.BodyType": {"value": 0.0000000000000000}, 
       "log.initial.Mars.Density": {"value": 3904.9872862933507349, "unit": u.kg / u.m ** 3}, 
       "log.initial.Mars.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.initial.Mars.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.initial.Mars.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.initial.Mars.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Mars.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.initial.Mars.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.initial.Mars.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.initial.Mars.Eccentricity": {"value": 0.0934123300000000}, 
       "log.initial.Mars.MeanMotion": {"value": 1.0586067017457037e-07, "unit": 1 / u.sec}, 
       "log.initial.Mars.OrbPeriod": {"value": 5.9353349046612397e+07, "unit": u.sec}, 
       "log.initial.Mars.SemiMajorAxis": {"value": 2.2793663724184332e+11, "unit": u.m}, 
       "log.initial.Mars.COPP": {"value": -0.0161452572136818}, 
       "log.initial.Mars.OrbAngMom": {"value": 3.4992832768748888e+39, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.initial.Mars.ArgP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.initial.Mars.MeanAnomaly": {"value": 0.3388116919678649, "unit": u.rad}, 
       "log.initial.Mars.Inc": {"value": 0.0322992376694091, "unit": u.rad}, 
       "log.initial.Mars.LongA": {"value": 0.8653087613069771, "unit": u.rad}, 
       "log.initial.Mars.MeanLongitude": {"value": 6.2038307709676923, "unit": u.rad}, 
       "log.initial.Mars.LongP": {"value": 5.8650190789998273, "unit": u.rad}, 
       "log.initial.Mars.TotOrbEnergy": {"value": -6.1982871601954894e+33, "unit": u.Joule}, 
       "log.initial.Mars.OrbPotEnergy": {"value": 0.0000000000000000, "unit": u.Joule}, 
       "log.initial.Mars.PositionXSpiNBody": {"value": 2.0803448047941187e+11}, 
       "log.initial.Mars.PositionYSpiNBody": {"value": -1.9601556499889069e+09}, 
       "log.initial.Mars.PositionZSpiNBody": {"value": -5.1582818987714539e+09}, 
       "log.initial.Mars.VelXSpiNBody": {"value": 1160.3534380638780021}, 
       "log.initial.Mars.VelYSpiNBody": {"value": 2.6297809759394786e+04}, 
       "log.initial.Mars.VelZSpiNBody": {"value": 522.4080685981735996}, 
       "log.initial.Mars.SpiNBodyInc": {"value": 0.0322992376694091, "unit": u.rad}, 
       "log.initial.Mars.SpiNBodyLongA": {"value": 0.8653087613069771, "unit": u.rad}, 
       "log.final.system.Age": {"value": 1.0000000000000000e+09, "unit": u.sec}, 
       "log.final.system.Time": {"value": 1.0000000000000000e+09, "unit": u.sec}, 
       "log.final.system.TotAngMom": {"value": 4.9405680894478570e+40, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.final.system.TotEnergy": {"value": -2.4289670619264795e+40, "unit": u.Joule}, 
       "log.final.system.PotEnergy": {"value": -2.4824390943032254e+40, "unit": u.Joule}, 
       "log.final.system.KinEnergy": {"value": 5.3472690917973792e+38, "unit": u.Joule}, 
       "log.final.system.DeltaTime": {"value": 5.9826503140891415e+04, "unit": u.sec}, 
       "log.final.Star.Mass": {"value": 1.9884154399559999e+30, "unit": u.kg}, 
       "log.final.Star.Radius": {"value": 6.3781000000000000e+09, "unit": u.m}, 
       "log.final.Star.RadGyra": {"value": 0.5000000000000000}, 
       "log.final.Star.BodyType": {"value": 0.0000000000000000}, 
       "log.final.Star.Density": {"value": 1.8295511047660136, "unit": u.kg / u.m ** 3}, 
       "log.final.Star.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.final.Star.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.final.Star.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.final.Star.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Star.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Star.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.final.Star.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.final.Star.Eccentricity": {"value": -1.0000000000000000}, 
       "log.final.Star.MeanMotion": {"value": -1.0000000000000000, "unit": 1 / u.sec}, 
       "log.final.Star.OrbPeriod": {"value": -1.0000000000000000, "unit": u.sec}, 
       "log.final.Star.SemiMajorAxis": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.final.Star.COPP": {"value": 0.0000000000000000}, 
       "log.final.Star.OrbAngMom": {"value": 2.5470730912583657e+35, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.final.Star.ArgP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.final.Star.MeanAnomaly": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.final.Star.Inc": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.final.Star.LongA": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.final.Star.MeanLongitude": {"value": -1.0000000000000000, "unit": u.rad}, 
       "log.final.Star.LongP": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.final.Star.TotOrbEnergy": {"value": -6.1982903940958106e+33, "unit": u.Joule}, 
       "log.final.Star.OrbPotEnergy": {"value": -1.0000000000000000, "unit": u.Joule}, 
       "log.final.Star.PositionXSpiNBody": {"value": -7.2996044877414301e+05}, 
       "log.final.Star.PositionYSpiNBody": {"value": 1.4037087234044124e+05}, 
       "log.final.Star.PositionZSpiNBody": {"value": 1.7226497595673525e+04}, 
       "log.final.Star.VelXSpiNBody": {"value": -0.0106736151355116}, 
       "log.final.Star.VelYSpiNBody": {"value": -0.1732430416842315}, 
       "log.final.Star.VelZSpiNBody": {"value": -0.0074053647089155}, 
       "log.final.Star.SpiNBodyInc": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.final.Star.SpiNBodyLongA": {"value": 0.0000000000000000, "unit": u.rad}, 
       "log.final.Mercury.Mass": {"value": 3.3026188579999998e+23, "unit": u.kg}, 
       "log.final.Mercury.Radius": {"value": 2.4428122999999998e+06, "unit": u.m}, 
       "log.final.Mercury.RadGyra": {"value": 0.5000000000000000}, 
       "log.final.Mercury.BodyType": {"value": 0.0000000000000000}, 
       "log.final.Mercury.Density": {"value": 5408.7829368350048753, "unit": u.kg / u.m ** 3}, 
       "log.final.Mercury.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.final.Mercury.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.final.Mercury.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.final.Mercury.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Mercury.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Mercury.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.final.Mercury.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.final.Mercury.Eccentricity": {"value": 0.2056310160605740}, 
       "log.final.Mercury.MeanMotion": {"value": 8.2668331941003352e-07, "unit": 1 / u.sec}, 
       "log.final.Mercury.OrbPeriod": {"value": 7.6004742803612044e+06, "unit": u.sec}, 
       "log.final.Mercury.SemiMajorAxis": {"value": 5.7908781130333809e+10, "unit": u.m}, 
       "log.final.Mercury.COPP": {"value": 0.0001191262369582}, 
       "log.final.Mercury.OrbAngMom": {"value": 8.9600009443842632e+38, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.final.Mercury.ArgP": {"value": 0.5092924420822773, "unit": u.rad}, 
       "log.final.Mercury.MeanAnomaly": {"value": 0.3482401512005951, "unit": u.rad}, 
       "log.final.Mercury.Inc": {"value": 0.1222430915966452, "unit": u.rad}, 
       "log.final.Mercury.LongA": {"value": 0.8431407868342485, "unit": u.rad}, 
       "log.final.Mercury.MeanLongitude": {"value": 1.7006733801171208, "unit": u.rad}, 
       "log.final.Mercury.LongP": {"value": 1.3524332289165257, "unit": u.rad}, 
       "log.final.Mercury.TotOrbEnergy": {"value": -6.1982903940958106e+33, "unit": u.Joule}, 
       "log.final.Mercury.OrbPotEnergy": {"value": -1.6243570020986443e+27, "unit": u.Joule}, 
       "log.final.Mercury.PositionXSpiNBody": {"value": -1.4284864976405479e+10}, 
       "log.final.Mercury.PositionYSpiNBody": {"value": 4.4617202242102600e+10}, 
       "log.final.Mercury.PositionZSpiNBody": {"value": 4.9562860799583359e+09}, 
       "log.final.Mercury.VelXSpiNBody": {"value": -5.6172594401872055e+04}, 
       "log.final.Mercury.VelYSpiNBody": {"value": -1.3055052637710956e+04}, 
       "log.final.Mercury.VelZSpiNBody": {"value": 4086.5414331345400569}, 
       "log.final.Mercury.SpiNBodyInc": {"value": 0.1222430915966452, "unit": u.rad}, 
       "log.final.Mercury.SpiNBodyLongA": {"value": 0.8431407868342485, "unit": u.rad}, 
       "log.final.Venus.Mass": {"value": 4.8673315899999993e+24, "unit": u.kg}, 
       "log.final.Venus.Radius": {"value": 6.0528168999999994e+06, "unit": u.m}, 
       "log.final.Venus.RadGyra": {"value": 0.5000000000000000}, 
       "log.final.Venus.BodyType": {"value": 0.0000000000000000}, 
       "log.final.Venus.Density": {"value": 5239.9810321605646095, "unit": u.kg / u.m ** 3}, 
       "log.final.Venus.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.final.Venus.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.final.Venus.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.final.Venus.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Venus.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Venus.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.final.Venus.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.final.Venus.Eccentricity": {"value": 0.0067654968650054}, 
       "log.final.Venus.MeanMotion": {"value": 3.2363398042051070e-07, "unit": 1 / u.sec}, 
       "log.final.Venus.OrbPeriod": {"value": 1.9414479588996157e+07, "unit": u.sec}, 
       "log.final.Venus.SemiMajorAxis": {"value": 1.0821018201844994e+11, "unit": u.m}, 
       "log.final.Venus.COPP": {"value": 0.0000000000000000}, 
       "log.final.Venus.OrbAngMom": {"value": 1.8444477639687901e+40, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.final.Venus.ArgP": {"value": 0.9555745530164241, "unit": u.rad}, 
       "log.final.Venus.MeanAnomaly": {"value": 4.0769944731867680, "unit": u.rad}, 
       "log.final.Venus.Inc": {"value": 0.0593149665160410, "unit": u.rad}, 
       "log.final.Venus.LongA": {"value": 1.3373158524982192, "unit": u.rad}, 
       "log.final.Venus.MeanLongitude": {"value": 6.3698848787014111, "unit": u.rad}, 
       "log.final.Venus.LongP": {"value": 2.2928904055146431, "unit": u.rad}, 
       "log.final.Venus.TotOrbEnergy": {"value": -6.1982903940958106e+33, "unit": u.Joule}, 
       "log.final.Venus.OrbPotEnergy": {"value": -3.4289361518136559e+28, "unit": u.Joule}, 
       "log.final.Venus.PositionXSpiNBody": {"value": 1.0815757144997708e+11}, 
       "log.final.Venus.PositionYSpiNBody": {"value": 8.2767987273171253e+09}, 
       "log.final.Venus.PositionZSpiNBody": {"value": -6.1349323261908092e+09}, 
       "log.final.Venus.VelXSpiNBody": {"value": -2813.6433684762664598}, 
       "log.final.Venus.VelYSpiNBody": {"value": 3.4759390241119749e+04}, 
       "log.final.Venus.VelZSpiNBody": {"value": 640.1253385111917851}, 
       "log.final.Venus.SpiNBodyInc": {"value": 0.0593149665160410, "unit": u.rad}, 
       "log.final.Venus.SpiNBodyLongA": {"value": 1.3373158524982192, "unit": u.rad}, 
       "log.final.Earth.Mass": {"value": 5.9721859999999998e+24, "unit": u.kg}, 
       "log.final.Earth.Radius": {"value": 6.3781000000000000e+06, "unit": u.m}, 
       "log.final.Earth.RadGyra": {"value": 0.5000000000000000}, 
       "log.final.Earth.BodyType": {"value": 0.0000000000000000}, 
       "log.final.Earth.Density": {"value": 5495.0385490920862139, "unit": u.kg / u.m ** 3}, 
       "log.final.Earth.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.final.Earth.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.final.Earth.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.final.Earth.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Earth.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Earth.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.final.Earth.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.final.Earth.Eccentricity": {"value": 0.0167040607934187}, 
       "log.final.Earth.MeanMotion": {"value": 1.9909898045534393e-07, "unit": 1 / u.sec}, 
       "log.final.Earth.OrbPeriod": {"value": 3.1558098855201557e+07, "unit": u.sec}, 
       "log.final.Earth.SemiMajorAxis": {"value": 1.4959770364503244e+11, "unit": u.m}, 
       "log.final.Earth.COPP": {"value": 0.0064749041550392}, 
       "log.final.Earth.OrbAngMom": {"value": 2.6606451888126111e+40, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.final.Earth.ArgP": {"value": 1.9946317274685996, "unit": u.rad}, 
       "log.final.Earth.MeanAnomaly": {"value": 4.2793278732695272, "unit": u.rad}, 
       "log.final.Earth.Inc": {"value": 0.0592958228450282, "unit": u.rad}, 
       "log.final.Earth.LongA": {"value": 6.0858921788479794, "unit": u.rad}, 
       "log.final.Earth.MeanLongitude": {"value": 6.0766664724065196, "unit": u.rad}, 
       "log.final.Earth.LongP": {"value": 1.7973385991369923, "unit": u.rad}, 
       "log.final.Earth.TotOrbEnergy": {"value": -6.1982903940958106e+33, "unit": u.Joule}, 
       "log.final.Earth.OrbPotEnergy": {"value": -1.5075691558866053e+27, "unit": u.Joule}, 
       "log.final.Earth.PositionXSpiNBody": {"value": 1.4648471525397711e+11}, 
       "log.final.Earth.PositionYSpiNBody": {"value": -3.5305272423862312e+10}, 
       "log.final.Earth.PositionZSpiNBody": {"value": -3.5065836080093104e+08}, 
       "log.final.Earth.VelXSpiNBody": {"value": 6486.6780379433712369}, 
       "log.final.Earth.VelYSpiNBody": {"value": 2.8796295693837215e+04}, 
       "log.final.Earth.VelZSpiNBody": {"value": 1751.8263488514717210}, 
       "log.final.Earth.SpiNBodyInc": {"value": 0.0592958228450282, "unit": u.rad}, 
       "log.final.Earth.SpiNBodyLongA": {"value": 6.0858921788479794, "unit": u.rad}, 
       "log.final.Mars.Mass": {"value": 6.3902390200000002e+23, "unit": u.kg}, 
       "log.final.Mars.Radius": {"value": 3.3931492000000002e+06, "unit": u.m}, 
       "log.final.Mars.RadGyra": {"value": 0.5000000000000000}, 
       "log.final.Mars.BodyType": {"value": 0.0000000000000000}, 
       "log.final.Mars.Density": {"value": 3904.9872862933507349, "unit": u.kg / u.m ** 3}, 
       "log.final.Mars.HZLimitDryRunaway": {"value": -1.0000000000000000, "unit": u.m}, 
       "log.final.Mars.HZLimRecVenus": {"value": -1.0000000000000000}, 
       "log.final.Mars.HZLimRunaway": {"value": -1.0000000000000000}, 
       "log.final.Mars.HZLimMoistGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Mars.HZLimMaxGreenhouse": {"value": -1.0000000000000000}, 
       "log.final.Mars.HZLimEarlyMars": {"value": -1.0000000000000000}, 
       "log.final.Mars.Instellation": {"value": -1.0000000000000000, "unit": u.kg / u.sec ** 3}, 
       "log.final.Mars.Eccentricity": {"value": 0.0934396712499110}, 
       "log.final.Mars.MeanMotion": {"value": 1.0585624532112576e-07, "unit": 1 / u.sec}, 
       "log.final.Mars.OrbPeriod": {"value": 5.9355830051584579e+07, "unit": u.sec}, 
       "log.final.Mars.SemiMajorAxis": {"value": 2.2794298912144470e+11, "unit": u.m}, 
       "log.final.Mars.COPP": {"value": -0.0161366204391491}, 
       "log.final.Mars.OrbAngMom": {"value": 3.4992868389766775e+39, "unit": (u.kg * u.m ** 2) / u.sec}, 
       "log.final.Mars.ArgP": {"value": 5.0000685133399099, "unit": u.rad}, 
       "log.final.Mars.MeanAnomaly": {"value": 5.6682441835193353, "unit": u.rad}, 
       "log.final.Mars.Inc": {"value": 0.0323141352804722, "unit": u.rad}, 
       "log.final.Mars.LongA": {"value": 0.8653182087158257, "unit": u.rad}, 
       "log.final.Mars.MeanLongitude": {"value": 11.5336309055750696, "unit": u.rad}, 
       "log.final.Mars.LongP": {"value": 5.8653867220557352, "unit": u.rad}, 
       "log.final.Mars.TotOrbEnergy": {"value": -6.1982903940958106e+33, "unit": u.Joule}, 
       "log.final.Mars.OrbPotEnergy": {"value": 0.0000000000000000, "unit": u.Joule}, 
       "log.final.Mars.PositionXSpiNBody": {"value": 8.5927370326019836e+10}, 
       "log.final.Mars.PositionYSpiNBody": {"value": -1.9293056341077945e+11}, 
       "log.final.Mars.PositionZSpiNBody": {"value": -6.1583914211678486e+09}, 
       "log.final.Mars.VelXSpiNBody": {"value": 2.3051619090617016e+04}, 
       "log.final.Mars.VelYSpiNBody": {"value": 1.1937602122560538e+04}, 
       "log.final.Mars.VelZSpiNBody": {"value": -317.0811459425559633}, 
       "log.final.Mars.SpiNBodyInc": {"value": 0.0323141352804722, "unit": u.rad}, 
       "log.final.Mars.SpiNBodyLongA": {"value": 0.8653182087158257, "unit": u.rad}, 
   } 
)
class Test_SS_SpiNBody(Benchmark): 
   pass 
