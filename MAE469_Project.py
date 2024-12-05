# MAE 469 - Introduction to Astrodynamics
# Group Members : Josh Hall, Collin Horne, Makenzie Karakula, Mike Smith, Zane Wilson
# Date : November 02, 2024 - PRESENT
# ==================================
# ---- IMPORT LIBRARIES ----
import math
import numpy as np
from datetime import datetime

from astropy import units as u

import two_body_util as util
import MAE469_ProjectLibrary as PROJ



# ---- PLANETARY ORBITAL ELEMENT INPUTS -----
# Note: a (AU) , e (n/a) , i (deg) , raan (deg) , omega (deg) , theta (deg)
# -- EARTH --
aE = 1.000000; eE = 0.01671; iE = 0.00005; raanE = -11.26064; omegaE = 114.20783; thetaE = -2.48284
# -- MARS --
aM = 1.523662; eM = 0.093412; iM = 1.85061; raanM = 49.57854; omegaM = 286.4623; thetaM = 19.41248
# -- JUPITER -- 
aJ = 5.203363; eJ = 0.048393; iJ = 1.3053; raanJ = 100.55615; omegaJ = -85.8023; thetaJ = 19.55053
# ---- GRAVITATIONAL CONSTANT W.R.T HELIOCENTRIC FRAME ----
mu = 1.0 # AU^3 / TU^2

# ---- CURRENT DATE PARAMETERS ----
Month = "Dec" # Input month as the first three letters, first capitalized.
Day = 25
Year = 2025
Hour = 8 # Uses military time 24HR system.
Min = 37
Unit = "TUsun" # Choose 'seconds' or 'TUsun' as result returned.

# ---- VARIOUS TABLE / PRINT PARAMETERS ----
raddeg = "degree" # The 'degree' or 'radian' for parameters.
table = "elements" # The 'full' prints all calculations and results, 'elements' only forwards results.
distance = "AU" # The distance unit for table printing. Mainly used for 'full' prints.
time = "TU" # The time unit for table printing. Mainly used for 'full' prints.

# ==== STEP 1 : Planetary Motion Model (r,v,theta) for December 25, 2025 at 8:37am UTC ====
# -- Calculate total time since J2000 Epoch. --
TOF = PROJ.TOFofInterest(Month,Day,Year,Hour,Min,Unit)
# -- Initial position and velocity vectors of planetary orbits from orbital elements. --
[repochE,vepochE] = PROJ.OrbitalElementsToRV(aE,eE,iE,raanE,omegaE,thetaE,mu,raddeg,table,distance,time)
[repochM,vepochM] = PROJ.OrbitalElementsToRV(aM,eM,iM,raanM,omegaM,thetaM,mu,raddeg,table,distance,time)
[repochJ,vepochJ] = PROJ.OrbitalElementsToRV(aJ,eJ,iJ,raanJ,omegaJ,thetaJ,mu,raddeg,table,distance,time)
# -- Calculate current heliocentric position and velocity vectors from epoch standards. --
[r2E,v2E] = PROJ.universalTOF_SCZ(repochE,vepochE,TOF,mu)
[r2M,v2M] = PROJ.universalTOF_SCZ(repochM,vepochM,TOF,mu)
[r2J,v2J] = PROJ.universalTOF_SCZ(repochJ,vepochJ,TOF,mu)
# -- Calculate orbital elements at current heliocentric position and velocity (find true anomaly). --
[e_scalarE2, aE, inclinationE2, TrueAnomalyE2, longAscendNodeE2, ArgOfPeriE2] = PROJ.rvToOrbitalElements(r2E,v2E,mu,table,distance,time)
[e_scalarM2, aM2, inclinationM2, TrueAnomalyM2, longAscendNodeM2, ArgOfPeriM2] = PROJ.rvToOrbitalElements(r2M,v2M,mu,table,distance,time)
[e_scalarJ2, aJ2, inclinationJ2, TrueAnomalyJ2, longAscendNodeJ2, ArgOfPeriJ2] = PROJ.rvToOrbitalElements(r2J,v2J,mu,table,distance,time)

print('\nHeliocentric Position and Velocity of Planetary Objects')
print('-------------------------------------------------------')
print('        |    Heliocentric Position Vector (r) %s    |   Heliocentric Velocity Vector (v) %s/%s    | True Anomaly (degrees)' % (distance,distance,time))
print('--------------------------------------------------------------------------------------------------------------------------')
print(' Earth  |        ',np.around(r2E,4), '        |         ',np.around(v2E,4), '         |       ', np.around(TrueAnomalyE2,4))
print('  Mars  |        ',np.around(r2M,4), '        |         ',np.around(v2M,4), '         |       ', np.around(TrueAnomalyM2,4))
print('Jupiter |        ',np.around(r2J,4), '        |         ',np.around(v2J,4), '         |       ', np.around(TrueAnomalyJ2,4),'\n')

# ==== STEP 2 : Determine Launch Date Between 2021-2030 For 190-Day Transfer from Earth to Mars ====
#PROJ.plotSynodicPeriod()
#[first_guess_departure, first_guess_arrival] = PROJ.reasonableLaunchSearch('1/1/21 12:00')
#print(first_guess_departure)
#[rDeparture,rArrival] = PROJ.launchDatePositions(first_guess_departure, first_guess_arrival)
#print(rDeparture)
#print(rArrival)

dt = (190*u.day).to(util.TU_SUN)
epoch = datetime.fromisoformat('2000-01-01 11:58:00.000')
earth_orbit = 200*u.km + 1*util.DU_EARTH
mars_orbit = 1000*u.km + 3389.5*u.km

start_date = datetime.fromisoformat('2021-01-01 12:00:00.000')
end_date = datetime.fromisoformat('2030-01-01 12:00:00.000')
#end_date = datetime.fromisoformat('2022-10-01 12:00:00.000')
step_size = 1*u.day

best_date, best_dv = PROJ.find_best_dv(start_date, end_date, step_size,
                                       repochE, vepochE, repochM, vepochM,
                                       mu, dt, earth_orbit, mars_orbit, epoch)

best_TOF = ((abs((best_date - epoch).total_seconds())*u.s).to(util.TU_SUN)).value

[r2E,v2E] = PROJ.universalTOF_SCZ(repochE,vepochE,best_TOF,mu)
[r2M,v2M] = PROJ.universalTOF_SCZ(repochM,vepochM,best_TOF,mu)   

PROJ.earth_to_mars_dv(r2E*util.AU_SUN, v2E*util.AUTU_SUN, 
                        r2M*util.AU_SUN, v2M*util.AUTU_SUN, 
                        dt, earth_orbit, mars_orbit,
                        debug=True)

[e_scalarE2, aE, inclinationE2, TrueAnomalyE2, longAscendNodeE2, ArgOfPeriE2] = PROJ.rvToOrbitalElements(r2E,v2E,mu,table,distance,time)
[e_scalarM2, aM2, inclinationM2, TrueAnomalyM2, longAscendNodeM2, ArgOfPeriM2] = PROJ.rvToOrbitalElements(r2M,v2M,mu,table,distance,time)

print('\nHeliocentric Position and Velocity of Planetary Objects at Transfer Start')
print('-------------------------------------------------------')
print('        |    Heliocentric Position Vector (r) %s    |   Heliocentric Velocity Vector (v) %s/%s    | True Anomaly (degrees)' % (distance,distance,time))
print('--------------------------------------------------------------------------------------------------------------------------')
print(' Earth  |        ',np.around(r2E,4), '        |         ',np.around(v2E,4), '         |       ', np.around(TrueAnomalyE2,4))
print('  Mars  |        ',np.around(r2M,4), '        |         ',np.around(v2M,4), '         |       ', np.around(TrueAnomalyM2,4))
print(best_date)
print(f"dV for transfer: {best_dv}")

PROJ.zane_plot(r2E, v2E, r2M, v2M)

