
from datetime import datetime

import numpy as np

from astropy import units as u

import MAE469_ProjectLibrary as PROJ
from orbits import Orbit
import two_body_util as util


raddeg = "degree" # The 'degree' or 'radian' for parameters.
table = "elements" # The 'full' prints all calculations and results, 'elements' only forwards results.
distance = "AU" # The distance unit for table printing. Mainly used for 'full' prints.
time = "TU" # The time unit for table printing. Mainly used for 'full' prints.

class Planet:
    def __init__(self, orbital_elements, epoch, meu):
        self.orbit_at_epoch = Orbit(orbital_elements=orbital_elements, meu=meu)
        self.epoch = epoch
        self.meu = meu

    def get_orbit_at_time(self, time):
        r_epoch = (self.orbit_at_epoch.r_vector).flatten()
        v_epoch = (self.orbit_at_epoch.v_vector).flatten().to(util.AUTU_SUN)
        print(r_epoch, v_epoch)
        dt = (abs((time - self.epoch).total_seconds())*u.s).to(util.TU_SUN)
        r2, v2 = util.time_of_flight_universal_var(r_epoch, 
                                                   v_epoch, 
                                                   dt, 
                                                   self.meu)
        print (r2, v2)
        orbit = Orbit(r_vector=r2, v_vector=v2.to(util.AUTU_SUN), meu=self.meu)
        return orbit
    
    def get_rv_theta(self, time):
        orbit = self.get_orbit_at_time(time)
        r = orbit.r_vector
        v = orbit.v_vector
        theta = orbit.orbital_elements.theta

        return (r, v, theta.to(u.deg))