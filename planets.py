
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
        r_epoch = (self.orbit_at_epoch.r_vector).value
        v_epoch = (self.orbit_at_epoch.v_vector).value
        dt = (abs((time - self.epoch).total_seconds())*u.s).to(util.TU_SUN)
        print(dt)
        [r2E,v2E] = PROJ.universalTOF_SCZ(r_epoch, 
                                          v_epoch, 
                                          dt.value, 
                                          self.meu.value)
        orbit = Orbit(r_vector=r2E*util.AU_SUN, 
                      v_vector=v2E*util.AUTU_SUN, 
                      meu=self.meu)
        return orbit