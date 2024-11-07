
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
        dt = (abs((time - self.epoch).total_seconds())*u.s).to(util.TU_SUN)
        r2, v2 = util.time_of_flight_universal_var(r_epoch, 
                                                v_epoch, 
                                                dt, 
                                                self.meu)
        orbit = Orbit(r_vector=r2, v_vector=v2.to(util.AUTU_SUN), meu=self.meu)

        return orbit
    
    def get_rv_theta(self, time):
        orbit = self.get_orbit_at_time(time)
        r = orbit.r_vector
        v = orbit.v_vector
        theta = orbit.orbital_elements.theta

        return (r, v, theta.to(u.deg))
    
class Planet_Proj:
    def __init__(self, a, e, i, raan, omega, theta, epoch, meu):
        self.a = a
        self.e = e
        self.i = i
        self.raan = raan
        self.omega = omega
        self.theta = theta
        self.epoch = epoch
        [repochE,vepochE] = PROJ.OrbitalElementsToRV(a, e, i, raan, omega, 
                                                     theta,meu,raddeg,table,
                                                     distance,time)
        self.r_epoch = repochE
        self.v_epoch = vepochE
        self.mu = meu
        
    def get_orbit_at_time(self, time):
        dt = (abs((time - self.epoch).total_seconds()))
        [r2,v2] = PROJ.universalTOF_SCZ(self.r_epoch,self.v_epoch,dt,self.mu)
        [e_scalar, a, inclination, TrueAnomaly, longAscendNode, ArgOfPeri] = PROJ.rvToOrbitalElements(r2,v2,self.mu,table,distance,time)

        return (r2, v2, ArgOfPeri)
