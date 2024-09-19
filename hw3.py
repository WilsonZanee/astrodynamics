import os

import numpy as np
from astropy import units as u
from astropy.units import imperial as imp

from orbits import Orbit, OrbitalElements, RadarObservation

import two_body_util as util

width = os.get_terminal_size().columns

# Question 1: Radar observation to Orbital Elements ---------------------------
print("\n" + "-"*width + "\n")
print("Question 1: Radar observation to Orbital Elements\n")

#range = 3*util.DU_EARTH
#range_rate = 0.4*util.DUTU_EARTH
#az = (30*u.deg).to(u.rad)
#az_rate = 0*u.rad/util.TU_EARTH
#el = (10*u.deg).to(u.rad)
#el_rate = 0.1*u.rad/util.TU_EARTH

range = 5*util.DU_EARTH
range_rate = 0.33*util.DUTU_EARTH
az = (25*u.deg).to(u.rad)
az_rate = 0.6*u.rad/util.TU_EARTH
el = (86*u.deg).to(u.rad)
el_rate = 1*u.rad/util.TU_EARTH

ref_zulu_sid_time = 0*u.deg
lat = -106.6*u.deg
dt = util.get_sec("18:00:00")*u.s


obs = RadarObservation(range, 
                       range_rate, 
                       az, 
                       az_rate, 
                       el, 
                       el_rate, 
                       lat=lat, 
                       ref_zulu_theta=ref_zulu_sid_time, 
                       dt=dt)

r, v = obs.get_vectors_in_geocentric_IJK()
print(f"r{r}")
print(f"v{v}")


#print(orbit.orbital_elements)

# Question 2:  ------------------------------------
print("\n" + "-"*width + "\n")
print("\n")

# Question 3:  -------------------------------------------------
print("\n" + "-"*width + "\n")
print("\n")


