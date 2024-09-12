import os

import numpy as np
from astropy import units as u

from orbits import Orbit, OrbitalElements
import two_body_util as util

width = os.get_terminal_size().columns

# Question 2: Earth Radar Site ----------------------------------------------------
print("\n" + "-"*width)
print("Question 2: Earth Radar Site\n")

r = np.array([2, 0.5, 1])*util.DU_EARTH
v = np.array([0.5, 0.5, -0.5])*util.DUTU_EARTH
meu = 1*util.MEU_EARTH

orbit = Orbit(r_vector=r, v_vector=v, meu=meu)
print(orbit.orbital_elements)

# Question 3: Earth Satellite r, v vectors ----------------------------------------
print("\n" + "-"*width)
print("Question 3: Earth Satellite r, v vectors\n")
meu = 1*util.MEU_EARTH
a = 6.392*util.DU_EARTH
e = 0.488   
i = np.deg2rad(63.5)*u.rad
raan = np.deg2rad(96.4)*u.rad
omega = np.deg2rad(246)*u.rad
theta = np.deg2rad(18)*u.rad
oe = OrbitalElements(a, e, i, raan, omega, theta)
orbit = Orbit(orbital_elements=oe, meu=meu)
print(orbit.r_vector)
print(orbit.v_vector)

