import os

import numpy as np
from astropy import units as u
from astropy.units import imperial as imp

from orbits import Orbit, OrbitalElements
import two_body_util as util

width = os.get_terminal_size().columns

# Question 2: Earth Radar Site ------------------------------------------------
print("\n" + "-"*width + "\n")
print("Question 2: Earth Radar Site\n")

r = np.array([2, 0.5, 1])*util.DU_EARTH
v = np.array([0.5, 0.5, -0.5])*util.DUTU_EARTH
meu = 1*util.MEU_EARTH

orbit = Orbit(r_vector=r, v_vector=v, meu=meu)
print(orbit.orbital_elements)

# Question 3: Earth Satellite r, v vectors ------------------------------------
print("\n" + "-"*width + "\n")
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
print(f"r: {orbit.r_vector}")
print(f"v: {orbit.v_vector.to(util.DUTU_EARTH, equivalencies=util.du_tu)}")

# Question 4: Unit Conversion -------------------------------------------------
print("\n" + "-"*width + "\n")
print("Question 4: Unit Conversion\n")

a_mars = (1.524*util.AU_SUN).to(u.km)
v_mars = (0.795*util.AUTU_SUN).to(u.km/u.s)
print(f"Mars Orbit Semi-Major Axis: {a_mars}")
print(f"Mars Orbit Velocity: {v_mars}")

r_voyager = (14154247350*imp.mile).to(util.AU_SUN)
v_voyager = (38027*imp.mile/u.hour).to(util.AUTU_SUN)
print(f"Voyager radius: {r_voyager}")
print(f"Voyager velocity: {v_voyager}")


